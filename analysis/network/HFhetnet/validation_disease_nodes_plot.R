## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-30
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
##   Validation script of disease prediction
## ---------------------------
source("analysis/utils/utils.R")
source("analysis/utils/utils_hetnet.R")

library(tidyverse)

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")
edge.list = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds")
auc_res= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/aurocs_multilayer_disease_pred_heart_genenet2022.rds")

Phe_dic= data %>% distinct(PheCode, Phenotype, category) %>% drop_na

# plot precision recall ---------------------------------------------------

df1 = map(auc_res, function(x){
  x[[1]]$pr$auc.integral
})%>% unlist()%>% enframe(name = "PheCode")%>% left_join(Phe_dic)

p.pr = df1 %>% ggplot(., aes(x= category, y= value))+
  geom_boxplot()+
  geom_jitter(alpha= 0.2)+
  coord_flip()+
  ggtitle("Precision Recall  AUC")


uff= df %>% filter(value >0.99)%>% pull(PheCode)

gd= process_gd(weight_cutoff = 0.29)
gd%>% arrange(nodeB)

x= gd%>% filter(nodeB %in% uff)

x%>% arrange(nodeB)


# full AUROC --------------------------------------------------------------

df2 = map(auc_res, function(x){
  x[[1]]$AUROC
})%>% unlist()%>% enframe(name = "PheCode")%>% left_join(Phe_dic)

p.fullAUROC = df2 %>% ggplot(., aes(x= category, y= value))+
  geom_boxplot()+
  geom_jitter(alpha= 0.2)+
  coord_flip()+
  ggtitle(paste0("full AUROC"))

df.cor= gd %>% filter(nodeB %in% disease_to_predict) %>% group_by(nodeB)%>%
  count() %>% rename(PheCode= nodeB) %>% left_join(df)

ggplot(df.cor, aes(x= value, y= n))+
  geom_point()
cor.test(df.cor$n, df.cor$value,method = "kendall")


# median rank -------------------------------------------------------------

df3 = map(auc_res, function(x){
  x[[1]]$median_rank_stat
})%>% unlist()%>% enframe(name = "PheCode")%>% left_join(Phe_dic)

p.median_rank= df3 %>% ggplot(., aes(x= category, y= value))+
  geom_boxplot()+
  geom_jitter(alpha= 0.2)+
  coord_flip()+
  ggtitle(paste0("median_rank/ n"))



# partial AUROC -----------------------------------------------------------

df5 = map(auc_res, function(x){
  x[[1]]$pAUROC
})%>% unlist()%>% enframe(name = "PheCode")%>% left_join(Phe_dic)

p.auroc= df5 %>% ggplot(., aes(x= category, y= value))+
  geom_boxplot()+
  geom_jitter(alpha= 0.2)+
  coord_flip()+
  ggtitle(paste0("partial AUROC (FPR 0.01)"))

# partial AUROC corrected------------------------------------------------------------------
df4 = map(auc_res, function(x){
  x[[1]]$pAUROC.object.corrected$auc
})%>% unlist()%>% enframe(name = "PheCode")%>% left_join(Phe_dic)

p.pauroc.c= df4 %>% ggplot(., aes(x= category, y= value))+
  geom_boxplot()+
  geom_jitter(alpha= 0.2)+
  coord_flip()+
  ggtitle(paste0("partial AUROC (FPR 0.01), corrected"))


cowplot::plot_grid(p.fullAUROC,p.median_rank, p.pauroc.c, p.pr)



# combine plots ------------------------------------------------------------
df = rbind(df1%>% mutate(metric = "AUC-PR"),
           df2 %>% mutate(metric = "AUROC"),
           df3 %>% mutate(metric = "median_rank"),
           df4 %>% mutate(metric = "pAUROC"))
p.all= df %>%
  filter(metric!= "pAUROC") %>%
  ggplot(., aes(x= metric, y= value))+
  geom_jitter(alpha= 0.2)+
  geom_violin(alpha= 0.8)+
  geom_boxplot(width=0.15)+
  scale_y_continuous(n.breaks= 10)+
  labs(x= "")+
  theme_bw()+
  theme(axis.text = element_text(size=11, color="black"))
p.all

df %>% group_by(metric)%>% summarise(median(value))
df %>% filter(category== "circulatory system")%>% group_by(metric)%>% summarise(median(value))

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/int_val_results.pdf",
    width = 4, height= 5)
unify_axis(p.all)
dev.off()

df %>% pivot_wider(names_from= metric, values_from= value)%>%
  ggplot(aes(x= AUROC, y=   `AUC-PR`))+
  geom_point()


df= df%>% mutate(category= ifelse(category== "NULL","symptoms", category ))
p.categories= df %>%filter(metric != "pAUROC")%>%
  ggplot(., aes(x= category, y= value))+
  facet_grid(rows= vars(metric), scales = "free")+
  #geom_violin()+
  geom_boxplot(width=0.5)+
  geom_jitter(alpha= 0.2)+
  scale_y_continuous(n.breaks= 8)+
  labs(x= "")+
  theme_bw()+
  theme(axis.text = element_text(size=11, color="black"),
        axis.text.x= element_text(angle= 60, hjust= 1))

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/int_val_results_category.pdf",
    width = 4, height= 5)
unify_axis(p.categories)
dev.off()


df%>%
  group_by(category, metric)%>%
  summarise(m.v= median(value))%>%
  pivot_wider(names_from = metric, values_from = m.v)%>%
  arrange(desc(AUROC))%>%
  ungroup()%>%
  mutate(AUROC.r= rank(desc(AUROC )),
         AUCPR.r= rank(desc(`AUC-PR`)),
         medianrank.r= rank(median_rank),
         s.r= AUROC.r+AUCPR.r+medianrank.r)%>%
  arrange(desc(s.r))


# assess correlation of prediciton performance with patient number --------

## geneset size
genecount= edge.list$disease_gene %>% group_by(nodeB)%>% count

df.2= df %>% left_join(genecount%>% dplyr::rename(PheCode= nodeB), by= "PheCode")

dis.cor= df.2 %>% group_by(metric)%>% summarise(rho= cor.test(n, value)$estimate,
                                       p= cor.test(n, value)$p.value)%>%
  mutate(test= "disease genes")


#patient frequency:
pat.counts= disease_frequencies(pids.list$hf_all, icd = data%>%
                                  filter(PheCode %in% V(edge.list$disease$heidelberg)$name))

pat.cor= df %>% left_join(pat.counts, by= "PheCode")%>%
  group_by(metric)%>%
  summarise(rho= cor.test(rel_freq, value)$estimate,
            p= cor.test(rel_freq, value)$p.value)%>%
  mutate(test= "disease frequency")

## mean geneset confidence

gd= edge.list$disease_gene %>% group_by(nodeB)%>% summarise(m.weight= median(weight),
                                       mean.weight= mean(weight))

complete= left_join(df, gd %>% dplyr::rename(PheCode= nodeB), by= "PheCode")

disgenet.cor= complete%>%
 ungroup() %>%
  group_by(metric)%>%
  summarise(rho= cor.test(m.weight, value)$estimate,
            p= cor.test(m.weight, value)$p.value)%>%
  mutate(test= "DisGeNET confidence")


p.cors= rbind(dis.cor, pat.cor,disgenet.cor)%>%
  mutate(significant = ifelse(p<0.05, "yes", "no"))%>%
  filter(metric != "pAUROC")%>%
  ggplot(., aes(x= metric, y= rho , size = -log10(p), col= significant))+
  facet_grid(rows = vars(test))+
  geom_point()+
  ylim(c(-0.4, 0.4))+
  geom_hline(yintercept = 0, color= "darkgrey")+
  labs(x= "")+
  theme_bw()+
  theme(axis.text = element_text(size=11, color="black"))


p.cors
pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/hetnet/int_val_correlations.pdf",
    width = 5, height= 5)
p.cors
dev.off()


rbind(dis.cor, pat.cor,disgenet.cor)%>%
  mutate(significant = ifelse(p<0.05, "yes", "no"))%>%
  filter(metric != "pAUROC")%>%
  ggplot(., aes(x= ))

# test category dependence --------------------------------------------------------------------


df2= df%>%
  pivot_wider(names_from = metric, values_from = value)%>%
  arrange(desc(AUROC))%>%
  ungroup()
df2
df

kruskal.test(AUROC ~ category, data = df2)
kruskal.test(median_rank  ~ category, data = df2)
kruskal.test(`AUC-PR`  ~ category, data = df2)


####

df.cor= gd %>% filter(nodeB %in% disease_to_predict) %>% group_by(nodeB)%>%
  count() %>% rename(PheCode= nodeB) %>% left_join(df)

ggplot(df.cor, aes(x= value, y= n))+
  geom_point()
cor.test(df.cor$n, df.cor$value,method = "kendall")


df.overlap= gd %>% filter(nodeB %in% disease_to_predict) %>% group_by(nodeB)%>%
  rename(PheCode= nodeB) %>% left_join(df)

gsets= split( df.overlap$nodeA, df.overlap$PheCode)
