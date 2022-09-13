## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-12-15
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------
library(RandomWalkRestartMH)
library(tidyverse)
library(ComplexHeatmap)

source("analysis/utils/utils_hetnet.R")

edge.list= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds")

ML.class= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")




ppi_net_multiplex= create.multiplex(edge.list$gene)

dd_net_multiplex= create.multiplex(edge.list$disease)

dg_fil= edge.list$disease_gene

PPI_Disease_Net <- create.multiplexHet(ppi_net_multiplex,
                                       dd_net_multiplex,
                                       as.data.frame(dg_fil[,c(1,2,3)]), "disease")

PPIHetTranMatrix <- compute.transition.matrix(PPI_Disease_Net)


# get gene ranking for hfpef and hfref nodes ------------------------------
ML.class= ML.class %>% mutate(hf = ifelse(importance>0 & estimate< -0.2, "hfpef",
                                          ifelse(importance>0 & estimate>0.2, "hfref", "ns")
                                          )
                              )

seed.hfpef= ML.class %>% filter(hf=="hfpef")%>% pull(PheCode)
seed.hfref= ML.class %>% filter(hf=="hfref")%>% pull(PheCode)

Phe_dic%>% filter(PheCode %in% seed.hfref)

# main_disease= c("585.3", "401.1", "440", "272.13", "250.2", "296.22", "496",
#                 "280.1","327.3","411.4")

g.hfpef <-
  Random.Walk.Restart.MultiplexHet(x= PPIHetTranMatrix,
                                   MultiplexHet_Object = PPI_Disease_Net,
                                   Multiplex1_Seeds= c(),
                                   Multiplex2_Seeds = seed.hfpef,
                                   r=0.8)

g.hfref <-
  Random.Walk.Restart.MultiplexHet(x= PPIHetTranMatrix,
                                   MultiplexHet_Object = PPI_Disease_Net,
                                   Multiplex1_Seeds= c(),
                                   Multiplex2_Seeds = seed.hfref,
                                   r=0.8)



g.hfpef= g.hfpef$RWRMH_Multiplex1 %>%
  rename(value= Score,
         name= NodeNames)%>% as_tibble()

g.hfref= g.hfref$RWRMH_Multiplex1 %>%
  rename(value= Score,
         name= NodeNames)%>% as_tibble()

# plots -------------------------------------------------------------------

sets= load_validation_genes(disgenet_value= 0.29)
sets$set_phe= NULL
sets$set_reheat= NULL
sets$set_reheat_up= sets$set_reheat_up[1:250]
val.set= unique(unlist(sets))
length(unique(unlist(sets)))

pef= sapply(sets, function(x){
  res = validate_results2(x, g.hfpef)
  c("PR_AUC"= res$pr$auc.integral,
    "AUROC" =res$roc$auc)
})
ref= sapply(sets, function(x){
  res = validate_results2(x, g.hfref)
  c("PR_AUC"= res$pr$auc.integral,
    "AUROC" =res$roc$auc)
})

ref["HF"]= rep("HFrEF",2)


rbind(pef, ref)%>% Heatmap()

# PLOT RANKS --------------------------------------------------------------

g.pef= g.hfpef %>% mutate(rank= rank(desc(value)))
g.ref= g.hfref %>% mutate(rank= rank(desc(value)))

# add a simple way to prioritize genes by multiplying RW probability with rank-difference and rank new results:
g.ranks= full_join(g.pef, g.ref, by= "name") %>%
  mutate(rank.diff= rank.y-rank.x) %>%
  rename(rank.hfpef = rank.x,
         rank.hfref= rank.y,
         RW.value.hfpef= value.x,
         RW.value.hfref= value.y) %>%
  mutate(hfref.prio= RW.value.hfref * -rank.diff,
         hfpef.prio= RW.value.hfpef * rank.diff)
# get top candidates:

saveRDS(g.ranks, file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks.rds")
g.ranks= readRDS( file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks.rds")

saveRDS(list(g.hfpef, g.hfref), file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks2.rds")
g.list= readRDS( file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks2.rds")

g.ranks %>% write.csv(., file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/predicted_HF_genes.csv")


top_n = 50

g.pef.top= g.ranks %>% arrange(desc(hfpef.prio))%>% slice(1:top_n) %>% pull(name)
g.ref.top= g.ranks %>% arrange(desc(hfref.prio))%>% slice(1:top_n) %>% pull(name)

p.full.ranks= g.ranks %>% ggplot(., aes(x= rank.hfpef, y= rank.hfref))+
  geom_point()

library(ggrepel)

## PLOT RW PROBs
p.ef = ggplot(g.pef, aes(x= rank,y= value))+
  geom_point()+
  theme_classic()+
  geom_vline(xintercept = 250, col= "darkgrey")+
  labs(x= "gene ranking hfpef",
       y= "RW probability")+
  ggtitle("HFpEF")

p.ef

p.ref = ggplot(g.ref, aes(x= rank,y= value))+
  geom_point()+
  theme_classic()+
  geom_vline(xintercept = 250, col= "darkgrey")+
  labs(x= "gene ranking hfref",
       y= "RW probability")+
  ggtitle("HFrEF")

p.ref

p.cutoff= plot_grid(p.ef, p.ref, labels = "AUTO")

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/hetnet/gene.rankings.RW.prob.pdf",
    height= 5,
    width= 8)
p.cutoff
dev.off()

############ plot pef and ref gene candidates

p.pef.top250 =
  g.ranks %>%
  filter(rank.hfpef<500)%>%
  mutate(label = ifelse(name %in% g.pef.top, name, ""),
         col= ifelse(rank.diff >500, "yes", "no")) %>%
  ggplot(., aes(x= rank.hfpef,
                y= rank.hfref,
                col= hfpef.prio))+
  geom_point()+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  #geom_hline( yintercept = log10(100), type = 2, color= "grey")+
  labs(x= "ranking hfpef",
       y= "ranking hfref",
       col ="RW probability * \n ranking difference")+
  #geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()+
  scale_color_gradient(low= "darkgrey", high= "blue")+
  theme(axis.text = element_text(size= 14))+
  ggtitle("HFpEF Ranking top 250")

p.pef.top250

p.pef.top100 = g.ranks %>%
  filter(rank.hfpef<200)%>%
  mutate(label = ifelse(name %in% g.pef.top, name, ""),
         col= ifelse(rank.diff >500, "yes", "no")) %>%
  ggplot(., aes(x= rank.hfpef,
                y= rank.hfref,
                col= col))+
  geom_point()+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  #geom_hline( yintercept = log10(100), type = 2, color= "grey")+
  labs(x= "ranking hfpef",
       y= "ranking hfref",
       col ="ranking difference \n >500")+
  geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()+
  theme(axis.text = element_text(size= 14))+
  ggtitle("HFpEF Ranking top 100")

p.pef.top250
p.pef.top100

p.pef.HFgenes = g.ranks %>%
  filter(rank.hfpef<250)%>%
  mutate(label = ifelse(name %in% val.set, name, ""),
         label2 = ifelse(!name %in% val.set, "other", "hf"),
         col= ifelse(rank.diff >500, "yes", "no")) %>%
  ggplot(. ,aes(x= rank.hfpef,
                y= rank.hfref,
                col= label2))+
  geom_point( size= 2)+
  scale_color_manual(values= cols.nice[c(1,4)])+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  labs(x= "ranking hfpef",
       y= "ranking hfref",
       col ="HF gene set")+
  geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()

p.pef.HFgenes

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/manuscript/main/HFpEF.gene.rankings.dg2.pdf",
    height= 8,
    width= 12)

p.pef.top100
p.pef.top250
p.pef.HFgenes

dev.off()

########## Plot same for ref

p.ref.top250 = g.ranks %>%
  filter(rank.hfref<250)%>%
  mutate(label = ifelse(name %in% g.ref.top, name, ""),
         col= ifelse(rank.diff < (-500), "yes", "no")) %>%
  ggplot(., aes(x= rank.hfref,
                y= rank.hfpef,
                col= col))+
  geom_point()+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  #geom_hline( yintercept = log10(100), type = 2, color= "grey")+
  labs(x= "ranking hfref",
       y= "ranking hfpef",
       col ="ranking difference \n >500")+
  geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()+
  theme(axis.text = element_text(size= 14))+
  ggtitle("HFrEF Ranking top 250")

p.ref.top100 = g.ranks %>%
  filter(rank.hfref<100)%>%
  mutate(label = ifelse(name %in% g.ref.top, name, ""),
         col= ifelse(rank.diff < (-500), "yes", "no")) %>%
  ggplot(., aes(x= rank.hfref,
                y= rank.hfpef,
                col= col))+
  geom_point()+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  #geom_hline( yintercept = log10(100), type = 2, color= "grey")+
  labs(x= "ranking hfref",
       y= "ranking hfpef",
       col ="ranking difference \n >500")+
  geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()+
  theme(axis.text = element_text(size= 14))+
  ggtitle("HFrEF Ranking top 100")

p.ref.top250
p.ref.top100

p.ref.HFgenes = g.ranks %>%
  filter(rank.hfref<250)%>%
  mutate(label = ifelse(name %in% val.set, name, ""),
         label2 = ifelse(!name %in% val.set, "other", "hf"),
         col= ifelse(rank.diff >(-500), "yes", "no")) %>%
  ggplot(. ,aes(x= rank.hfref,
                y= rank.hfpef,
                col= label2))+
  geom_point( size= 2)+
  scale_color_manual(values= cols.nice[c(1,4)])+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  labs(x= "ranking hfref",
       y= "ranking hfpef",
       col ="HF gene set")+
  geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()

p.ref.HFgenes

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/manuscript/main/HFrEF.gene.rankings.dg2.pdf",
    height= 8,
    width= 12)

p.ref.top100
p.ref.top250
p.ref.HFgenes

dev.off()


p.ref.HFgenes2 = g.ranks %>%
  filter(rank.hfref<250)%>%
  mutate(label = ifelse(abs(rank.diff)<50, name, ""),
         label2 = ifelse(!name %in% val.set, "other", "hf"),
         col= ifelse(rank.diff >(-500), "yes", "no")) %>%
  ggplot(. ,aes(x= rank.hfref,
                y= rank.hfpef,
                col= label2))+
  geom_point( size= 2)+
  scale_color_manual(values= cols.nice[c(1,4)])+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  labs(x= "ranking hfref",
       y= "ranking hfpef",
       col ="HF gene set")+
  geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()

p.ref.HFgenes2
# compare network modules -------------------------------------------------


monet= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/gene.modules.MONET.M1.rds")

names(monet)= paste0("c.", seq(1:length(monet)))

monet_length= lapply(monet, length)

names.5= names(monet_length[unlist(map(monet_length,function(x)(x>5)))])

enframe(monet, name= "cluster", value= "name")%>%
  unnest(name)

after_clustering =
  g.ranks %>%
  left_join(enframe(monet, name= "cluster", value= "name")%>%
              unnest(name))%>% print(n=100)

after_clustering%>%
  arrange(desc(hfpef.prio))

t= after_clustering%>% group_by(cluster)%>% mutate(mean_diff= median(hfpef.prio))%>%
  filter(cluster %in% names.5)

t%>% arrange(desc(mean_diff)) %>% distinct(cluster, mean_diff) %>% print(n=100)


after_clustering %>% filter(cluster==  "c.373")%>% print(n=100)


t%>% distinct(cluster, mean_diff) %>% arrange(desc(mean_diff))%>% print(n=500
                                                                        )
after_clustering%>% filter(grepl("ANGP", name))
g.ranks %>%
  arrange(desc(hfpef.prio)) %>% slice(1:50)




