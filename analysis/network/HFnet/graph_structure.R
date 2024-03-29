## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-05-05
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## characterize graph structure of the HFnet
## ---------------------------

library(tidyverse)
library(igraph)
library(ComplexHeatmap)
library(psych)

source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_network.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils.R")

dd.net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids.rds")
link.data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")


phe_dic= link.data$phe_dic
## CAVE: when calculating shortest path with igraph, or other metrics that rely on edge weights.
# weights can be interpreted as "costs", i.e. high weight is high distance between two nodes.
# create weight inverted net
x= E(dd.net)$weight

E(dd.net)$weight.inv = abs(max(x)-x)+0.00001
range(E(dd.net)$weight)
range(E(dd.net)$weight.inv)

# calculate df with node centralitiy measures  ----------------------------
centrality.df= enframe(igraph::degree(dd.net), value= "degree", name= "PheCode")
centrality.df$size= V(dd.net)$size
centrality.df$strength= strength(dd.net)
centrality.df= centrality.df %>% left_join(phe_dic)

#closeness (average shortest path length)
closeness.= closeness(dd.net, normalized =T ,
                      weights = (E(dd.net)$weight.inv))

cc= transitivity(dd.net, type = "weighted")
btwn= betweenness(dd.net, directed = F, weights= E(dd.net)$weight.inv, normalized = T)
stren= strength(dd.net)

centrality.df$strength= stren
centrality.df$cc= cc
centrality.df$btw= btwn
centrality.df$closeness= closeness.


assortativity_nominal(dd.net, types = as.numeric(as.factor(V(dd.net)$group_louv)), directed = F)
assortativity_nominal(dd.net, types = as.numeric(as.factor(V(dd.net)$group_cat)), directed = F)
assortativity_degree(dd.net)

centrality.df%>% arrange(desc(closeness))%>% print(n=100)


# plot relation to size ---------------------------------------------------

corr.plot = centrality.df[, c( "size", "degree", "strength", "cc", "btw", "closeness")]
p.centrality.measures= pairs.panels(corr.plot, method= "pearson",
             )

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/node_feature_corrs.pdf",
    height = 6,
    width = 6)
pairs.panels(corr.plot, method= "pearson",breaks = 25, hist.col = cols.nice[2]
)
dev.off()


# ANOVA  ------------------------------------------------------------------

#calculate ANOVAs p-value
sapply(c( "size", "degree", "strength", "cc", "btw", "closeness"), function(x){
  print(x)
  #anova= aov(formula = as.formula(paste0(x, " ~ category")), data =centrality.df)
  kruskal.test(formula = as.formula(paste0(x, " ~ category")), data =centrality.df)$p.value
  #summary(anova)[[1]][["Pr(>F)"]][1]

})


summary(anova)[[1]][["Pr(>F)"]][1]

cor.test(df$size, df$strength)
p.cent= centrality.df %>% filter(category!= "NULL")%>%
  pivot_longer(cols= c(closeness, btw, degree, size, cc),
               names_to = "metric",
               values_to = "value")%>%
  ggplot(., aes(x= category, y= value))+
    geom_boxplot()+
    facet_grid(rows= vars(metric), scales= "free_y")+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 45, hjust = 1))


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/node_feature_cats.pdf",
    height = 9,
    width = 7)
unify_axis(p.cent)+labs(y= "")
dev.off()

hist(centrality.df$degree, breaks= 50)

# calculate rank of cardiac - non cardiac links  --------------------------

#get all cardio nodes
V.df=  igraph::as_data_frame(dd.net, what = "vertices")

cardio_nodes=V.df%>%
  filter(group_cat== "circulatory system")%>%
  pull(name)


df= lapply(cardio_nodes, function(x){
  t= V.df %>% filter(name %in% names(neighbors(graph = dd.net, v = x)))
  t1 = as.data.frame(prop.table(table(t$group_cat)))
  t1%>% mutate(node= x)
})%>% do.call(rbind, .)%>% as_tibble()

df= lapply(V(dd.net)$name, function(x){
  t= V.df %>% filter(name %in% names(neighbors(graph = dd.net, v = x)))
  t1 = as.data.frame(prop.table(table(t$group_cat)))
  t1%>% mutate(node= x)
})%>% do.call(rbind, .)%>% as_tibble()

df %>%
  rename(PheCode = node )%>% left_join(phe_dic)%>% filter(category==Var1)%>%
  ggplot(., aes(x= category , y= Freq))+
  geom_boxplot()+
  coord_flip()


df%>%mutate(cl= as.factor(ifelse(node  %in% c("hfpef", "hfref"), 1, 2)),
            label= ifelse(node  %in% c("hfpef", "hfref"), node , "other"))%>%
  ggplot(., aes(x= Var1, y= Freq))+
  geom_boxplot()+
  geom_jitter(mapping = aes(), alpha= 0.2)+
  scale_color_manual(values= c( "blue", "red","black"))+
  theme_minimal()+
  #geom_label_repel(aes(label= label ), box.padding = 2, max.overlaps = 200,show.legend = FALSE)+
  theme(axis.text.x = element_text(angle=40, hjust=1))+
  labs(x= "",
       y= "% of neighbors in category")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


df%>% filter(node %in% c("hfpef", "hfref"))%>%
  ggplot(., aes(x= Var1, y= node, fill = Freq))+
  theme(axis.text.x = element_text(angle=40, hjust=1))+
  geom_tile()+
  coord_equal()


df%>%
  filter(Var1== "circulatory system")%>%
  ggplot(., aes(x= Freq))+
  geom_histogram(bins = 30)+
  geom_vline()

table(df%>% filter(Var1== "circulatory system")%>%
  group_by(node)%>%
  mutate(x= Freq<0.241 )%>% pull(x))

df.m= df%>% filter(Var1== "circulatory system")
m= mean(df.m$Freq)
sds= sd(df.m$Freq)

pnorm(q = 0.241, mean= m, sd= sds, lower.tail = T)
pnorm(q = 0.33, mean= m, sd= sds, lower.tail = T)
19/63
33/49

df.m%>% filter(Freq<0.241)


# check highes ranked nodes -----------------------------------------------

write.csv(centrality.df,  "output/Supp_tables/Supplement_table_1_centrality_rankings.csv")
centrality.df %>% arrange(desc(btw))%>% print(n=400)
centrality.df %>% arrange(desc(degree))



# check main diseases -----------------------------------------------------


# check main diseases
main_disease= c("585.3", "401.1", "272.13", "250.2", "296.22", "496",
                "280.1","327.3","411.4")

centrality.df %>%
  mutate(r.size= rank(desc(size)),
         r.degree= rank(desc(degree)),
         r.closeness = rank(desc(closeness )),
         r.btw = rank(desc(btw ))
         )%>%
  filter(PheCode %in% main_disease)

phe_dic= phe_dic%>% mutate(Phenotype= ifelse(PheCode == "280.1",
                                             "Iron deficiency anemia",
                                             Phenotype
))

node_props= map(main_disease, function(x){
  df= getProps_local(net = dd.net, directed = F,node = x) %>%
    mutate(node= x)
})%>% do.call(rbind, .)

order.ph= node_props %>%
  filter(measure %in% c("degree", "betweenness", "closeness"))%>%
  left_join(phe_dic%>% rename(node= PheCode))%>% group_by(Phenotype)%>%
  summarise(m.r= mean(rank))%>%
  arrange(desc(m.r))%>% pull(Phenotype)


p.node.feat=
  node_props %>%
  filter(measure %in% c("degree", "betweenness", "closeness"))%>%
  left_join(phe_dic%>% rename(node= PheCode)) %>%
  mutate(Phenotype= factor(Phenotype, levels= order.ph))%>%
  ggplot(., aes(y= Phenotype, x= measure, fill = rank))+
  geom_tile(col = "black")+
  geom_text(aes(label = round(rank, 1)), size= 2.5) +
  scale_fill_gradient(low= "red", high = "grey", limits= c(1, 283))+
  coord_equal()+
  theme_bw()+
  theme(axis.text.x= element_text(angle = 60, hjust= 1))

unify_axis(p.node.feat)

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/topcentral_disease.pdf",
    height = 6.5,
    width = 4)
unify_axis(p.node.feat)
dev.off()


