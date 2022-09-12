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

source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_network.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils.R")

dd.net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids.rds")
link.data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")
dd.net= net
phe_dic= link.data$phe_dic
## CAVE: when calculating shortest path with igraph, or other metrics that rely on edge weights.
# weights can be interpreted as "costs", i.e. high weight is high distance between two nodes.
# create weight inverted net
dd.net.inv = dd.net
#E(dd.net.inv)$weight= (1/(E(dd.net.inv)$weight-1))+1
x= E(dd.net.inv)$weight
1/(1-x)

E(dd.net.inv)$weight = 1/(1-x)



# calculate df with node centralitiy measures  ----------------------------

centrality.df= enframe(igraph::degree(dd.net), value= "degree", name= "PheCode")
centrality.df$size= V(dd.net)$size
centrality.df$strength= strength(dd.net)
centrality.df= centrality.df %>% left_join(phe_dic)

#closeness (average shortest path length)
closeness.= closeness(dd.net.inv, normalized = T,
                      weights = E(dd.net.inv)$weight)

cc= transitivity(dd.net, type = "weighted")
btwn= betweenness(dd.net.inv, directed = F, weights= E(dd.net)$weight, normalized = T)
stren= strength(dd.net)

centrality.df$strength= stren
centrality.df$cc= cc
centrality.df$btw= btwn
centrality.df= centrality.df %>%  left_join(phe_dic)
centrality.df$closeness= closeness.

df2= centrality.df%>% drop_na()#%>% mutate(category= factor(category))

assortativity_nominal(dd.net, types = V(dd.net)$group_cat, directed = F)
?assortativity_nominal
# plot relation to size ---------------------------------------------------
library(psych)
corr.plot = centrality.df[, c( "size", "degree", "strength", "cc", "btw", "closeness")]
p.centrality.measures= pairs.panels(corr.plot, method= "pearson",
             )
pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/node_feature_corrs.pdf",
    height = 6,
    width = 6)
pairs.panels(corr.plot, method= "pearson",
)
dev.off()

#calculate ANOVAs p-value
sapply(c( "size", "degree", "strength", "cc", "btw", "closeness"), function(x){
  print(x)
  anova= aov(formula = as.formula(paste0(x, " ~ category")), data =centrality.df)
  summary(anova)[[1]][["Pr(>F)"]][1]
})

summary(anova)[[1]][["Pr(>F)"]][1]

cor.test(df$size, df$strength)
p.degree=
  ggplot(centrality.df, aes(x= size, y= degree, col = category))+
  geom_point()+
  labs(x= "log10 count")
p.degree

p.degree=
  ggplot(centrality.df, aes(x= size, y= closeness, col = category))+
  geom_point()+
  labs(x= "log10 count")
p.degree


ggplot(centrality.df, aes(x= size, y= closeness))+
  geom_point()+
  theme(axis.text.x = element_text(angle= 60, hjust= 1))


p.dis.cats= plot_grid(
  ggplot(centrality.df, aes(x= category, y= closeness))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle= 60, hjust= 1))+
    labs(x= ""),

 ggplot(centrality.df, aes(x= category, y= size))+
   geom_boxplot()+
   theme(axis.text.x = element_text(angle= 60, hjust= 1))+
   labs(x= "",
        y= ),

ggplot(centrality.df, aes(x= category, y= degree))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle= 60, hjust= 1))+
  labs(x= ""),

ggplot(centrality.df, aes(x= category, y= btw))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle= 60, hjust= 1))+
  labs(x= ""),

ggplot(centrality.df, aes(x= category, y= cc))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle= 60, hjust= 1))+
  labs(x= ""),

nrow = 3
)
pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/node_feature_cats.pdf",
    height = 10,
    width = 8)
p.dis.cats
dev.off()


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

write.csv(centrality.df,  "output/Supp_tables/node_features.csv")
centrality.df %>% arrange(desc(btw))%>% print(n=400)
centrality.df %>% arrange(desc(degree))



# check main diseases -----------------------------------------------------


# check main diseases
main_disease= c("585.3", "401.1", "272.13", "250.2", "296.22", "496",
                "280.1","327.3","411.4", "hfpef", "hfref")
phe_dic %>% filter(PheCode %in% main_disease)

phe_dic= phe_dic%>% mutate(Phenotype= ifelse(PheCode == "280.1",
                                             "Iron deficiency anemia",
                                             Phenotype
))

node_props= map(main_disease, function(x){
  df= getProps_local(net = dd.net, directed = F,node = x) %>%
    mutate(node= x)
})%>% do.call(rbind, .)

p.node.feat=
  node_props %>%
  filter(measure %in% c("degree", "betweenness", "closeness"))%>%
  left_join(phe_dic%>% rename(node= PheCode)) %>%
  mutate(rank= 352-rank)%>%
  ggplot(., aes(y= Phenotype, x= measure, fill = rank))+
  geom_tile()+
  scale_fill_gradient(low= "white", high = "red")+
  coord_equal()+
  theme_bw()+
  theme(axis.text.x= element_text(angle = 60, hjust= 1))

p.node.feat

