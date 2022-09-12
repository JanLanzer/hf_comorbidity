## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-01-28
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## create networks with top 8000 edges from tetrachoric and partial tetrachoric coefficients
## and compare network transitivity
## ---------------------------

library(tidyverse)
library(igraph)
library(ggpubr)

source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_network.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_hetnet.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils.R")


data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")
Phe_dic= data$phe_dic

### 1. Tetrachoric
tet= data$links%>%
  filter(fisher.p.adj<0.05)%>%
  arrange(desc(corr.tet))%>%
  top_n(.,n = 6000, wt = corr.tet)

tet= data$links%>%
  filter(fisher.p.adj<0.05, corr.tet>0)%>%
  arrange(desc(corr.tet))

dd.net= quick_base_net(tet, data$pid, data$data, "corr.tet")

### 2. Partial Tetrachoric
ptet= data$links%>%
  #filter(fisher.p.adj<0.01)%>%
  arrange(desc(pcorr.corpor))%>%
  top_n(.,n = 6000, wt = pcorr.corpor)

ptet= data$links%>%
  filter(rho_part_lower>0)%>%
  arrange(desc(pcorr.corpor))

dd.pnet= quick_base_net(links = ptet,pids =  data$pid, data = data$data, "pcorr.corpor")


### 3. compare
graph_centrality= function(dd.net, phe_dic){

  # create weight inverted net
  dd.net.inv = dd.net
  #E(dd.net.inv)$weight= (1/(E(dd.net.inv)$weight-1))+1

  plot(E(dd.net.inv)$weight, E(dd.net)$weight)
  # correlated measures?
  #degree and size
  df= enframe(igraph::degree(dd.net), value= "degree")
  df$size= V(dd.net)$size
  df$strength= strength(dd.net)

  #closeness (average shortest path length)
  cl= enframe(closeness(dd.net.inv, normalized = T,
                        weights = NULL), value= "closeness")

  cc= transitivity(dd.net, type = "local",weights = E(dd.net)$weight)
  names(cc)= V(dd.net)$name
  cc= enframe(cc, value= "cc")
  #cc2= transitivity(dd.net, type = "local")


  btwn= betweenness(dd.net.inv, directed = F, weights= NULL, normalized = T)
  btwn = enframe(btwn, value= "btwn")
  #merge
  df = left_join(df, cl)%>% left_join(cc)%>% left_join(btwn)%>%
    rename(PheCode= name)%>% left_join(phe_dic)
}


pt.c= graph_centrality(dd.net = dd.pnet, phe_dic = data$phe_dic)
t.c= graph_centrality(dd.net,phe_dic = data$phe_dic)

pt.c= pt.c %>% filter(PheCode %in% t.c$PheCode)
t.c= t.c %>% filter(PheCode %in% pt.c$PheCode)

wilcox.test(pt.c$cc,  t.c$cc, paired = T)
cc= transitivity(dd.pnet, type = "global")
cc2= transitivity(dd.net, type = "global")

boxplot(  pt.c$cc,  t.c$cc, names = c("partial rho", "rho"))

boxplot(  pt.c$btwn ,  t.c$btwn , names = c("partial rho", "rho"))

boxplot(  pt.c$closeness ,  t.c$closeness , names = c("partial rho", "rho"))

boxplot(  pt.c$degree ,  t.c$degree , names = c("partial rho", "rho"))


plot.df= rbind(pt.c %>% mutate(rho="partial tet"),
      t.c %>% mutate(rho="tet"))

my_comparisons <- list( c(" partial.tet", "tet") )

p.local.cc.comp= ggplot(plot.df , aes(x= rho, y= cc))+
  geom_violin()+
  geom_jitter(alpha= 0.2)+
  geom_boxplot( width= 0.3)+
  theme_classic()+
  labs(y= "local cluster coefficient")+
  ylim(c(0, 1.2))+
  stat_compare_means(label = "p.signif", label.y = 1.1)

p.local.cc.comp

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/CC_comparison.pdf",
    width= 3,
    heigh= 3)
p.local.cc.comp
dev.off()



p.btwn= ggplot(plot.df , aes(x= rho, y= btwn))+
  geom_violin()+
  geom_jitter(alpha= 0.2)+
  geom_boxplot( width= 0.3)+
  theme_classic()+
  labs(y= "local btwn")+
  #ylim(c(0, 1.2))+
  stat_compare_means(label = "p.signif")

p.btwn

p.degree= ggplot(plot.df , aes(x= rho, y= degree))+
  geom_violin()+
  geom_jitter(alpha= 0.2)+
  geom_boxplot( width= 0.3)+
  theme_classic()+
  labs(y= "local btwn")+
  #ylim(c(0, 1.2))+
  stat_compare_means(label = "p.signif")

p.btwn


p.btwn= ggplot(plot.df , aes(x= rho, y= closeness))+
  geom_violin()+
  geom_jitter(alpha= 0.2)+
  geom_boxplot( width= 0.3)+
  theme_classic()+
  labs(y= "local closeness")+
  #ylim(c(0, 1.2))+
  stat_compare_means(label = "p.signif")
p.btwn


ggplot(plot.df , aes( x= degree))+
  geom_histogram(bins = 50)+
  facet_grid(rows= vars(rho), scales = "free_y")

  geom_violin()+
  geom_jitter(alpha= 0.2)+
  geom_boxplot( width= 0.3)+
  theme_classic()+
  labs(y= "local closeness")+
  #ylim(c(0, 1.2))+
  stat_compare_means(label = "p.signif")
