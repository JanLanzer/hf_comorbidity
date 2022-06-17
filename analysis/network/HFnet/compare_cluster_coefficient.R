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

source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_network.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_hetnet.R")
data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")


### 1. Tetrachoric

tet= data$links%>%
  filter(fisher.p.adj<0.01)%>%
  arrange(desc(corr.tet.sirt))%>%
  top_n(.,n = 8000, wt = corr.tet.sirt)

dd.tet= table_to_links(tet,
                   p.val =  1,
                   weight_dd = 0.0,
                   weight_col = "corr.tet.sirt")


# remove duplicates and select largest connected component:
dd.tet= simplify_link_table(dd.tet)
dd.tet= dd.tet %>% rename(disease1= nodeA,
                  disease2= nodeB)

dd.net= quick_base_net(dd.tet, pids, data$data, "corr.tet.sirt")

### 2. Partial Tetrachoric
ptet= .links%>%
  filter(fisher.p.adj<0.01)%>%
  arrange(desc(pcorr.corpor))%>%
  top_n(.,n = 8000, wt = pcorr.corpor)

dd.ptet= table_to_links(ptet,
                       p.val =  1,
                       weight_dd = 0.0,
                       weight_col = "pcorr.corpor")


# remove duplicates and select largest connected component:
dd.ptet= simplify_link_table(dd.ptet)
dd.ptet= dd.ptet %>% rename(disease1= nodeA,
                          disease2= nodeB)

dd.pnet= quick_base_net(links = dd.ptet,pids =  pids, data = data$data, "pcorr.crpor")


### 3. compare
graph_centrality= function(dd.net, phe_dic){

  # create weight inverted net
  dd.net.inv = dd.net
  E(dd.net.inv)$weight= 1-E(dd.net.inv)$weight
  E(dd.net.inv)$weight= 1

  # correlated measures?
  #degree and size
  df= enframe(degree(dd.net), value= "degree")
  df$size= V(dd.net)$size
  df$strength= strength(dd.net)

  #closeness (average shortest path length)
  cl= enframe(closeness(dd.net.inv, normalized = T,
                        weights = NULL), value= "closeness")

  cc= transitivity(dd.net, type = "local")
  names(cc)= V(dd.net)$name
  cc= enframe(cc, value= "cc")
  #cc2= transitivity(dd.net, type = "local")


  btwn= betweenness(dd.net.inv, directed = F, weights= NULL, normalized = T)
  btwn = enframe(btwn, value= "btwn")
  #merge
  df = left_join(df, cl)%>% left_join(cc)%>% left_join(btwn)%>%
    rename(PheCode= name)%>% left_join(phe_dic)
}

pt.c= graph_centrality(dd.pnet, phe_dic = data$phe_dic)
t.c= graph_centrality(dd.net,phe_dic = data$phe_dic)

wilcox.test(pt.c$cc,  t.c$cc, paired = T)
cc= transitivity(dd.pnet, type = "global")
cc2= transitivity(dd.net, type = "global")

boxplot(  pt.c$cc,  t.c$cc, names = c("partial.rho", "rho"))

plot.df= rbind(pt.c %>% mutate(rho="partial tet"),
      t.c %>% mutate(rho="tet"))

p.local.cc.comp= ggplot(plot.df , aes(x= rho, y= cc))+
  geom_violin()+
  geom_jitter(alpha= 0.2)+
  theme_classic()+
  labs(y= "local cluster coefficient")

p.local.cc.comp

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/CC_comparison.pdf",
    width= 3,
    heigh= 3)
p.local.cc.comp
dev.off()

