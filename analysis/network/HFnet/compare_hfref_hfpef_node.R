## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-12-09
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## compare hfref and hfpef within the diseae net
## ---------------------------

library(tidyverse)
library(igraph)
library(ComplexHeatmap)

source("~/GitHub/HF_gene_mining/R/utils/utils_network.R")
source("~/GitHub/HF_gene_mining/R/utils/utils.R")

dd.net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids.rds")

link.data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")



# compare neighbors -------------------------------------------------------


hfpef = link.data$links %>%
  filter(fisher.p.adj<0.01, (pcorr.corpor)>0.0) %>%
  filter(grepl("hfpef", disease2)) %>% arrange(desc(pcorr.corpor)) %>%
  dplyr::select(disease1, dis1_phenotype, dis2_phenotype, everything())# %>% pull(disease1)

hfref = link.data$links %>%
  filter(fisher.p.adj<0.01, (pcorr.corpor)>0.0) %>%
  filter(grepl("hfref", disease2)) %>% arrange(desc(pcorr.corpor)) %>%
  dplyr::select(disease1, dis1_phenotype, dis2_phenotype, everything())# %>% pull(disease1)




# compare node characteristics --------------------------------------------

#define main comorbidities to compare:
phe_dic= data%>%
  distinct(PheCode, Phenotype, category)%>%
  drop_na%>%
  add_row(PheCode= "hfpef",
          Phenotype= "HFpEF",
          category= "circulatory system") %>%
  add_row(PheCode= "hfref",
          Phenotype= "HFrEF",
          category= "circulatory system")



hfpef.centr= getProps_local(net = dd.net, directed = F,node =  "hfpef") %>% mutate(hf= "hfpef")
hfref.centr= getProps_local(dd.net, F, "hfref")%>% mutate(hf= "hfref")

p.node.feat= rbind(hfpef.centr, hfref.centr) %>%
  ggplot(., aes(x= measure, y= rank, fill = hf))+
  geom_col(position = "dodge")+
  scale_fill_manual(values= cols.nice)+
  theme_minimal()

p.node.feat

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/node.characteristics.pdf",
    height = 3,
    width=5)
p.node.feat
dev.off()


# check main diseases
main_disease= c("585.3", "401.1", "272.13", "250.2", "296.22", "496",
                "280.1","327.3","411.4", "hfpef", "hfref")

main_disease= c("585.3", "401.1", "272.13", "250.2", "296.22", "496",
                "280.1","327.3","411.4")
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
  left_join(phe_dic%>% rename(node= PheCode)) %>%

  ggplot(., aes(x= Phenotype, y= rank, fill = Phenotype))+
  facet_grid(vars(rows= measure))+
  geom_col(position = "dodge")+
  labs(fill = "",
       x= "Disease",
       y= "Ranking in HFnet")+
  scale_fill_manual(values= col.set)+
  theme_minimal()+
  theme(axis.text.x = element_blank())

p.node.feat


## shortest path from hfpef/hfref to major comorbidities
# create weight inverted net
dd.net.inv = dd.net
hist(E(dd.net)$weight)
#E(dd.net.inv)$weight= (1/(E(dd.net.inv)$weight-1))+1
E(dd.net.inv)$weight= 1-(E(dd.net.inv)$weight)
hist(E(dd.net.inv)$weight)

all_shortest_paths(from = "hfref",
                   graph =  dd.net.inv,
                   to= main_disease,
                   weights = E(dd.net)$weight

                   )

net_dist_uw= distances(dd.net,
                    v = main_disease,
                    to = main_disease,
                    algorithm = "unweighted",
                    )


rownames(net_dist_uw) =
  colnames(net_dist_uw) =
  phe_dic[match(rownames(net_dist_uw), phe_dic$PheCode),]%>% pull(Phenotype)
colors = structure(c("white", cols.nice[]),names= c(1:6) )# black, red, green, blue
colors = structure( c(1:6), c("white", cols.nice[]))
p.uw= Heatmap(net_dist_uw,show_row_dend = F, name= "distance", col = c("grey", "darkred", "red"))
print(p.uw)

## weighted distance=
net_dist_w= distances(dd.net.inv,
                       v = main_disease,
                       to = main_disease,
                       algorithm = "dijkstra",
)


rownames(net_dist_w) =
  colnames(net_dist_w) =
  phe_dic[match(rownames(net_dist_w), phe_dic$PheCode),]%>% pull(Phenotype)

p.w= Heatmap(net_dist_w,show_row_dend = F, cluster_rows = T , cluster_columns = T,
        name= "distance", col = c("blue", "red"))

p.w

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/main_comor_distances.pdf")
print(p.w)
print(p.uw)
dev.off()


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/main_comor_centrality.pdf",
    width= 5,
    heigh= 4)
print(p.node.feat)
dev.off()


# cluster disease net -----------------------------------------------------

net= dd.net

#lovain
module.net= modularize(net)

module.net$h.map


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/node.module.categories.pdf",
    height = 5,
    width=7)
module.net$h.map
dev.off()

g= module.net$fullnet
memberships= V(g)$group_louv
names(memberships) = V(g)$name

mems= memberships[c("hfpef", "hfref")]
node_df= igraph::as_data_frame(g, what= "vertices") %>% as_tibble()

p.box.plot=
  node_df %>%
  mutate(group_louv= factor(group_louv))%>%
  ggplot(., aes(x= group_louv, y= size))+
  geom_boxplot()+
  geom_jitter(alpha= 0.2)+
  labs(y="log10 count ")+
  theme_minimal()
p.box.plot


# distance between modules ------------------------------------------------

# 1st approach do the mean distance between all nodes of a module
net2 = module.net$fullnet
E(net2)$weight= 1-E(net2)$weight
clusts= unique(V(net)$group_louv)

clust_dists= sapply(c("hfpef", "hfref"), function(y){
  sapply(clusts, function(x){

    mean(distances(net2, y, V(net)$name[V(net2)$group_louv==x])  )
  })
})

rownames(clust_dists)= clusts
s= Heatmap(clust_dists, col = c("red", "black"))
print(s)

#  save Nodes list of clusters -------------------------------------------------------------------

library(WriteXLS)

node_df.split = split(node_df, f = node_df$group_louv)
map(names(node_df.split),function(x){
  write.csv(node_df.split[[x]],
            paste0("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/supplement/tables/hfnet_cluster/cluster_", x, ".csv"))
})




# plot clusters  ----------------------------------------------------------


g.hfpef= induced_subgraph(g, V(g)$group_louv == mems["hfpef"])
g.hfref= induced_subgraph(g, V(g)$group_louv == mems["hfref"])

as_data_frame(g.hfpef, "edges") %>% filter(to == "hfpef")



quickvisnet_clust= function(net, main = "net", save_comp = F){
  require(visNetwork)

  data <- toVisNetworkData(net)

  data$edges= data$edges %>%
    dplyr::mutate(width = weight*10) %>%
    #mutate(arrows = "to") %>%
    as_tibble
  # mutate(from = dis1_phenotype,
  #        to = dis2)

  data$nodes = data$nodes %>%
    as_tibble  %>%
    mutate(label = description,
           group = group_cat,
           value = size)

  network =visNetwork(data$nodes, data$edges, main = main,
                      height = "1000px", width = "100%") %>%
    visEdges(shadow = F,
             color = list(color = "darkgrey", highlight = "black")) %>%
    visPhysics(stabilization = T) %>%
    visEdges(smooth = T) %>%
    #addFontAwesome()%>%
    visLegend(useGroups = F,
              #addNodes = df_color,
              ncol = 1, stepY=1000) %>%
    #visHierarchicalLayout(direction = "LR", levelSeparation = 100)%>%
    visOptions(highlightNearest =  list(enabled = TRUE, algorithm = "hierarchical",
                                        degree = list(from = 1, to = 1)), collapse = TRUE) %>%
    visInteraction(dragNodes = F,
                   dragView = F,
                   zoomView = T)%>%
    visNodes(font = list(strokeWidth= 5,
                       color= "black",
                       size= 20
                     ))

  if(save_comp==T){
    network =network %>%
      visPhysics(stabilization = T) %>%
      visEdges(smooth = T) %>%
      visInteraction(dragNodes = F)
    return(network)
  }
  return(network)
}
hfpef.net = quickvisnet_clust(g.hfpef, "HFpEF cluster", F)
hfref.net= quickvisnet_clust(g.hfref, "HFrEF cluster", T)
#
#
# c(neighbors(g, "hfref"), neighbors(g, "hfpef"))
# V(g)[V(g)$name %in% c("hfref", "hfpef")]
# g.hf= induced_subgraph(g,
#                        c(neighbors(g, "hfref"),
#                          neighbors(g, "hfpef"),
#                          V(g)[V(g)$name %in% c("hfref", "hfpef")]))
#
# x= quickvisnet_clust(g.hf, "HF with 1st order neighbors", save_comp = T)



