## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-12-08
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## compare HF net with morbinet
## ---------------------------
library(ComplexHeatmap)
library(igraph)
library(tidyverse)

source("~/GitHub/HF_gene_mining/R/utils/utils_hetnet.R")
source("~/GitHub/HF_gene_mining/R/utils/utils_network.R")

link.data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")
hf_net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
edgelist= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds")
Phedic= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/icd10_phewas_dictionary.rds")
hpo_net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/hpo_net.rds" )
hpo_net_full = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/hpo_net_full.rds" )

Phedic= link.data$phe_dic
# we compare  -------------------------------------------------------------

# 1 HF_net (normal)
# 2 HF_net (OR cut off)
# 3 hpo net
# 4 hpo net reduced
# 5 morbinet
# 6 morbinet reduced

# 1 HF_net the way it is implemented
hf_nodes= c("hfpef", "hfref")
hf_real = hf_net

hf_real= induced_subgraph(hf_real,
                          !V(hf_real)$name %in% hf_nodes)

# 2 HF_net with OR
hf_or = table_to_links( link.data$links ,
                     p.val = 0.55,
                     weight_dd = 2,
                     weight_col = "odds.ratio")
hf_or= hf_or %>% filter(!nodeA %in% hf_nodes,
                  !nodeB %in% hf_nodes)

hf_or= simplify_link_table(hf_or) %>%
  graph_from_data_frame(., directed = F)


#3) full hpo net
hpo_net_full
hpo_net_full =  simplify_link_table(as_data_frame(hpo_net_full))%>%
  graph_from_data_frame(., directed = F)

# 4) reduced hpo net
hpo_net =  simplify_link_table(as_data_frame(hpo_net))%>%
  graph_from_data_frame(., directed = F)


#5) full morbinet
morbinet= process_morbinet()
morbinet= simplify_link_table(morbinet)
morbinet= graph_from_data_frame(morbinet, directed = F)

#6) reduced morbinet
morbinet_r= induced_subgraph(
  morbinet,
  V(morbinet)$name %in% V(hf_real)$name)

#7) reduced hf_net to morbinet
hf_real_r= induced_subgraph(
  hf_real,
  V(hf_real)$name %in% V(morbinet_r)$name)

hf_real_r2= induced_subgraph(
  hf_real,
  V(hf_real)$name %in% V(hpo_net)$name)
# compare sizes  ----------------------------------------------------------

net.list= list("hfnet"= hf_real,
         "hf_or"= hf_or,
         "hfnet_morbi_set"= hf_real_r,
         "hfnet_hpo_set"= hf_real_r2,
          "hpo_full"= hpo_net_full,
          "hpo_subset"= hpo_net,
         "morbinet_full"= morbinet,
         "morbinet_subset"= morbinet_r)

net.sizes= sapply(net.list, function(x){
  tibble("nodes"= length(V(x)),
         "edges"= length(E(x))
  )
})


p.network.sizes=
  net.sizes%>% as.data.frame() %>%rownames_to_column("feature") %>%
  pivot_longer(-feature)%>% mutate(value = unlist(value))%>%
  ggplot(., aes(x= name, y= value))+
  facet_grid(vars(cols= feature),  scales="free")+
  geom_col(fill = cols.nice[1])+
  theme_light()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle= 60, hjust= 1))



# compare jaccard overlap: ------------------------------------------------

vertex.sim=
  sapply(net.list, function(x){
  map(net.list, function(y){
    jaccard_index(x, y)

  })%>% unlist()
})

edge.sim =
  sapply(net.list, function(x){
    map(net.list, function(y){
      jaccard_index(x, y, "edge")

    })%>% unlist()
  })


v.map = Heatmap(vertex.sim, name = "jaccard",
        col = cols.nice[c(3,2)],
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", vertex.sim[i, j]), x, y, gp = gpar(fontsize = 8))
        },
        show_row_dend = F)


e.map = Heatmap(edge.sim, name = "jaccard",
                col = cols.nice[c(4,1)],
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.1f", edge.sim[i, j]), x, y, gp = gpar(fontsize = 8))
                },show_row_dend = F

                )


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/disease_net_comparison_heatmaps.pdf",
    width= 5,
    height = 5)
print(v.map)
print(e.map)
dev.off()

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/disease_net_comparison_size.pdf",
    width= 4,
    height = 5)
print(p.network.sizes)
dev.off()


# us

# Show edges conserved between both:  -------------------------------------

g1= net.list$hf_or
g2= net.list$morbinet_subset
g1_edg= igraph::as_data_frame(g1, what= "edges") %>% select(from,to) %>% as_tibble
g2_edg = igraph::as_data_frame(g2, what= "edges")%>% select(from,to)%>% as_tibble

# create an edgeID that is the same in both networks (sort the 2 nodes and than paste)
g1_IDs= map2(g1_edg$from, g1_edg$to, function(x,y){
  paste(sort(c(x,y)),collapse="_")
})
g1_edg= g1_edg %>% mutate(edgeID = unlist(g1_IDs))

g2_IDs= map2(g2_edg$from, g2_edg$to, function(x,y){
  paste(sort(c(x,y)),collapse="_")
})
g2_edg= g2_edg %>% mutate(edgeID = unlist(g2_IDs))

# now calculate jaccard based on this edge id union and intersect
u= union(g1_edg$edgeID, g2_edg$edgeID)
i= intersect(g1_edg$edgeID, g2_edg$edgeID)

length(g2_edg$edgeID)
length(i)
# hf net edges that intersect:
g1_edg %>% filter(edgeID %in% i) %>%
  left_join(Phedic%>% rename(from= PheCode,
                              Phenotype_from= Phenotype)%>% select(-category)) %>%
  left_join(Phedic%>% rename(to= PheCode,
                              Phenotype_to= Phenotype)%>% select(-category)) %>% print(n=500)

# hf net special edges:
g1_edg %>% filter(!edgeID %in% i) %>%
  left_join(Phedic%>% rename(from= PheCode,
                              Phenotype_from= Phenotype)%>% select(-category)) %>%
  left_join(Phedic%>% rename(to= PheCode,
                              Phenotype_to= Phenotype)%>% select(-category)) %>% print(n=500)
