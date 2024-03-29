## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-08
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## script to modify the way we compare hfref hfpef,
## instead of creating a network for each cohort, I will create hfpef and hfref nodes based and analyze their topology within a larger hf_net
## ---------------------------


source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_network.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_hetnet.R")

library(tidyverse)
library(cowplot)
library(ggrepel)
library(igraph)
library(RandomWalkRestartMH)

.links= readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_table_hf_cohort_fil.rds")
link.data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")

pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")

hpo_net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/hpo_net_full.rds" )
HFnet= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
hpo_net= induced_subgraph(hpo_net, vids = V(hpo_net)$name %in% V(HFnet)$name)
# test with gene pred: ----------------------------------------------------

go =process_GO_layer("MF")
goBP= process_GO_layer("BP")
omni= process_omnipath()

ppi= process_human_ppi_multiverse()
ppi3= process_ppi_lit()
ppi = rbind(ppi, ppi3)

dg= process_gd(weight_cutoff = 0.29)

dd2= simplify_link_table(igraph::as_data_frame(hpo_net))
dd= simplify_link_table(igraph::as_data_frame(HFnet))

phecodes.hpo= unique(c(dd2$nodeA, dd2$nodeB))
phecodes.hdnet= unique(c(dd$nodeA, dd$nodeB))

saveRDS(list(phecodes.hdnet, phecodes.hpo),
        "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/phecodes.dd_hetnet.rds")

gene_networks= list("goMF"= go,
                    "goBP" = goBP,
                    "PPI"= ppi,
                    "Omnipath"= omni)


# create subset of cardiac expressed genes --------------------------------

# add three layers of evidence for gene expression in heart:
gtex= get_GTEX_heart_genes()
hprot= get_human_heart_proteome()
reheat= as_tibble(read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/METArank_March2020.csv"))
g.list= list("GTEX"= gtex, "Proteome"= hprot, "ReHeaT"= reheat%>% pull(gene))

#plot:
ggVennDiagram::ggVennDiagram(g.list)

#union
hg= unique(c(gtex, hprot,reheat%>% pull(gene)))
length(hg)

#subset layers:
gene_networks_heart= lapply(gene_networks, function(x){
  x %>% filter(nodeA %in% hg & nodeB %in% hg)
})

map(gene_networks, dim)

map(gene_networks_heart, dim)


# compute multiplex net ---------------------------------------------------

#create two networks:

gene_graphs= map(gene_networks_heart, function(x) (graph_from_data_frame(x, directed = F)))

dd_net = graph_from_data_frame(dd, directed = F)
dd_net2 = graph_from_data_frame(dd2, directed = F)

ppi_net_multiplex= create.multiplex(gene_graphs)

dd_net_multiplex= create.multiplex(list("heidelberg"= dd_net,
                                        "hpo"= dd_net2))

genes= map(seq(1:length(gene_graphs)), function(x){
  V(ppi_net_multiplex[[x]])$name
})%>% unlist() %>% unique
length(genes)

diseases= map(seq(1:2), function(x){
  V(dd_net_multiplex[[x]])$name
})%>% unlist() %>% unique



dg_fil= dg %>% filter(nodeA %in% genes,
                      nodeB %in% diseases) %>%
  dplyr::select(nodeA, nodeB, weight)

length(unique(dg_fil$nodeA))
length(unique(dg_fil$nodeB))

saveRDS(list("gene"= gene_graphs,
             "disease"= list("heidelberg"= dd_net,
                             "hpo"= dd_net2),
             "disease_gene"= dg_fil),
        "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds"
)


PPI_Disease_Net <- create.multiplexHet(ppi_net_multiplex,
                                       dd_net_multiplex,
                                         as.data.frame(dg_fil[,c(1,2,3)]), "disease")

PPIHetTranMatrix <- compute.transition.matrix(PPI_Disease_Net)
