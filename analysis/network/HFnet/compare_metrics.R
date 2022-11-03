## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-08-29
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
##compare pairwisemetrics
## ---------------------------

.links= readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_table_hf_cohort_fil.rds")
link.data= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")


get_edges_by_metric= function(fullnet,
                              cutoff_metric){


edge.module_net = as_tibble(igraph::as_data_frame(fullnet, what = "edges"))
nodes.module_net = as_tibble(igraph::as_data_frame(fullnet, what = "vertices"))

edge_IDs= map2(edge.module_net$from, edge.module_net$to, function(x,y){
  paste(sort(c(x,y)),collapse="_")
}) %>% unlist()

}

p.rho= quick_base_net(.links%>% filter(rho_part_lower>0),#, fisher.p.adj<0.05),
                        pids,
                        data,
                        weight_col ="pcorr.corpor") %>% get_edges_by_metric()

rho= quick_base_net(.links%>% filter(corr.tet>0, fisher.p.adj<0.05),
                    pids,
                    data_comb,
                    weight_col ="pcorr.corpor") %>% get_edges_by_metric()

OR= quick_base_net(.links%>% filter(odds.ratio>1, fisher.p.adj<0.1),
                    pids,
                    data_comb,
                    weight_col ="pcorr.corpor") %>% get_edges_by_metric()

OR2= quick_base_net(.links%>% filter(corr.phi>0),
                   pids,
                   data_comb,
                   weight_col ="pcorr.corpor") %>% get_edges_by_metric()

library(ggVennDiagram)

ggvenn::ggvenn(list("part"= p.rho, "tet"= rho))


p1= ggplot(.links %>%
         mutate(edge= ifelse(rho_part_lower>0, "y", "n")), aes(x= pcorr.corpor, y= corr.tet, col = edge))+
  geom_point(size= 0.1)

p2= ggplot(.links%>%
         mutate(edge= ifelse((corr.tet>0 & fisher.p.adj<0.05), "y", "n")), aes(x= pcorr.corpor, y= corr.tet,col = edge))+
  geom_point(size= 0.1)


cowplot::plot_grid(p1, p2)
