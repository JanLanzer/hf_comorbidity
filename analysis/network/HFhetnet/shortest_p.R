## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-06-15
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Get shortest paths for visualization

library(tidyverse)
library(igraph)
edges= readRDS("~/Downloads/multilayer_edge_list.rds")
edges= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds")

gene_df = map(edges$gene, function( x){igraph::as_data_frame(x,what= "edges")})
gene_df = do.call(what = rbind, gene_df)
dis_df = map(edges$disease, function( x){igraph::as_data_frame(x,what= "edges")})
dis_df2 = do.call(what=rbind, dis_df)
colnames(edges$disease_gene) = c("from", "to", "weight")

df= rbind(dis_df2[, 1:3], gene_df[,1:3],edges$disease_gene)
rownames(df)= NULL
dim(df)
write_csv(df, "data/hetnet_edges.csv")
write_csv(df, "data/hetnet_edges.csv")

graph= graph_from_data_frame(df[,1:2], directed = F)

#genes= read_csv(file ="data/predicted_HF_genes.csv")
comp_feat= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")
y= comp_feat%>% arrange((estimate))%>% slice(1:20)%>% pull(PheCode)
genes= read_csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/predicted_HF_genes.csv")
x= genes %>% arrange(desc(hfpef.prio))%>% slice(1:50)%>% pull(gene)

s= all_shortest_paths(graph, from = y, to = x)

unlist(s$res)

ss= igraph::induced_subgraph(graph, v= names(unlist(s$res)))



plot(ss)
