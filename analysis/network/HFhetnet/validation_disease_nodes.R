## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-24
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##   perform internal validation, by removing on disease-geneset association and try predict it back.
##
## ---------------------------

#source("~/GitHub/RWH_analysis/scripts/network_scripts/create_networ/create_network_links_v2.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_network.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_hetnet.R")

library(tidyverse)
library(cowplot)
library(ggrepel)
library(igraph)
library(RandomWalkRestartMH)

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids.rds")
phe_dic= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/icd10_phewas_dictionary.rds")
edge.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds")

## create multiplex networks:
ppi_net_multiplex= create.multiplex(edge.list$gene)

dd_net_multiplex= create.multiplex(edge.list$disease)

disease= dd_net_multiplex$Pool_of_Nodes
genes = ppi_net_multiplex$Pool_of_Nodes

dg = edge.list$disease_gene

#get diseases with more than 2 genes
disease_to_predict= dg %>% group_by(nodeB)%>% count %>% filter(n>2)%>% pull(nodeB)

#plot the size of disgenet gene sets to predict:

dg_fil = dg %>% filter(nodeB %in% disease_to_predict)

p_dg_size= dg_fil %>% group_by(nodeB)%>% count%>%mutate(dum= "DisGeNET")%>%
  ggplot(., aes(y= n, x= dum))+
  geom_violin()+
  geom_boxplot(width= 0.5)+
  scale_y_log10()+
  labs(x= "",
       y= "gene set size")+
  theme_minimal()


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/hetnet/disgenet_geneset_size.pdf",
    height= 3,
    width= 2)
unify_axis(p_dg_size)
dev.off()

# get diesease similarity map -------------------------------------------------------------------------

library(qdapTools)
library(pheatmap)


jacc.dist.phecodes= get_Disgenet_overlaps(dg_fil)

p.disgenet.jacc= Heatmap(jacc.dist.phecodes,
        name= "Jaccard Index",
        show_row_dend = F,
        show_column_names = F,
        show_row_names = F)

# p.disgenet.jacc= pheatmap(as.matrix(jacc.dist.phecodes),
#                           show_rownames = F,
#                           show_colnames = F,
#                           legend_labels = "Jaccard Index")

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/hetnet/disgenet_disease_distance.pdf",
    height= 10,
    width= 10)
print(p.disgenet.jacc)
dev.off()

dist1= as.matrix(jacc.dist.phecodes)
hist(dist1)
saveRDS(dist1, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/disgnete_disease_distance_matrix.rds")


# now by using a cut off we can select diseases to be removed:

cutoff = 0.5

#get the exclusion vecotor for every disease:
l.vec= map(rownames(dist1[,]), function(x){
  disease_to_remove= names(dist1[x, dist1[x, ]<cutoff])
  (length(disease_to_remove))

})

names(l.vec)= rownames(dist1)

excluded_disease_per_disease= sort(l.vec%>% unlist()) %>%
  enframe(name = "PheCode") %>%
  left_join(phe_dic%>% distinct(PheCode, Phenotype))%>%
  arrange(desc(value))%>%
  print(n=100)

p.hist= excluded_disease_per_disease %>%
  ggplot(., aes(x= value))+
  geom_histogram(bins= 20)+
  theme_minimal()+
  labs(x= paste("count of exlcuded diseases at jaccard.distance >", cutoff))

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/hetnet/disgenet_disease_distance_hist.pdf",
    height= 5,
    width= 5)
p.hist
dev.off()


# main auc.run=  ----------------------------------------------------------
node.list.genes=lapply(names(edge.list$gene), function(x){
  V(edge.list$gene[[x]])$name
})
names(node.list.genes)=  names(edge.list$gene)

auc_res= map(disease_to_predict, function(x){

  test_genes= dg_fil %>% filter(nodeB == x) %>% pull(nodeA)

  if(length(test_genes)< 1){
    return(NULL)
  }

  # use distance matrix to remove disease that share a lot of genes)
  disease_to_remove= names(dist1[x, dist1[x, ]<cutoff])
  disease_to_remove= unique(c(disease_to_remove, x))

  #remove critical diseases from the bipartite mapping:
  dg_2= dg_fil %>% filter(!nodeB %in% disease_to_remove)


  #merge all network_layers=
  PPI_Disease_Net <- create.multiplexHet(ppi_net_multiplex,
                                         dd_net_multiplex,
                                         as.data.frame(dg_2[,c(1,2,3)]), "disease")

  # calculate transitionmatrix
  PPIHetTranMatrix <- compute.transition.matrix(PPI_Disease_Net)

  # run RW with query disease as seed
  RWRH_PPI_Disease_Results <-
    Random.Walk.Restart.MultiplexHet(x= PPIHetTranMatrix,
                                     MultiplexHet_Object = PPI_Disease_Net,
                                     Multiplex1_Seeds= c(),
                                     Multiplex2_Seeds = x,
                                     r=0.8)

  g.test= RWRH_PPI_Disease_Results$RWRMH_Multiplex1 %>%
    dplyr::rename(value= Score,
           name= NodeNames)%>% as_tibble()

  test_res= validate_results_partial(set= test_genes, gene_results = g.test, FPR_cutoff = 0.02)

  # add a characterization of gene set prediction
  dg_nodes = unique(dg_2$nodeA)

  gene_coverage_by_layer= enframe(test_genes, value = "name", name= "num")%>%
    dplyr::select(-num)%>%
    mutate(dg= ifelse(name %in% dg_nodes, 1, 0),
           goMF= ifelse(name %in% node.list.genes$goMF, 1, 0),
           PPI= ifelse(name %in% node.list.genes$PPI, 1, 0),
           Omni= ifelse(name %in% node.list.genes$Omi, 1, 0),
           goBP= ifelse(name %in% node.list.genes$goBP, 1, 0)
           )%>%
    pivot_longer(names_to = "genelayer", values_to = "member", -name)%>%
    group_by(genelayer)%>%
    summarise(sum= sum(member)) %>%
    mutate(percentage= sum/length(test_genes))

  return(list("test_results" =test_res,
              "dis_removed"= disease_to_remove,
              "test_geneset_coverage"= gene_coverage_by_layer
              ))

})

names(auc_res)= disease_to_predict


saveRDS(auc_res,"T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/aurocs_multilayer_disease_pred_heart_genenet2022.rds")





