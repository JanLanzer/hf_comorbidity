## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-14
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##  Randomize layers to assess importance of each layer
##
## ---------------------------

library(tidyverse)
library(igraph)
library(RandomWalkRestartMH)
# Load Network layers and data --------------------------------------------

# D-D

source("analysis/utils/utils_hetnet.R")

edges =
  readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds")

distance_disease=
  readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/disgnete_disease_distance_matrix.rds")

# Prepare data objects for the tests: -------------------------------------

# the set of diseases, for which we will aim to recover genes in the the hetnet:
set.seed(20)

dg= edges$disease_gene %>% filter(nodeB %in% V(edges$disease$heidelberg)$name)

disease_sets = sample(unique(dg$nodeB), replace = F)

# new ---------------------------------------------------------------------

get_disease_pred_results= function(disease_to_predict,
                                   ppi_net_multiplex,
                                  dd_net_multiplex,
                                  dg_fil,
                                  cutoff= 0.7,
                                  dis_distance){

  #filter disgenet:
  dg_fil= dg_fil %>% filter(nodeA %in% ppi_net_multiplex$Pool_of_Nodes,
                            nodeB %in% dd_net_multiplex$Pool_of_Nodes)


  auc_res= map(disease_to_predict, function(x){
    print(x)

    test_genes= dg_fil %>% filter(nodeB == x) %>% pull(nodeA)

    if(length(test_genes)< 1){
      return(NULL)
    }

    # use distance matrix to remove disease that share a lot of genes)
    disease_to_remove= names(dis_distance[x, dis_distance[x, ]<cutoff])
    disease_to_remove= unique(c(disease_to_remove, x))

    #remove critical diseases from the bipartite mapping:
    dg_2= dg_fil %>%
      filter(!nodeB %in% disease_to_remove)


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
                                       r=0.85)

    g.test= RWRH_PPI_Disease_Results$RWRMH_Multiplex1 %>%
      rename(value= Score,
             name= NodeNames)%>% as_tibble()

    test_res= validate_results_partial(set= test_genes, gene_results = g.test, FPR_cutoff = 0.02)

    # # add a characterization of gene set prediction
    # dg_nodes = unique(dg_2$nodeA)
    # gene_coverage_by_layer=
    #   enframe(test_genes, value = "name", name= "num")%>%
    #   select(-num)%>%
    #   mutate(dg= ifelse(name %in% dg_nodes, 1, 0),
    #          goMF= ifelse(name %in% node.list.genes$goMF, 1, 0),
    #          PPI= ifelse(name %in% node.list.genes$PPI, 1, 0),
    #          Omni= ifelse(name %in% node.list.genes$Omi, 1, 0),
    #          goBP= ifelse(name %in% node.list.genes$goBP, 1, 0)
    #   )%>%
    #   pivot_longer(names_to = "genelayer", values_to = "member", -name)%>%
    #   group_by(genelayer)%>%
    #   summarise(sum= sum(member)) %>%
    #   mutate(percentage= sum/length(test_genes))

    return(test_res)

  })
}



 cutoff_distance= 0.7


#  call for each layer ----------------------------------------------------


#### DD-comorbidity net
randomized_dd= function(edges, disease_sets){

  xdd =randomize_links_by_layer(as_data_frame(edges$disease$heidelberg))
   # we only need to pass different versions of the randomized dd net:

  ppi_net_multiplex= create.multiplex(edges$gene)

  dg= edges$disease_gene %>% filter(nodeB %in% unique(c(xdd$real$from, xdd$real$to)))

  dis_distance= get_Disgenet_overlaps(dg)

  res_dd=
      lapply(xdd, function(x){

        dnet= graph_from_data_frame(x)

        dd_net_multiplex= create.multiplex(list("heidelberg"= dnet))#,"barcelona"= dd_net2)) #to be decided

        aurocs. = get_disease_pred_results(disease_to_predict = disease_sets,
                                 ppi_net_multiplex,
                                 dd_net_multiplex = dd_net_multiplex,
                                 dg_fil = dg,
                                 cutoff = cutoff_distance,
                                 dis_distance = dis_distance

                                )


    })

  return(res_dd)
}

#### D-G layer

randomized_dg =function(edges, disease_sets){

  ppi_net_multiplex= create.multiplex(edges$gene)

  dd_net_multiplex= create.multiplex(list("heidelberg"= edges$disease$heidelberg,
                                          "hpo"= edges$disease$hpo))#,"barcelona"= dd_net2)) #to be decided

  dg= edges$disease_gene %>% filter(nodeB %in% dd_net_multiplex$Pool_of_Nodes)

  dis_distance= get_Disgenet_overlaps(dg)

  xdg= randomize_links_by_layer(dg)

  res_dg= lapply(xdg, function(x){

    x = x %>% filter(nodeA %in% ppi_net_multiplex$Pool_of_Nodes,
                     nodeB %in% dd_net_multiplex$Pool_of_Nodes)

    aurocs. = get_disease_pred_results(disease_to_predict = disease_sets,
                                         ppi_net_multiplex,
                                         dd_net_multiplex = dd_net_multiplex,
                                         dg_fil = x,
                                         cutoff = cutoff_distance,
                                        dis_distance = dis_distance

      )

    })


  return(res_dg)
}

#### G-G layer
randomized_gg= function(edges, disease_sets){

  edges.obj= map(edges$gene,function(x) (randomize_links_by_layer(as_data_frame(x))))

  dd_net_multiplex= create.multiplex(list("heidelberg"= edges$disease$heidelberg,
                                          "hpo"= edges$disease$hpo)
                                     )#,"barcelona"= dd_net2)) #to be decided

  dg= edges$disease_gene %>% filter(nodeB %in% dd_net_multiplex$Pool_of_Nodes)

  dis_distance= get_Disgenet_overlaps(dg)

  # disease_sets=disease_sets[5]
  # layer= names(edges.obj)[1]
  # x= names(edges.obj$goMF)[1]
  res_gg= map(names(edges.obj), function(layer){

    print(layer)

    res_gg= lapply(names(edges.obj[[layer]]), function(x){

      # use the modified edge list to create a single perturbed gene layer:
      gg= edges.obj[[layer]][[x]]%>%
        graph_from_data_frame(., directed = F)

      # use the other unperturbed gene layers
      new.gglist= edges$gene[names(edges$gene) != layer]

      #merge both:
      new.gglist[[layer]]= gg

      ppi_net_multiplex=  create.multiplex(new.gglist)

      aurocs. = get_disease_pred_results(disease_to_predict = disease_sets,
                                         ppi_net_multiplex,
                                         dd_net_multiplex = dd_net_multiplex,
                                         dg_fil = dg_fil,
                                         cutoff = cutoff_distance,
                                         dis_distance = dis_distance
      )

    })


  })

  saveRDS(res_gg, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/randomization_layer_gg.rds")

 return(res_gg)
}

##

  dd_res= randomized_dd(edges, disease_sets)
  saveRDS(dd_res, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/randomization_layer_dd2.rds")

  dg_res= randomized_dg(edges, disease_sets)
  saveRDS(dg_res, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/randomization_layer_dg.rds")


  dg_res2= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/randomization_layer_dg.rds")


# plot dd -----------------------------------------------------------------
library(ggpubr)

plot_feat= function(feat= "pr", dd_res){

  listed.aurocs= map(dd_res, function(x){
    map(x, function(y){
      y[[feat]]
    })%>% unlist()
  })

  dd.df = enframe(listed.aurocs)%>% unnest(value) %>% mutate(disease= rep(disease_sets, 3))

  p.AUROCS.paired=
    ggpaired(dd.df%>% filter(name != "random.degree"), x = "name", y = "value",
             color = "black", line.color = "gray", line.size = 0.1  , linetype =1     )+
    stat_compare_means(paired = TRUE, method = "wilcox.test")


}


p1= plot_feat(feat= "AUROC", dd_res)+ggtitle("AUROC")
p2= plot_feat(feat= "pAUROC", dd_res)+ggtitle("partial AUROC")
p3= plot_feat(feat= "median_rank_stat", dd_res)+ggtitle("median_rank")

dd.df = enframe(listed.aurocs)%>% unnest(value) %>% mutate(disease= rep(disease_sets, 3))

dd_res$random.degree[[1]]$pAUROC.object.corrected$auc
listed.aurocs= map(dd_res, function(x){
  map(x, function(y){
    y$pr$auc.integral
  })%>% unlist()
})

dd.df = enframe(listed.aurocs)%>% unnest(value) %>% mutate(disease= rep(disease_sets, 3))

p4=
  ggpaired(dd.df%>% filter(name != "random.degree"), x = "name", y = "value",
           color = "black", line.color = "gray", line.size = 0.1  , linetype =1     )+
  stat_compare_means(paired = TRUE, method = "wilcox.test")+ggtitle("PR_ROC")

p4

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/hetnet/HFnet_randomization.pdf",
    height= 5, width = 7)

plot_grid(p1,p3,p4, nrow= 1)
dev.off()

