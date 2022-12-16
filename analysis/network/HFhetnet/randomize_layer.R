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
library(qdapTools)

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
cutoff_distance= 0.5
# new ---------------------------------------------------------------------

get_disease_pred_results= function(disease_to_predict,
                                   ppi_net_multiplex,
                                  dd_net_multiplex,
                                  dg_fil,
                                  cutoff= 0.5,
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
      dplyr::rename(value= Score,
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






#  call for each layer ----------------------------------------------------


#### DD-comorbidity net
randomized_dd= function(edges, disease_sets, nperm= 10){

  ppi_net_multiplex= create.multiplex(edges$gene)

  dg= edges$disease_gene %>% filter(nodeB %in% V(edges$disease$heidelberg)$name)

  dis_distance= get_Disgenet_overlaps(dg)
  disease_sets= rownames(dis_distance)

  perm_res= map(seq(1:nperm), function(perm){

    set.seed(perm)
    print(perm)
    xdd =randomize_links_by_layer(igraph::as_data_frame(edges$disease$heidelberg), rewireprob= 1)
    #xdd= xdd[c("random.degree", "random.complete")]
    xdd= xdd["random.complete"]
    res_dd= lapply(xdd, function(x){

          dnet= igraph::graph_from_data_frame(x)

          dd_net_multiplex= create.multiplex(list("heidelberg"= dnet))#,"barcelona"= dd_net2)) #to be decided

          aurocs. = get_disease_pred_results(disease_to_predict = disease_sets,
                                  ppi_net_multiplex,
                                   dd_net_multiplex = dd_net_multiplex,
                                   dg_fil = dg,
                                   cutoff = cutoff_distance,
                                   dis_distance = dis_distance

                                  )
          return(aurocs.)
          })
    return(res_dd)
  })
  return(perm_res)

   # we only need to pass different versions of the randomized dd net:



}

randomized_dd2= function(edges){

  ppi_net_multiplex= create.multiplex(edges$gene)

  dg= edges$disease_gene %>% filter(nodeB %in% unique(V(edges$disease$heidelberg)$name))

  dis_distance= get_Disgenet_overlaps(dg)

  #disease_sets= rownames(dis_distance)

  perm_res= map(seq(0, 1, 0.2), function(perm){
    set.seed(perm)
    print(perm)

    xdd =randomize_links_by_layer(igraph::as_data_frame(edges$disease$heidelberg), rewireprob= perm)

    dnet= igraph::graph_from_data_frame(xdd$random.complete)

    dd_net_multiplex= create.multiplex(list("heidelberg"= dnet))#,"barcelona"= dd_net2)) #to be decided

    aurocs. = get_disease_pred_results(disease_to_predict = disease_sets,
                                         ppi_net_multiplex,
                                         dd_net_multiplex = dd_net_multiplex,
                                         dg_fil = dg,
                                         cutoff = cutoff_distance,
                                         dis_distance = dis_distance
                                       )
    names(aurocs.)= disease_sets
    return(aurocs.)
    })

  names(perm_res)= seq(0, 1, 0.2)
  return(perm_res)
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

##call
r.res= randomized_dd2(edges)

saveRDS(r.res, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/randomization_layer_dd3.rds")

#dd_res= randomized_dd(edges,disease_sets =  disease_sets[1:2],nperm= 1)
  saveRDS(dd_res, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/randomization_layer_dd3.rds")

  dg_res= randomized_dg(edges, disease_sets)
  saveRDS(dg_res, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/randomization_layer_dg.rds")


  dg_res2= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/randomization_layer_dg.rds")


# plot dd -----------------------------------------------------------------

library(ggpubr)

plot_feat= function(feat= "pr", r.res){

  if(feat=="pr"){
    listed.aurocs= map(r.res, function(x){
      map(x, function(y){
        y$pr$auc.integral
      })%>% unlist()
    })
  }else{
    listed.aurocs= map(r.res, function(x){
      map(x, function(y){
        y[[feat]]
      })%>% unlist()
    })
  }

  names(r.res[[1]])

  df= map(names(listed.aurocs), function(x){
    print(x)
  enframe(listed.aurocs[[x]]) %>%
       mutate(res= x,
              feat= feat)

    })

  dd.df= do.call(rbind, df)
  #
  # ggplot(dd.df, aes(x= res, y= value))+
  #   geom_boxplot()
  #
  # p.AUROCS.paired=
  #   ggpaired(dd.df%>% filter(name != "random.degree"), x = "name", y = "value",
  #            color = "black", line.color = "gray", line.size = 0.1  , linetype =1     )+
  #   stat_compare_means(paired = TRUE, method = "wilcox.test")


}

p1= plot_feat(feat= "AUROC", r.res)
#p2= plot_feat(feat= "pAUROC", r.res)
p2= plot_feat(feat= "pr", r.res)
p3= plot_feat(feat= "median_rank_stat", r.res)
df.p= rbind(p1,p2,p3)

saveRDS(df.p,"T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/new_df.rds" )

## pair data with real results:
auc_res= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/aurocs_multilayer_disease_pred_heart_genenet2022.rds")

plot_feat2= function(feat= "pr", r.res2){

  if(feat=="pr"){
    listed.aurocs= map(r.res2, function(x){
      x$test_results$pr$auc.integral
    })
  }else{
    listed.aurocs= map(r.res2, function(x){
      x$test_results[[feat]]
    })
  }

  enframe(listed.aurocs) %>%
      unnest(value)%>%
      mutate(res= "HFnet+HPOnet",
             feat= feat)

  }

p1= plot_feat2(feat= "AUROC", r.res2 = auc_res)
p2= plot_feat2(feat= "pr", auc_res)
p3= plot_feat2(feat= "median_rank_stat", auc_res)

df.p = rbind(df.p, p1, p2,p3)

saveRDS(df.p,"T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/gene_pred_res.rds" )

df.p= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/gene_pred_res.rds" )
dis= unique(df.p %>% filter(res== "HFnet+HPOnet")%>% pull(name))
unique(df.p$res)
df.p$res= str_replace_all("0.", "rewire HFnet 0.", string = df.p$res)
#df.p$res= str_replace_all("0", "HFnet", string = df.p$res)
df.p = df.p %>%
  mutate(feat = ifelse(feat=="pr", "PR-AUC", feat),
         res= ifelse(res== "0", "HFnet", res),
         res= ifelse(res== "1", "rewired HFnet", res))%>%
  filter(res %in% c("HFnet+HPOnet", "HFnet", "rewired HFnet"))%>%
  mutate(res= factor(res, levels= c("HFnet+HPOnet", "HFnet", "rewired HFnet")))


p.all= ggplot(df.p%>% filter(name %in% dis),
              aes(x= res, y= value))+
  geom_jitter(color= "darkgrey")+
  geom_boxplot(width= 0.3,outlier.shape = NA)+
  facet_grid(cols= vars(feat), space = "free", scales= "free")+
  theme_bw()+
  #scale_y_continuous(limits= c(0,0.3))+
  labs(x= "rewire probability")+
  theme(axis.text = element_text(angle= 45, hjust= 1))

p.all

x= "PR-AUC"

p1= map(unique(df.p$feat), function(x){

  df= df.p%>% filter(name %in% dis & feat == x)
  if(x== "AUROC"){
    li= quantile(df$value, c(0.01, 1))
    #li= c(0.31, 1.2)
  }else{
    li= quantile(df$value, c(0,0.99))
    #li= c(0, )
  }

  my_comparisons <- list( c("HFnet", "HFnet+HPOnet"), c("HFnet", "rewired HFnet"))

  p.all= ggplot(df, aes(x= res, y= value))+
    geom_jitter(aes(color= res))+
    geom_boxplot(width= 0.3,outlier.shape = NA)+
    theme_bw()+
    scale_color_manual(values= c(cols.nice))+
    scale_y_continuous(limits= li)+
    labs(x= "",
         y= x)+
    theme(legend.position = "none",
          axis.text.x = element_text(angle= 60, hjust= 1, size= 11))
  unify_axis(p.all)
})

cowplot::plot_grid(plotlist = p1, ncol = 3)

df%>% group_by(res)%>% summarise(median(value))

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/hetnet/HFnet_randomization2022.pdf",
    width= 7,
    height= 5)

cowplot::plot_grid(plotlist = p1, ncol = 3)

dev.off()

p.AUROCS.paired=df.p%>%
  filter(res %in% c("HFnet", "rewired HFnet"),
         feat == "AUROC",
         name %in% dis)%>%
   ggpaired(.,x = "res", y = "value",
            color = "black", line.color = "gray", line.size = 0.1  , linetype =1     )+
   stat_compare_means(paired = TRUE, method = "wilcox.test")

stat.test <- df.p %>%
  group_by(feat) %>%
  wilcox_test(value ~ res, p.adjust.method = "BH", paired = T)
# Remove unnecessary columns and display the outputs
stat.test




p.AUROCS.paired

df.p%>%
  filter(res %in% c("0", "0.2"),
         feat == "AUROC")%>%
        # !name %in% xs)%>%
  ggplot( aes(x= res, y= value))+
  geom_jitter(color= "darkgrey")+
  geom_boxplot(width= 0.3)


dd.df = enframe(listed.aurocs)%>% unnest(value) %>% mutate(disease= rep(disease_sets, 3))


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

