## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-12-15
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------
library(RandomWalkRestartMH)
library(tidyverse)
library(ComplexHeatmap)
library(ggrepel)
library(circlize)
library(cowplot)
library(igraph)


source("analysis/utils/utils_hetnet.R")

edge.list= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds")

#ML.class= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")
comp_feat= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")

feat.vector= comp_feat%>% filter(estimate!= 0) %>% arrange(desc(abs(estimate))) %>%
  pull(PheCode)

comp_feat= comp_feat%>% filter(PheCode  %in% feat.vector[1:100])
#comp_feat= comp_feat%>% filter(PheCode  %in% feat.vector[1:50])
ML.class= comp_feat %>% mutate(hf = ifelse( estimate< 0, "hfpef",
                                            ifelse(estimate>0, "hfref", "ns")
                                            )
                               )

# get gene ranking for hfpef and hfref nodes ------------------------------



ppi_net_multiplex= create.multiplex(edge.list$gene)

dd_net_multiplex= create.multiplex(edge.list$disease)

dg_fil= edge.list$disease_gene

PPI_Disease_Net <- create.multiplexHet(ppi_net_multiplex,
                                       dd_net_multiplex,
                                       as.data.frame(dg_fil[,c(1,2,3)]), "disease")

PPIHetTranMatrix <- compute.transition.matrix(PPI_Disease_Net)



seed.hfpef= ML.class %>% filter(hf=="hfpef")%>% pull(PheCode)
#tau.hfpef= ML.class %>% filter(hf=="hfpef")%>% pull(estimate)%>% abs()
#tau.hfpef2= tau.hfpef/sum(tau.hfpef)
seed.hfref= ML.class %>% filter(hf=="hfref")%>% pull(PheCode)
#tau.hfref= ML.class %>% filter(hf=="hfref")%>% pull(estimate)
#tau.hfref2= tau.hfref/sum(tau.hfref)

g.hfpef <-
  Random.Walk.Restart.MultiplexHet(x= PPIHetTranMatrix,
                                   MultiplexHet_Object = PPI_Disease_Net,
                                   Multiplex1_Seeds= c(),
                                   Multiplex2_Seeds = seed.hfpef,
                                   #tau2= tau.hfpef2,
                                   r=0.8)

g.hfref <-
  Random.Walk.Restart.MultiplexHet(x= PPIHetTranMatrix,
                                   MultiplexHet_Object = PPI_Disease_Net,
                                   Multiplex1_Seeds= c(),
                                   Multiplex2_Seeds = seed.hfref,
                                   #tau2= tau.hfref2,
                                   r=0.8)



g.hfpef= g.hfpef$RWRMH_Multiplex1 %>%
  rename(value= Score,
         name= NodeNames)%>% as_tibble()

g.hfref= g.hfref$RWRMH_Multiplex1 %>%
  rename(value= Score,
         name= NodeNames)%>% as_tibble()

#
# PLOT RANKS --------------------------------------------------------------

g.pef= g.hfpef %>% mutate(rank= rank(desc(value)))
g.ref= g.hfref %>% mutate(rank= rank(desc(value)))

# add a simple way to prioritize genes by multiplying RW probability with rank-difference and rank new results:
g.ranks= full_join(g.pef, g.ref, by= "name") %>%
  mutate(rank.diff= rank.y-rank.x) %>%
  rename(rank.hfpef = rank.x,
         rank.hfref= rank.y,
         RW.value.hfpef= value.x,
         RW.value.hfref= value.y) %>%
  mutate(hfref.prio= RW.value.hfref * -rank.diff,
         hfpef.prio= RW.value.hfpef * rank.diff)

#get gene names
gene.names= read.csv(url("https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"),
                     sep= "\t")%>%
  as_tibble()%>%
  rename(gene= Approved.symbol,
         gene.name= Approved.name)%>%
  select(gene, gene.name)

g.ranks= g.ranks%>%
  rename(gene= name)%>%
  left_join(gene.names)%>%
  mutate(rank.hfpef.prio= rank(desc(hfpef.prio)),
         rank.hfref.prio= rank(desc(hfref.prio)))%>%
  select(gene, gene.name, everything())

# save results

saveRDS(g.ranks, file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks.rds")
g.ranks= readRDS( file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks.rds")

saveRDS(list("hfpef"= g.hfpef,
             "hfref"= g.hfref), file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks2.rds")

g.list= readRDS( file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks2.rds")

# as .csv
g.ranks %>% write.csv(., file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/predicted_HF_genes.csv")
g.ranks %>% arrange(desc(hfpef.prio))%>% filter(grepl("PPA", gene))

# use HF genesets for comparison -------------------------------------------------------------------

## use HF genesets
sets= load_validation_genes(disgenet_value= 0.29)
names(sets)= c("DisGeNET", "PheWAS", "Kegg_DCM", "ReHeaT", "r2", "Cardiomyopathy", "Top_single_variants", "Top_common_variants")
sets= sets[-5]
val.set= unique(unlist(sets))
length(unique(unlist(sets)))

#calc overlap
df= sapply(sets, function(x){
  sapply(sets, function(y){
    #length(intersect(x,y))/length(union(x,y))
    x = prop.table(table(x %in% y))[1]
    as.numeric(x)
  })
})
rownames(df) = colnames(df)
diag(df)= NA
df= 1-df

hmap.lap= Heatmap(df,
        row_names_side = "left",
        name= "% geneset (column)\nin geneset (row)",
        cluster_rows = F,
        cluster_columns = F)
pdf(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/hetnet/HF_geneset_overlap.pdf",
    width= 5,
    height= 4)
hmap.lap
dev.off()



pef= sapply(sets, function(x){
  res = validate_results2(set = x, gene_results = g.list[[1]])
  c("PR_AUC"= res$pr$auc.integral,
    "AUROC" =res$roc$auc)
})
ref= sapply(sets, function(x){
  res = validate_results2(x, g.list[[2]])
  c("PR_AUC"= res$pr$auc.integral,
    "AUROC" =res$roc$auc)
})

#ref$HF= rep("HFrEF",2)



###
roc.df%>%
  ggplot(aes(x= value , y= as.factor(validation_genes)))+
  geom_point()+
  geom_boxplot()+
  #theme(axis.text. = element_blank())+
  coord_flip()

###
df= rbind(pef, ref)
#rownames(df) = c("HFpEF_PR_AUC", "HFpEF_AUROC",  "HFrEF_PR_AUC", "HFrEF_AUROC" )
rownames(df) = c("HFpEF", "HFpEF",  "HFrEF", "HFrEF" )

col_fun = colorRamp2(c(0, 0.1, 0.5, 1), c("white", "blue", "red", "darkred"))
df%>% t() %>% Heatmap(name = "value")

map_genes=Heatmap(t(df),
                  name = "AUC",
                  #col = col_fun,
                  top_annotation = HeatmapAnnotation(foo = anno_block(labels = c("AUC-PR", "AUROC") , gp = gpar(fill = 2))),
                  column_km = 2,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.3f", t(df)[i, j]), x, y, gp = gpar(fontsize = 10))},
                  cluster_columns = F,
                  show_row_dend = F
)

pdf(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/HF_geneset_recovery.pdf",
    width= 4.1,
    height= 4)
print(map_genes)
dev.off()




#do wilcox paired
df= rbind(pef, ref)
rownames(df) = c("HFpEF_PR_AUC", "HFpEF_AUROC",  "HFrEF_PR_AUC", "HFrEF_AUROC" )

df2= t(df)%>%
  as.data.frame()%>%
  rownames_to_column("set") %>% as_tibble()%>%
  pivot_longer(c(HFpEF_AUROC , HFrEF_AUROC), names_to = "name", values_to = "AUC")


library(ggpubr)
p.AUROCS.paired=
  ggpaired(df2,
           x= "name", id= "set", y = "AUC",
           color = "black", line.color = "gray", line.size = 0.5  ,
           linetype =1  ,point.size = 3
            )+
  geom_point(aes(color= set))+
  stat_compare_means(paired = TRUE, method = "wilcox.test")+
  labs(x= "Condition",
       y= "AUROC")

p.AUROCS.paired

## the AUROC might be repeated for the hfpef prio vector

do_auc_on_vec= function(gg, col, sets){
  gg= gg %>% rename("value" := !!(col))

  res= sapply(sets, function(x){
    res = validate_results2(x, gg)
    c("PR_AUC"= res$pr$auc.integral,
      "AUROC" =res$roc$auc)
  })

}

pef= do_auc_on_vec(gg,"hfpef.prio", sets)
ref= do_auc_on_vec(gg,"hfref.prio", sets)



# add randomized comorbidity set prediction to calculate z-AUROC ----------

#we use teh same hfhetnet representation:

seed.hfpef= ML.class %>% filter(hf=="hfpef")%>% pull(PheCode)

seed.hfref= ML.class %>% filter(hf=="hfref")%>% pull(PheCode)

sets= load_validation_genes(disgenet_value= 0.29)
names(sets)= c("DisGeNET", "PheWAS", "Kegg_DCM", "ReHeaT", "r2", "Cardiomyopathy", "Top_single_variants", "Top_common_variants")
sets= sets[-5]


set.seed(5)
get_AUROC_dist= function(seed1= seed.hfpef,
                         nperm= 10,
                         sets){


  df1= lapply(1:nperm, function(p){

    diseases= PPI_Disease_Net$Multiplex2$Pool_of_Nodes
    diseases= diseases[!diseases %in% c(seed1)]

    seeds= sample(diseases, length(seed1), replace = F)

    g.rank <-
      Random.Walk.Restart.MultiplexHet(x= PPIHetTranMatrix,
                                       MultiplexHet_Object = PPI_Disease_Net,
                                       Multiplex1_Seeds= c(),
                                       Multiplex2_Seeds = seeds,
                                       #tau2= tau.hfpef2,
                                       r=0.8)

    g.rank= g.rank$RWRMH_Multiplex1 %>%
      rename(value= Score,
             name= NodeNames)%>% as_tibble()

    FIND= sapply(sets, function(set){

      res = validate_results2(set = set,
                              gene_results = g.rank)

      c("PR_AUC"= res$pr$auc.integral,
        "AUROC" =res$roc$auc)

    } )
  })

  #names(df1)= paste0("rep.", 1:nperm)
  do.call(rbind, df1)%>%
    as.data.frame()%>%
    rownames_to_column("metric")
}

r= get_AUROC_dist(seed1 = seed.hfpef,
                  sets= sets,
                  nperm= 1000)

r2= get_AUROC_dist(seed1 = seed.hfref,
                  sets= sets,
                  nperm= 1000)


calc_z= function(distro, real_score){
  real_score = real_score %>%
    as.data.frame() %>%
    rownames_to_column("metric")

  distro= as_tibble(distro)%>%
    mutate(metric= str_replace_all(metric, ".[0-9]", ""))

  df= sapply(unique(real_score$metric), function(x){
    sapply(colnames(real_score)[-1], function(y){

      X= distro%>% filter(metric == x) %>% pull(!!as.symbol(y))
      Y= real_score %>% filter(metric == x) %>% pull(!!as.symbol(y))

      Z= (Y- mean(X)) / sd(X)

      })
  })
  df %>% as.data.frame()%>% rownames_to_column("set")

}

p.map=rbind(calc_z(r, pef)%>% mutate("comor"= "HFpEF"),
            calc_z(r2, ref)%>% mutate("comor"= "HFrEF")
            )%>%
  pivot_longer(cols = c(PR_AUC, AUROC) , names_to= "metric", values_to = "value")


p.zs= p.map%>% ggplot(., aes(x= set, y= value, col = comor))+
  geom_point(size = 2)+
  geom_hline(yintercept = c(1.65, -1.65, -1.96,1.96), col = "darkgrey")+
  coord_flip()+
  facet_grid(rows= vars(metric))+
  labs(x= "HF gene set",
       y= "z-score",
       color= "comorbidity\nprofile")+
  theme_bw()+
  theme(axis.text = element_text(size= 11, color= "black"))+
  scale_y_continuous(breaks= c(1.65, -1.65, 0,-1.96,1.96))
p.zs

pdf(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/hetnet/HF_geneset_zscores.pdf",
    width= 5.5,
    height= 3.5)
p.zs
dev.off()

Heatmap(p.map)


g.ref= g.ranks %>% select(gene, hfref.prio)%>% arrange(hfref.prio)%>%
  rename(name= gene,
         value= hfref.prio)%>%
  mutate(value= abs(value))
g.pef= g.ranks %>% select(gene, hfpef.prio)%>% arrange(hfpef.prio)%>%
  rename(name= gene,
         value= hfpef.prio)%>%
  mutate(value= abs(value))

pef= sapply(sets, function(x){
  res = validate_results2(set = x, gene_results = g.pef)
  c("PR_AUC"= res$pr$auc.integral,
    "AUROC" =res$roc$auc)
})
ref= sapply(sets, function(x){
  res = validate_results2(x, g.ref)
  c("PR_AUC"= res$pr$auc.integral,
    "AUROC" =res$roc$auc)
})
