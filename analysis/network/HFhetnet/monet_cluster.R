## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-06-17
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## interpret monet cluster and measure distance to hf nodes
## ---------------------------

library(tidyverse)
library(RandomWalkRestartMH)

M1= readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/Monet_M1_cluster.rds")
K1= readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/Monet_K1_cluster.rds")
R1= readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/Monet_R1_cluster.rds")
source("analysis/utils/utils_hetnet.R")
source("analysis/utils/utils_network.R")

edge.list= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds")
phedic= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/icd10_phewas_dictionary.rds")%>% distinct(PheCode, Phenotype)


# create hetnet obj
ppi_net_multiplex= create.multiplex(edge.list$gene)

dd_net_multiplex= create.multiplex(edge.list$disease)

dg_fil= edge.list$disease_gene

PPI_Disease_Net <- create.multiplexHet(ppi_net_multiplex,
                                       dd_net_multiplex,
                                       as.data.frame(dg_fil[,c(1,2,3)]), "disease")

PPIHetTranMatrix <- compute.transition.matrix(PPI_Disease_Net)


# loop over gene clusters and measure distance to hfpef and hfref---------------
M1= K1

names(M1)= paste0("cluster.", 1:length(M1))
#take cluster with 4 genes and more
M1.f= M1[lengths(M1)>1]

unlist(M1.f)[!unlist(M1.f) %in% PPI_Disease_Net$Multiplex1$Pool_of_Nodes]
M1.f$cluster.393

res. = map(M1.f, function(x){
  x= x[x%in%PPI_Disease_Net$Multiplex1$Pool_of_Nodes]
  g.hfpef <-
    Random.Walk.Restart.MultiplexHet(x= PPIHetTranMatrix,
                                     MultiplexHet_Object = PPI_Disease_Net,
                                     Multiplex1_Seeds= x,
                                     Multiplex2_Seeds = c(),
                                     r=0.8)


})

#rank transform disease rankings for every gene cluster
df= map(res., function(x){
  x$RWRMH_Multiplex2 %>% mutate(rank= rank(desc(Score), "average"))%>% as_tibble()%>%
    rename(PheCode=NodeNames)%>%
    left_join(phedic)
})

#get hfpef and hfref cluster ranks
df2= sapply(df, function(x){
  c("hfpef"= x %>% filter(PheCode== "hfpef")%>% pull(rank),
    "hfref"= x %>% filter(PheCode== "hfref")%>% pull(rank)
  )
})%>% t()%>%
  as.data.frame() %>%
  rownames_to_column("cluster")%>%as_tibble()

df2%>%
  pivot_longer(names_to = "HFtype", -cluster)%>%
  mutate(cluster= factor(cluster, levels= names(M1.f)))%>%
  ggplot(aes(x= cluster, y= value, color= HFtype, group = HFtype))+
  geom_point()+
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1 ,vjust= 0.5))

df2%>% mutate(rank.diff =hfpef-hfref)%>%
  mutate(cluster= factor(cluster, levels= names(M1.f)))%>%
  ggplot(aes(x= cluster, y= rank.diff))+
  geom_point()+
  theme_bw()+
  geom_hline(yintercept = 0)+
  theme(axis.text.x = element_text(angle= 90, hjust= 1 ,vjust= 0.5))

M1.f$cluster.19

clust.hfpef= df2%>% mutate(rank.diff =hfpef-hfref)%>% arrange((rank.diff))%>% slice(1:25)%>% pull(cluster)
M1.f[clust.hfpef]

df2%>%
  pivot_longer(names_to = "HFtype", -cluster)%>%
  filter(cluster %in% clust.hfpef)%>%
  mutate(cluster= factor(cluster, levels= (clust.hfpef)))%>%
  ggplot(aes(x= cluster,  y= value, color= HFtype, group = HFtype))+
  geom_point()+
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1 ,vjust= 0.5))



clust.hfpef2= df2%>% arrange(hfpef)%>% slice(1:25)%>% pull(cluster)
M1.f[clust.hfpef2]

df2%>%
  pivot_longer(names_to = "HFtype", -cluster)%>%
  filter(cluster %in% clust.hfpef2)%>%
  mutate(cluster= factor(cluster, levels= (clust.hfpef2)))%>%
  ggplot(aes(x= cluster,  y= value, color= HFtype, group = HFtype))+
  geom_point()+
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1 ,vjust= 0.5))



# cluster interpret: ------------------------------------------------------

msigDB= readRDS(file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/databases/Genesets_Dec19.rds")

##ORA func

GSE_analysis = function(geneList,Annotation_DB){
  library(dplyr)
  library(tidyr)
  library(tibble)

  geneList = geneList[geneList %in% unique(unlist(Annotation_DB))]

  ResultsDF = matrix(0,nrow = length(Annotation_DB),ncol = 5)
  rownames(ResultsDF) = names(Annotation_DB)
  colnames(ResultsDF) = c("GenesInPathway","GenesInList","GeneNames","p_value","corr_p_value")

  DB_genecontent = length(unique(unlist(Annotation_DB)))

  GenesDB = DB_genecontent
  SelectedGenes = length(geneList)

  for(gset in rownames(ResultsDF)){
    GP = length(((Annotation_DB[[gset]])))
    GL = length(intersect(Annotation_DB[[gset]],geneList))

    ResultsDF[gset,"GenesInList"] = GL
    ResultsDF[gset,"GenesInPathway"] = GP
    ResultsDF[gset,"GeneNames"] = paste(intersect(Annotation_DB[[gset]],geneList),collapse = ",")
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)
  }

  ResultsDF[,"corr_p_value"] = p.adjust(ResultsDF[,"p_value"],method = "BH")
  ResultsDF = data.frame(ResultsDF,stringsAsFactors = F)
  ResultsDF = ResultsDF[order(ResultsDF[,"p_value"]),]

  ResultsDF = ResultsDF %>%
    rownames_to_column("gset") %>%
    mutate_at(c("GenesInPathway","GenesInList",
                "p_value","corr_p_value"),
              as.numeric) %>%
    dplyr::arrange(corr_p_value,GenesInList)

  return(ResultsDF)

}


df.genelists= M1.f[clust.hfpef]

#msig_sigs= map(df.genelists, function(x){
  lapply(names(msigDB[[1]]), function(y){
    #print(x)
    #print(msigDB[[1]][[y]])
    GSE_analysis(x,msigDB[[1]][[y]])

  })%>% do.call(rbind,. )%>% as_tibble()
}) %>% do.call(rbind, .) %>% as_tibble()

#names(df.genelists)= paste0("clust.", 1:length(df.genelists))
msig_sigs= map(names(df.genelists), function(x){
  GSE_analysis(geneList = df.genelists[[x]],
               Annotation_DB = msigDB[[5]])%>%
    mutate(cluster= x)

})%>% do.call(rbind,. )%>% as_tibble()

gsets= msig_sigs%>% group_by(gset)%>% filter(corr_p_value<0.01)%>% pull(gset)

msig_sigs%>% filter(cluster == "cluster.11")%>% arrange(corr_p_value)

msig_sigs %>% group_by(cluster) %>% top_n(., n = 2, wt= -corr_p_value)

ggplot(msig_sigs%>% filter(gset %in% gsets),
       aes(x= gset, y= cluster, fill = -log10(p_value)))+
  geom_tile()+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1 ,vjust= 0.5))
#}) %>% do.call(rbind, .) %>% as_tibble()

msig_sigs%>% filter(corr_p_value<0.05)



