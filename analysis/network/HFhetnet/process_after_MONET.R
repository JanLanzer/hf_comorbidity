## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2021-12-17
##
## Copyright (c) Jan D. Lanzer, 2021
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
##  read MONET results
library(tidyverse)

uw_modules= read_table("data/2022-06-15-154247__M1__result-modules__gene_edges_uw.tsv", col_names = F)
uw_modules= read_table("data/", col_names = F)

uw_modules= read_delim("data/2022-06-17-113941__K1__result-modules__gene_edges_uw.tsv", col_names = F)
length(colnames(uw_modules))
df = unite(data= uw_modules[,3:length(colnames(uw_modules))], col= "geneset",  sep= ",")
x= df[2,]
df.genelists= apply(df, MARGIN = 1, function(x){
  #print(x)
  listed= str_split(x, ",")%>% unlist()
  listed= str_split(listed, "\t")%>% unlist()
  listed= listed[listed != "NA"]
  return(listed)
})

df.genelists[10:20]


table(unlist(map(df.genelists, length)))
hist(unlist(map(df.genelists, length)))


saveRDS(df.genelists, "data/Monet_K1_cluster.rds")



# K1 ------------------------------------------------------------------------------------------

uw_modules= read_delim("data/2022-06-17-113941__K1__result-modules__gene_edges_uw.tsv", col_names = F)

length(colnames(uw_modules))
df = unite(data= uw_modules[,3:length(colnames(uw_modules))], col= "geneset",  sep= ",")

df.genelists= apply(df, MARGIN = 1, function(x){
  print(x)
  listed= str_split(x, ",")%>% unlist()
  listed= listed[listed != "NA"]
  return(listed)
})

df.genelists[10:20]


table(unlist(map(df.genelists, length)))
hist(unlist(map(df.genelists, length)))


saveRDS(df.genelists, "data/Monet_K1_cluster.rds")

uw_modules= read_delim("data/2022-06-17-134633__R1__result-modules__gene_edges_uw.tsv", col_names = F)
length(colnames(uw_modules))
df = unite(data= uw_modules[,3:length(colnames(uw_modules))], col= "geneset",  sep= ",")
x= df[2,]
df.genelists= apply(df, MARGIN = 1, function(x){
  #print(x)
  listed= str_split(x, ",")%>% unlist()
  listed= str_split(listed, "\t")%>% unlist()
  listed= listed[listed != "NA"]
  return(listed)
})

df.genelists[10:20]


table(unlist(map(df.genelists, length)))
hist(unlist(map(df.genelists, length)))


saveRDS(df.genelists, "data/Monet_R1_cluster.rds")

# interpret cluster  --------------------------------------------------------------------------

msigDB= readRDS("/home/jan/R-projects/sc_hfpef/data/prior_knowledge/Genesets_Dec19.rds")


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

msigDB
df.genelists[[1]]

msig_sigs= map(df.genelists[1:5], function(x){
  lapply(names(msigDB[[1]]), function(y){
    #print(x)
    #print(msigDB[[1]][[y]])
    GSE_analysis(x,msigDB[[1]][[y]])

  })%>% do.call(rbind,. )%>% as_tibble()
}) %>% do.call(rbind, .) %>% as_tibble()

names(df.genelists)= paste0("clust.", 1:length(df.genelists))
msig_sigs= map(names(df.genelists[1:50]), function(x){
  GSE_analysis(geneList = df.genelists[[x]],
               Annotation_DB = msigDB[[1]])%>%
    mutate(cluster= x)

    })%>% do.call(rbind,. )%>% as_tibble()


ggplot(msig_sigs, aes(x= gset, y= cluster, fill = -log10(p_value)))+
  geom_tile()+
  coord_flip()
#}) %>% do.call(rbind, .) %>% as_tibble()



