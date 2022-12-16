## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-06-21
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## use GWAS results for gene prioritization
## ---------------------------

library(decoupleR)
library(fgsea)
library(readr)
library(tidyverse)

## diverse GWAS
gwas1= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/gwas_rebecca/Aung_et_al_for_Pascal.sum.genescores.txt.gz")
gwas2= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/gwas_rebecca/HERMES_gwas_for_pascal.sum.genescores.txt.gz")

gwas3= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/hfpef_snps_with_nearest_gene_unique.txt")
gwas4= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/hfref_snps_with_nearest_gene_unique.txt")

gwas.pasc.pef= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/hfpef_rs_snp_p.sum.genescores.txt")
gwas.pasc.ref= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/hfref_rs_snp_p.sum.genescores.txt")

gwas.pasc.pef2= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/hfpef_rs_snp_p_large_window.sum.genescores.txt")
gwas.pasc.ref2= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/hfref_rs_snp_p_large_window.sum.genescores.txt")

gwas.lvedvi= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/MRI_lvedvi_snps.sum.genescores.txt")
gwas.lvef= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/MRI_lvef_snps.sum.genescores.txt")
gwas.lveSvi= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/MRI_lvesvi_snps.sum.genescores.txt")
gwas.svi= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/MRI_svi_snps.sum.genescores.txt")

gwas.hits= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/41467_2020_15823_MOESM10_ESM(1).csv", sep= ";")
## gene hits

gg= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks.rds")

## comorbidity gwas:
gwas.bmi= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/BMI_snps_p.sum.genescores.txt")
gwas.nicm= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/NICM_snps_p.sum.genescores.txt")
gwas.dm= read_tsv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/HF_subtype_gwas/T2D_EUR_snps_p.sum.genescores.txt")


##### first mapping (not pascal, but minimal p)

get_min_p = function(x){
  x%>%
    group_by(V8)%>%
    summarise(m.p= min(Pvalue))%>%
    arrange(m.p) %>%
    mutate(logp= -log10(m.p))
}

gwas.hfpef= as_tibble(gwas3) %>% arrange(Pvalue)
gwas.hfref= as_tibble(gwas4) %>% arrange(Pvalue)

gwas.hfpef2 = get_min_p(gwas.hfpef)%>% rename(gene_symbol = V8)
gwas.hfref2 = get_min_p(gwas.hfref)%>% rename(gene_symbol = V8)

## network prediciton gene sets
gc.hfpef= gg %>% arrange(desc(hfpef.prio)) %>% slice(1:50) %>% pull(name)
gc.hfref= gg %>% arrange(desc(hfref.prio)) %>% slice(1:50) %>% pull(name)

gsets= list("hfpef"= gc.hfpef,
            "hfref"= gc.hfref)

###
perform_fgsea= function(gwas,
                        gwas.name,
                        gsets,
                        ...){
  set.seed(3)
  gwas = gwas %>% arrange(desc(logp))
  stats= gwas$logp
  names(stats)= gwas$gene_symbol


  gsea.res=fgseaMultilevel(pathways = gsets,
                       stats = stats,
                      # nperm = 1000,
                      scoreType = "pos"
  )

  #gsea.res %>% mutate(gwas= gwas.name)

}


pef.res= perform_fgsea(gwas.hfpef2, "hfpef",gsets = gsets) %>% mutate(gwas= "hfpef")
ref.res= perform_fgsea(gwas.hfref2, "hfref",gsets = gsets, nperm = 2000) %>% mutate(gwas= "hfref")

df= rbind(pef.res, ref.res)

ggplot(df, aes(x= gwas, y= pathway, fill = -log10(pval)))+
  geom_tile()+
  scale_fill_gradient(low= "white", high= "red")

saveRDS(df,"T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/gwas.enrich.results.rds" )



# test pascal maps --------------------------------------------------------


pef.res= perform_fgsea(gwas.pef, "hfpef", gsets = gsets) %>% mutate(gwas= "hfpef")
ref.res= perform_fgsea(gwas.ref, "hfref",gsets = gsets, nperm = 1000) %>% mutate(gwas= "hfref")

gwas.pef= gwas.pasc.pef2%>% mutate(logp= -log10(pvalue))
gwas.ref= gwas.pasc.ref2%>% mutate(logp= -log10(pvalue))
pef.res= perform_fgsea(gwas.pef, "hfpef", gsets = gsets) %>% mutate(gwas= "hfpef")
ref.res= perform_fgsea(gwas.ref, "hfref",gsets = gsets, nperm = 10000) %>% mutate(gwas= "hfref")

pef.res$leadingEdge
ref.res$leadingEdge

# test gene cut offs ------------------------------------------------------
test_cutoffs= function(sqc= seq(25,500, 25),
                       gwas1,
                       gwas2,
                       name1="hfpef",
                       name2="hfref",
                       gg){

  l= map(sqc, function(x){

    #gc.hfpef= gg %>% arrange(desc(hfpef.prio)) %>% slice(1:x) %>% pull(name)
    #gc.hfref= gg %>% arrange(desc(hfref.prio)) %>% slice(1:x) %>% pull(name)
    gc.hfpef= gg %>% arrange(desc(RW.value.hfpef)) %>% slice(1:x) %>% pull(gene)
    gc.hfref= gg %>% arrange(desc(RW.value.hfref)) %>% slice(1:x) %>% pull(gene)

    gsets= list("hfpef"= gc.hfpef,
                "hfref"= gc.hfref)

    pef.res= perform_fgsea(gwas1 %>% filter(Status != "DAVIES_FAIL_FAREBROTHER_FAIL"),
                           name1,gsets = gsets, nperm= 2000) %>%
      mutate(gwas= name1, cutoff= x)
    ref.res= perform_fgsea(gwas2%>%filter(Status != "DAVIES_FAIL_FAREBROTHER_FAIL"),
                           name2, gsets=gsets, nperm = 2000) %>%
      mutate(gwas= name2, cutoff= x)

    df= rbind(pef.res, ref.res)
  })

  do.call(rbind, l) %>% as_tibble() %>% mutate(log10p= -log10(pval))%>%
    ggplot(., aes(x= cutoff, y= log10p, col= pathway))+
    facet_grid(rows= vars(gwas))+
    geom_point()+
    geom_path()+
    geom_hline(yintercept = -log10(0.05))+
    geom_hline(yintercept = -log10(0.1), lty= 2)



}


p1= test_cutoffs(gwas.hfpef =gwas.pasc.pef2%>% mutate(logp= -log10(pvalue)),
             gwas.hfref = gwas.pasc.ref2%>% mutate(logp= -log10(pvalue)),
             gg= gg)



p2= test_cutoffs(gwas.hfpef =gwas.bmi%>% mutate(logp= -log10(pvalue)),
             gwas.hfref = gwas.dm%>% mutate(logp= -log10(pvalue)),
             gg= gg)

p3= test_cutoffs(gwas.hfpef =gwas.nicm%>% mutate(logp= -log10(pvalue)),
                 gwas.hfref = gwas.dm%>% mutate(logp= -log10(pvalue)),
                 gg= gg)
p4= test_cutoffs(sqc= seq(10, 200, 10) ,
                 name1= "LVEDVI",
                 name2 = "LVEF",
                 gwas1 = gwas.lvedvi%>% mutate(logp= -log10(pvalue))%>%
                   filter(!Status %in% c("DAVIES_FAIL_FAREBROTHER_FAIL", "NOT_RUN")),
                 gwas2 = gwas.lvef%>% mutate(logp= -log10(pvalue))%>%
                   filter(!Status %in% c("DAVIES_FAIL_FAREBROTHER_FAIL", "NOT_RUN")),
                 gg= gg)

p5= test_cutoffs(sqc= seq(10, 200, 10) ,
                 name1= "svi",
                 name2 = "lvesvi",
                 gwas1 = gwas.svi%>% mutate(logp= -log10(pvalue))%>%
                   filter(!Status %in% c("DAVIES_FAIL_FAREBROTHER_FAIL", "NOT_RUN")),
                 gwas2 = gwas.lveSvi%>% mutate(logp= -log10(pvalue))%>%
                   filter(!Status %in% c("DAVIES_FAIL_FAREBROTHER_FAIL", "NOT_RUN")),
                 gg= gg)



cowplot::plot_grid(p4,p5)
g.ranks %>%
  ggplot(., aes(x= reorder(name, -hfpef.prio), y= hfpef.prio))+
  geom_point()+
  geom_hline()



g.ranks %>%
  filter(hfpef.prio > 0.001)
  ggplot(., aes(x= reorder(name, -hfref.prio), y= hfref.prio))+
  geom_point()


# get candidates by overlap -----------------------------------------------

  gwas.bmi = gwas.bmi%>% mutate(logp= -log10(pvalue))
  gwas= gwas.bmi %>% distinct(logp, gene_symbol)
  set.seed(3)
  gwas = gwas %>% arrange(desc(logp))
  stats= gwas$logp
  names(stats)= gwas$gene_symbol

  stats



  alpha= 0.0001



  mat= gg%>% select(gene, RW.value.hfpef, RW.value.hfref, hfpef.prio, hfref.prio)%>% as.data.frame()
  set= mat %>% filter(duplicated(gene)) %>% pull(gene)
  mat= as.matrix(mat %>% filter(!gene %in% set)%>%column_to_rownames("gene"))


  f.gwas= function(alpha=0.01, gwas){
    gwas%>% filter(pvalue<alpha)%>% pull(gene_symbol)

  }

  l= list(# "HFrEF"= gwas.pasc.ref2,
           #"HFREF2"= gwas.pasc.ref,
           #"HFPEF2"= gwas.pasc.pef,
        #"HFpEF"= gwas.pasc.pef2,
        "BMI"= gwas.bmi,
        "DM"= gwas.dm,
        "NICM"= gwas.nicm,
        "LVEF"= gwas.lvef,
        "LVEVID" = gwas.lvedvi,
        "SVI"= gwas.svi,
        "SVil"= gwas.lveSvi
  )

  l2= map(l, function(x){
    f.gwas(gwas = x, alpha= 0.001)
  })

  map(l, function(x){
    hist(x$pvalue)
  })

  hfpef.leading.edge= c("MYH7"  ,
    "THBD"  ,
    "XDH"   ,
    "LOX" ,
    "DAB2IP" ,
    "SMAD6"  ,
    "BRAF"  ,
    "MERTK",
    "NUAK1" ,
    "TGFB2" ,
    "CBX7"  ,
    "GATA5",
    "ATF6"  ,
    "GATA3" ,
    "CRLF1" )


  l2.df= enframe(l2, name = "source", value= "target")%>% unnest(target)%>% mutate(mor= 1)
  l2.df %>% group_by(source)%>% count()

  dec.res.ora= decoupleR::run_ora(mat = mat, network = l2.df,n_up = 100,n_background = 21000,minsize = 0)
  dec.res.ora %>% print(n=100)
  d.res= decoupleR::run_ulm(mat = mat, network = l2.df)

  d.res %>% print(n=100)

  d.res%>%
    filter(condition %in% c("hfpef.prio", "hfref.prio"))%>%
    mutate(padj= p.adjust(p_value),
           sig= ifelse(padj<0.01, "***",
                       ifelse(padj<0.05, "**",
                              ifelse(padj<0.1, "*", "")))
           )%>%
    ggplot(aes(x= source, y= score, fill =condition, alpha=sig))+
    geom_bar(stat="identity", width=.5, position = "dodge")
    geom_col(aes(fill = condition))

  hfpef.set= gg%>% slice_max(order_by = hfpef.prio, n= 100)%>% pull(gene)
  hfref.set= gg%>% slice_max(order_by = hfref.prio, n= 100)%>% pull(gene)

  gwas.lvef%>% filter(gene_symbol %in% hfpef.set)%>% arrange(pvalue)
  gwas.lvef%>% filter(gene_symbol %in% hfref.set)%>% arrange(pvalue)

  gwas.lvedvi%>% filter(gene_symbol %in% hfpef.set)%>% arrange(pvalue)%>% print(n=100)
  gwas.lvedvi%>% filter(gene_symbol %in% hfref.set)%>% arrange(pvalue)

gwas.lvedvi %>% group_by(Status)%>% count()


gwas.svi%>% filter(gene_symbol %in% hfpef.set)%>% arrange(pvalue)
gwas.svi%>% filter(gene_symbol %in% hfref.set)%>% arrange(pvalue)

gwas.lveSvi%>% filter(gene_symbol %in% hfpef.set)%>% arrange(pvalue)
# enrich predicted genes --------------------------------------------------

unique(gwas.lvedvi$Status)

  mat= gwas.lvedvi%>% filter(Status != "DAVIES_FAIL_FAREBROTHER_FAIL")%>%
    mutate(lvedvi= -log10(pvalue))%>%
    select(gene_symbol, lvedvi)

  mat2= gwas.lvef%>%filter(Status != "DAVIES_FAIL_FAREBROTHER_FAIL")%>%
  mutate(lvef= -log10(pvalue))%>%
    select(gene_symbol, lvef)

  mat= mat %>% inner_join(mat2)

  mat= as.matrix(mat%>%column_to_rownames("gene_symbol"))


  hfpef.set= gg%>% slice_max(order_by = hfpef.prio, n= 100)%>% pull(gene)
  hfref.set= gg%>% slice_max(order_by = hfref.prio, n= 100)%>% pull(gene)
  hfpef.set= hfpef.leading.edge

  net2= rbind(enframe(hfpef.set)%>% mutate(name= "hfpef"),
              enframe(hfref.set)%>% mutate(name= "hfref")) %>%
    mutate(mor= 1,
           likelihood= 1)%>%
    rename(source= name,
           target= value)




  hfpef.set= gg%>% slice_max(order_by = hfpef.prio, n= 100)%>%
    select(gene, hfpef.prio)%>%
    mutate(name= "hfpef")%>%
    rename(mor= hfpef.prio)

  hfref.set= gg%>% slice_max(order_by = hfref.prio, n= 100)%>%
    select(gene, hfref.prio)%>%
    mutate(name= "hfref")%>%
    rename(mor= hfref.prio )

  net2= rbind(hfpef.set, hfref.set)%>%# mutate(mor= 1)%>%
    rename(source= name,
           target= gene) %>%
    mutate(mor= mor+1)


  decoupleR::run_mlm(mat = mat, network = net2)
  #decoupleR::run_fgsea(mat = mat, network = net2)
  decoupleR::run_wmean(mat = mat, network = net2)
  decoupleR::run_wmean(mat = mat, network = net2)

# permutate pvals ---------------------------------------------------------


gwas.hfpef2

NES= map(c(1:100), function(x){
set.seed(x)
df2 <- gwas.hfpef2
df2$logp = sample(df2$logp,replace = F)
ref.res= perform_fgsea(df2, "hfref",gsets= gsets,  nperm = 1000) %>% mutate(gwas= "hfref")

df2 <- gwas.hfref2
df2$logp = sample(df2$logp,replace = F)
pef.res= perform_fgsea(df2, "hfpef",gsets= gsets,  nperm = 1000) %>% mutate(gwas= "hfref")
c(ref.res$NES, pef.res$NES)

})




# do GSE ------------------------------------------------------------------


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
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=14000, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)
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

hfpef.set= gg%>% slice_max(order_by = hfpef.prio, n= 200)%>% pull(gene)
hfref.set= gg%>% slice_max(order_by = hfref.prio, n= 200)%>% pull(gene)

sets= split(l2.df$target, l2.df$source)

sets2= split(gwas.hits$Nearest.Gene, gwas.hits$ï..Trait)
GSE_analysis(hfpef.set, sets)%>% as_tibble()
GSE_analysis(hfref.set, sets)%>% as_tibble()

GSE_analysis(hfpef.set, sets2[-1])%>% as_tibble()
GSE_analysis(hfref.set, sets2[-1])%>% as_tibble()
