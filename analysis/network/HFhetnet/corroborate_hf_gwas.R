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
## gene hits

gg= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks.rds")


##### first mapping (not pascal, but minimal p)

gwas.hfpef= as_tibble(gwas3) %>% arrange(Pvalue)
gwas.hfref= as_tibble(gwas4) %>% arrange(Pvalue)


get_min_p = function(x){
  x%>%
    group_by(V8)%>%
    summarise(m.p= min(Pvalue))%>%
    arrange(m.p) %>%
    mutate(logp= -log10(m.p))
}

gwas.hfpef2 = get_min_p(gwas.hfpef)%>% rename(gene_symbol = V8)
gwas.hfref2 = get_min_p(gwas.hfref)%>% rename(gene_symbol = V8)


## network prediciton gene sets
gc.hfpef= gg %>% arrange(desc(hfpef.prio)) %>% slice(1:100) %>% pull(name)
gc.hfref= gg %>% arrange(desc(hfref.prio)) %>% slice(1:100) %>% pull(name)

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
ref.res= perform_fgsea(gwas.ref, "hfref",gsets = gsets, nperm = 1000) %>% mutate(gwas= "hfref")


# test gene cut offs ------------------------------------------------------
test_cutoffs= function(sqc= seq(10,500, 50),
                       gwas.hfpef,
                       gwas.hfref,
                       gg){

  l= map(sqc, function(x){

    gc.hfpef= gg %>% arrange(desc(hfpef.prio)) %>% slice(1:x) %>% pull(name)
    gc.hfref= gg %>% arrange(desc(hfref.prio)) %>% slice(1:x) %>% pull(name)

    gsets= list("hfpef"= gc.hfpef,
                "hfref"= gc.hfref)

    pef.res= perform_fgsea(gwas.hfpef, "hfpef",gsets = gsets, nperm= 2000) %>% mutate(gwas= "hfpef", cutoff= x)
    ref.res= perform_fgsea(gwas.hfref, "hfref", gsets=gsets, nperm = 2000) %>% mutate(gwas= "hfref", cutoff= x)

    df= rbind(pef.res, ref.res)
  })

  do.call(rbind, l) %>% as_tibble() %>% mutate(log10p= -log10(pval))%>%
    ggplot(., aes(x= cutoff, y= log10p, col= pathway))+
    facet_grid(rows= vars(gwas))+
    geom_point()+
    geom_path()+
    geom_hline(yintercept = -log10(0.05))


}


test_cutoffs(gwas.hfpef =gwas.pasc.pef2%>% mutate(logp= -log10(pvalue)),
             gwas.hfref = gwas.pasc.ref2%>% mutate(logp= -log10(pvalue)),
             gg= gg)


test_cutoffs(gwas.hfpef =gwas.pasc.pef%>% mutate(logp= -log10(pvalue)),
             gwas.hfref = gwas.pasc.ref%>% mutate(logp= -log10(pvalue)),
             gg= gg)

g.ranks %>%
  ggplot(., aes(x= reorder(name, -hfpef.prio), y= hfpef.prio))+
  geom_point()+
  geom_hline()



g.ranks %>% filter(hfpef.prio > 0.001)
  ggplot(., aes(x= reorder(name, -hfref.prio), y= hfref.prio))+
  geom_point()


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


