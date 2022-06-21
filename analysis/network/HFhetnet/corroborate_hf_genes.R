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
##  interpret and compare g.ranks
## ---------------------------

library(tidyverse)
g.ranks=readRDS( file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks.rds")
source("analysis/utils/utils_hetnet.R")

# test_drugs --------------------------------------------------------------

drg = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/drug_targets_processed.rds")

solid_drgs= drg%>% group_by(nodeB) %>% count %>% filter(n>2) %>% pull(nodeB)
xy= split.data.frame(drg%>% select(nodeA, nodeB), f = drg$nodeB)
xy= lapply(xy, function(x){x %>% pull(nodeA)})

test_drugs= map(xy[solid_drgs], function(x){
  # print(x)
  # print(length(x))
  # if(length(x)<5){
  #   return(NULL)
  # }
  res_= validate_results2(x, g.hfpef)
  res_$roc$auc
})

names(test_drugs)= solid_drgs

hist(unlist(test_drugs))
pef.drugs= enframe(unlist(test_drugs)) %>% arrange(desc(value))

test_drugs.ref= map(xy[solid_drgs], function(x){
  # print(x)
  # print(length(x))
  # if(length(x)<5){
  #   return(NULL)
  # }
  res_= validate_results2(x, g.hfref)
  res_$roc$auc
})


names(test_drugs.ref)= solid_drgs

hist(unlist(test_drugs.ref))
ref.drugs= enframe(unlist(test_drugs.ref)) %>% arrange(desc(value))

drugs= left_join(pef.drugs, ref.drugs, by= "name")%>%
  mutate(auc.diff= value.x-value.y)
drugs%>%
  ggplot(., aes(x= value.x, y= value.y))+
  geom_point()

drugs%>% arrange(desc(auc.diff))
validate_results2(set= hfpef, g.hfref)$roc
validate_results2(set= hfpef, g.hfpef)$roc


# test sc sigs ------------------------------------------------------------

### quick test for sc up regulated

fibs= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/genesets/cells_DEA_fibs.csv", sep = ";")%>%
  as_tibble
end= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/genesets/cells_DEA_end.csv", sep = ";") %>%
  as_tibble

gene_translate= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/gene_translate.rds")
gene_translate = gene_translate %>%
  rename(gene= MGI.symbol) %>%
  select(gene, Gene.name)

fibs= fibs %>% left_join(gene_translate)
end= end %>% left_join(gene_translate)


### simple enrich:
g.list= readRDS( file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks2.rds")
g.hfpef= g.list[[1]]
g.hfref= g.list[[2]]
top_fib_genes= fibs$Gene.name[1:100]
top_end_genes= end$Gene.name[1:100]

fibs.pef= validate_results_partial(top_fib_genes, g.hfpef)
fibs.ref= validate_results_partial(top_fib_genes, g.hfref)

end.pef= validate_results_partial(top_end_genes, g.hfpef)
end.ref= validate_results_partial(top_end_genes, g.hfref)


g.ranks %>% filter(grepl("ANGP", name))

## plot differences:
library(ggrepel)

fibs.pef$pr$auc.integral
fibs.ref$pr$auc.integral

fibs.ref$median_rank_stat
fibs.pef$median_rank_stat

## check fib signatures:
fib.hfpef= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/supplement_table3_hfpef.csv", sep = ";")%>%
  mutate(study= "hfpef")
fib.angII= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/supplement_table3_angii.csv", sep = ";")%>%
  mutate(study= "angii")
fib.mi= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/supplement_table3_mi.csv", sep = ";")%>%
  mutate(study= "mi")

df= rbind(fib.angII, fib.hfpef, fib.mi)%>% as_tibble() %>% rename(gene= ï..gene )%>% left_join(gene_translate)
df= df %>% drop_na(Gene.name)
gsets= split(df$Gene.name, df$study)


validate_fgsea= function(g.ranks, col= "hfref.prio", gsets){

  stats= g.ranks %>% arrange(desc({{col}})) %>% pull({{col}})
  names(stats)= g.ranks %>% arrange(desc({{col}})) %>% pull(name)
  fgseaSimple(pathways = gsets, stats = stats, nperm = 1000, scoreType = "pos")
}

validate_fgsea(g.ranks, "RW.value.hfpef", gsets = sets)

validate_fgsea(g.ranks, "RW.value.hfpef", gsets = gsets)
validate_fgsea(g.ranks, "RW.value.hfref", gsets = gsets)



# ReHeaT ------------------------------------------------------------------
sets= load_validation_genes(disgenet_value= 0.29, top_reheat = 500)
sets$set_phe= NULL
sets$set_reheat= NULL
sets$set_reheat_up= NULL
val.set= unique(unlist(sets))
length(unique(unlist(sets)))

validate_fgsea(g.ranks, "hfpef.prio", gsets = list("r"= val.set))
validate_fgsea(g.ranks, "hfref.prio", gsets = list("r"= val.set))
gsets= list("HF genes"= val.set,
            "ReHeaT"= sets$set_reheat)

gsea.pef= validate_fgsea(g.ranks, "RW.value.hfpef", gsets = gsets)%>% mutate(ranking= "HFpEF")
gsea.ref= validate_fgsea(g.ranks, "RW.value.hfref", gsets = gsets)%>% mutate(ranking= "HFrEF")

rbind(gsea.pef, gsea.ref)%>%
  ggplot(., aes(x= ranking, y= pathway, fill = -log10(pval)))+
  geom_tile()+
  scale_fill_gradient(low= "white", high = "red")

val.pef= validate_results_partial(val.set, g.hfpef)
val.ref= validate_results_partial(val.set, g.hfref)

res1= val.pef
res2= val.ref
compare_val_re= function(res1,
                         res2,
                         name1= "HFpEF",
                         name2= "HFrEF",
                         gset_name){
  #auroc
  df= rbind(c("AUROC", res1$AUROC, res2$AUROC),
        c("median_rank", res1$median_rank_stat, res2$median_rank_stat),
        c("AUC-PR", res1$pr$auc.integral, res2$pr$auc.integral)
  )
  colnames(df)= c("metric", name1, name2)

  as_tibble(df) %>% pivot_longer(-metric)%>% mutate(value= as.numeric(value))%>%
    mutate(gset= gset_name)


}
df=
rbind(compare_val_re(val.pef, val.ref, gset_name = "hf_set"),
      compare_val_re(fibs.pef, fibs.ref, gset_name = "fib_set"),
      compare_val_re(end.pef,end.ref, gset_name = "end_set")
)

ggplot(df , aes(x= metric, y= value, fill= name))+
  facet_grid(rows= vars(gset))+
  geom_col(position = "dodge", width= 0.5)


