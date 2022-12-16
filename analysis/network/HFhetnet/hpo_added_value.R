## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-06-10
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## test for differences between  HPO net and HFnet
## ---------------------------


# -------------------------------------------------------------------------

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")

auc_res= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/aurocs_multilayer_disease_pred.rds")
random.df= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/randomization_layer_dd_df.rds")
Phe_dic= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/icd10_phewas_dictionary.rds")

Phe_dic= Phe_dic %>% distinct(PheCode, Phenotype) %>% drop_na


## these results are from prediction of full network (with HPO layer):

withHPO = map(auc_res, function(x){
  x[[1]]$pr$auc.integral
})%>% unlist()

## these results are from prediction without the HPO layer:
withoutHPO= random.df$pr%>% filter(name== "real") %>% pull(value)

boxplot(withHPO, withoutHPO, names = c("HFnet+HPO", "HFnet"))
wilcox.test(withHPO, withoutHPO,
            paired = F
)

## auroc

## these results are from prediction of full network (with HPO layer):

withHPO = map(auc_res, function(x){
  x[[1]]$AUROC
})%>% unlist()

## these results are from prediction without the HPO layer:
withoutHPO= random.df$auroc%>% filter(name== "real") %>% pull(value)

boxplot(withHPO, withoutHPO, names = c("HFnet+HPO", "HFnet"))
wilcox.test(withHPO, withoutHPO,
            paired = F
)




# test whether disgenet confidence associates with prediciton perf --------

gd= process_gd(weight_cutoff = 0.29)


gd= gd%>% group_by(nodeB)%>% summarise(m.weight= median(weight),
                                   mean.weight= mean(weight))

df1 = map(auc_res, function(x){
  x[[1]]$pr$auc.integral
})%>% unlist()%>% enframe(name = "PheCode")%>% left_join(Phe_dic)
df2 = map(auc_res, function(x){
  x[[1]]$AUROC
})%>% unlist()%>% enframe(name = "PheCode")%>% left_join(Phe_dic)


complete= left_join(df2, gd %>% rename(PheCode= nodeB))

complete %>%ggplot(., aes(x=mean.weight, y= value))+
  geom_point()+
  geom_smooth()
