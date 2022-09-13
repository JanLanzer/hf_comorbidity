## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-09-13
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## this script will take the distinctive disease features and run different
## parameter cut offs to assess AUROC drops via LR
## ---------------------------

library(tidyverse)

comp_feat= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")



# redo  -------------------------------------------------------------------

source("analysis/utils/utils_classifier_ML.R")

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/ICD10_labeled_phe.rds")

hf = read.csv("T:fsa04/MED2-HF-Comorbidities/data/RWH_March2020/levinson_comorbidities_in_hf_patients_2020-03-25.csv",
              sep = ";",
              na.strings=c("","NA")) %>% as_tibble

pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids.rds")

phecodes= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/top300_disease.rds")

pids= c(pids.list$hfpef,pids.list$hfref)

# gener1al preprocessing includes removing HF variables from the data frame and filtering:
data = preprocess(data, pids= pids)
data= data %>% filter(PheCode %in% phecodes)

#### full
mod_df= table_to_model_frame(data)
mod_df= add_response_variable(mod_df)

rownames(mod_df)

set.seed(20)

# do a universal cv_split
cv_splits <- vfold_cv(mod_df, strata = hf)



# select testing grid -----------------------------------------------------
hist(comp_feat$estimate, breaks= 40)

test.cutoffs= seq(0.01,1.4, 0.03)

auros=  map(test.cutoffs, function(x){
  print(x)

  feats= comp_feat%>% filter(abs(estimate)>x)%>% pull(PheCode)

  sub.data= mod_df[,c("hf", feats)]

  res= do.elasticnet1(sub.data, penalty = 0, ratio = 0)

  res$roc_num$.estimate

})

data.frame(test.cutoffs , unlist(auros))%>% plot()

nfeat=  map(test.cutoffs, function(x){
  print(x)

  feats= comp_feat%>% filter(abs(estimate)>x)%>% pull(PheCode)

  length(feats)

})%>% unlist()

data.frame(test.cutoffs , nfeat)%>% plot()


feats= comp_feat%>% filter(abs(estimate)>0.5)%>% pull(PheCode)
