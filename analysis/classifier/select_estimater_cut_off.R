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

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")

phecodes= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/topfreq_disease.rds")

cross.all= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/hfpef_hfref_fit_newID_noicd2022.rds")

hyp.param = get_best_hyper(cross.all, "lr")
rm(cross.all)



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

#lr_full= perform_full_fit(model = "lr", data, tfidf= F, hyp.param = hyp.param)

# select testing grid -----------------------------------------------------

feat.vector= comp_feat%>% filter(estimate!= 0) %>% arrange(desc(abs(estimate))) %>%
   pull(PheCode)

ddf= lapply(feat.vector[2:length(feat.vector)], function(x){
  print(x)

  feats= feat.vector[1:match(x,feat.vector )]

  sub.data= mod_df[,c("hf", feats)]

  res= do.LR(model_df = sub.data,
        cv_splits = cv_splits,
        penalty = hyp.param$penalty,
        mix = 0,  #we want all parameters included
        seed = 2    )

  res%>%
    mutate(nParam= length(feats))


})
ddf= ddf %>% do.call(rbind, .)

saveRDS(ddf, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/feat.selection.rds")

ddf %>% filter(nParam== 100)
p= ddf %>%
  ggplot(aes(y= mean, x= nParam))+
  geom_point(size = 0.9)+
  theme_minimal()+
  theme(panel.border = element_rect(size= 1, fill = NA))+
  labs(y= "CV mean AUROC",
       x= "n Parameter")



comp_feat %>% filter(PheCode %in% feat.vector[1:75])

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/classifier_cutoff.pdf",
    width= 5,
    height= 5
)
unify_axis(p)
dev.off()
# old version with training error (no CV) ---------------------------------


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
