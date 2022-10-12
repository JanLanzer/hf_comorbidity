## use classifier_functions.R
## to test different feature spaces for classifier performance and Dimension reduction


#library(cowplot,lib.loc = "C:/Program Files/R/R-4.0.0/library")
library(cowplot, lib.loc = .libPaths()[1])
library(tidymodels)
library(tidyverse)
library(qdapTools)


# loadings ----------------------------------------------------------------
source("analysis/utils/utils_classifier_ML.R")

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/ICD10_labeled_phe.rds")

hf = read.csv("T:fsa04/MED2-HF-Comorbidities/data/RWH_March2020/levinson_comorbidities_in_hf_patients_2020-03-25.csv",
              sep = ";",
              na.strings=c("","NA")) %>% as_tibble

pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids.rds")
phecodes= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/top300_disease.rds")
map(pids.list, length)
pids= c(pids.list$hfpef,pids.list$hfref)

# gener1al preprocessing includes removing HF variables from the data frame and filtering:
data = preprocess(data, pids= pids)
data= data %>% filter(PheCode %in% phecodes)
# second preprocessing removes disease with lower count than 10
#data= preprocess2(data, pids = pids, filter_abs = 10)

# 1) use cross sectional comorbidity profiles ---------------------------------------------------------------

#### full
mod_df= table_to_model_frame(data)
mod_df= add_response_variable(mod_df)

rownames(mod_df)

set.seed(20)

# do a universal cv_split
cv_splits <- vfold_cv(mod_df, strata = hf, pool = 0.1)

cross.all = wrap_ml(mod_df, cv_splits )
#cross.all.dim= get_dim_reductions(df, pids.list)

saveRDS(cross.all, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/hfpef_hfref_fit_newID_noicd.rds")
cross.all= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/hfpef_hfref_fit_newID_noicd.rds")

# # 2) USE TF IDF SPACE -----------------------------------------------------
#
# #### full
#  tfidf.df= table_to_tfidf_frame(data)
#  tfidf.df= add_response_variable(tfidf.df)
#  tfidf.df= tfidf.df[,colSums(tfidf.df %>% select(-hf))>0]
# #
# # ####
# cv_splits <- vfold_cv(tfidf.df, strata = hf)
# res_lr= do.elasticnet(tfidf.df, cv_splits)
#
#  tfidf.all = wrap_ml(tfidf.df, cv_splits )
# tfidf.all.dim= get_dim_reductions(df, pids.list)
# #
# saveRDS(tfidf.all, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/hfpef_hfref_fit_tfidf_newID.rds")
# tfidf.all= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/hfpef_hfref_fit_tfidf_newID.rds")

# 4) plot models -------------------------------------------------------------

## 1) compare different feature spaces for best model performance:
all_res= list(cross= cross.all)#,

              #time= time.all)

models= c("rf", "lr")

auc_res= sapply(names(all_res), function(x){
  map(models, function(y){
    all_res[[x]][[y]]$best_mods %>% arrange(desc(mean)) %>% slice(1) %>%
      pull(mean)

  })
})

cross.all$rf


p.compare.feat = as.data.frame(auc_res) %>% mutate(model= models,
                                                   model= factor(model, levels= c("rf", "lr"))) %>%  pivot_longer(-model)%>%
  mutate(value= unlist(value)) %>%
  ggplot(., aes(x= name, y= value, color = model))+
  geom_point(size= 3)+
  scale_color_manual(values=c("red", "black"))+
  theme_minimal()+
  scale_y_continuous(limits = c(0.6, 0.9))+
  labs(x= "feature space",
       y= "mean AUROC (10CV)")


pdf(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/hfpef_hfref.pdf",
    width= 6,
    height = 6)
p.compare.feat
dev.off()


###### perform the full fit:

# fitting and variable interpretation -------------------------------------

## cross space
model= "rf"

# FIT
## RF
#tfidf space:
hyp.param = get_best_hyper(cross.all, "lr")
lr_full= perform_full_fit(model = "lr",
                          data,
                          tfidf = F,
                          #pids = c(pids.list$hfref, pids.list$hfpef),
                          hyp.param = hyp.param)
#cross space:
hyp.param = get_best_hyper(cross.all, "rf")
rf_full_cross= perform_full_fit(model = "rf",
                          data,
                          tfidf = F,
                          #pids = c(pids.list$hfref, pids.list$hfpef),
                          hyp.param = hyp.param)


saveRDS(rf_full_cross,"T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/full_fit_rf.rds")
rf_full_cross= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/full_fit_rf.rds")

## LR
#cross space:
hyp.param = get_best_hyper(cross.all, "lr")
lr_full= perform_full_fit(model = "lr", data, tfidf= F, hyp.param = hyp.param)

lr_full$variables = lr_full$variables%>%
  arrange(desc(abs(estimate))) %>% left_join(Phe_dic)

saveRDS(lr_full, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/full_fit_lr_cross.rds")
lr_full = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/full_fit_lr_cross.rds")

lr_full$fit$fit$fit$fit$df
lr_full$roc_plot

rf_full_cross$variables
rf_full_cross$conf.matrix

# plot variables

#1 compare the two aurocs
roc_data= get_roc_df(rf_full_cross$fit,lr_full$fit)

p.training.AUROC= roc_data  %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = model)) +
  geom_path(lwd = 1.5, alpha = 0.8, aes(col= model)) +
  geom_abline(lty = 3) +
  coord_equal()


#compare top important features

comp_feat= rf_full_cross$variables %>% left_join(lr_full$variables)
comp_feat= comp_feat %>%  left_join(data%>% distinct(PheCode, Phenotype, category))

#lr based disease:
comp_feat %>% filter(importance<1, abs(estimate)>0.1)
comp_feat %>% filter(importance>5, abs(estimate)>0) %>% arrange(desc(abs(estimate))) %>% print(n=100)
cor.test(comp_feat$importance, abs(comp_feat$estimate))

comp_feat %>% filter(importance<1)

ggplot(data= comp_feat, aes(x= abs(estimate), y= importance))+
  geom_point()
library(ggrepel)
comp_feat = comp_feat %>% mutate(label = ifelse(importance>2 & abs(estimate)>0.1, Phenotype, ""))
ggplot(data= comp_feat, aes(x= estimate, y= importance, color = category))+
  geom_point()+
  geom_text_repel(aes(label= label, color= category),
                  box.padding = 1,
                  max.overlaps = 80,
                  show.legend = FALSE,force = 0.5)+
  #scale_y_log10()+
  labs(color = "cat")+
  theme_minimal()

comp_feat%>%
  mutate(hf= ifelse(estimate<0, "hfpef", "hfref"))%>%
  filter(abs(estimate)>0.1)%>%
  ggplot(., aes(x= reorder(Phenotype, abs(estimate)), y= abs(estimate)))+
  geom_point()+
  geom_col(width=0.1)+
  facet_grid(cols= vars(hf), scales = "free_x" )+
  theme(axis.text.x = element_text(angle= 90, hjust= 1))

hist(comp_feat$estimate, breaks= 50)

hist(comp_feat$importance, breaks= 50)

boxplot(log10(comp_feat$importance))

ggplot(rf_full_cross$variables%>% mutate(x= "rf"), aes(y= importance,x= x))+
  geom_jitter()+
  geom_boxplot()


saveRDS(comp_feat, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")
comp_feat= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")
# calculate the top 25% of most importante features
mod = rf_full_cross$variables
cut_off= quantile(mod$importance)[4]
mod %>% filter(importance> cut_off)
median(rf_full_cross$variables$importance)


