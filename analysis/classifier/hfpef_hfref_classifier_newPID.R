## use classifier_functions.R
## to test different feature spaces for classifier performance and Dimension reduction


#library(cowplot,lib.loc = "C:/Program Files/R/R-4.0.0/library")
library(cowplot, lib.loc = .libPaths()[1])
library(tidymodels)
library(tidyverse)
library(qdapTools)
library(ggrepel)

# loadings ----------------------------------------------------------------
source("analysis/utils/utils_classifier_ML.R")
source("analysis/utils/utils.R")

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")
phecodes= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/topfreq_disease.rds")

hf = read.csv("T:fsa04/MED2-HF-Comorbidities/data/RWH_March2020/levinson_comorbidities_in_hf_patients_2020-03-25.csv",
              sep = ";",
              na.strings=c("","NA")) %>% as_tibble


map(pids.list, length)
pids= c(pids.list$hfpef,pids.list$hfref)

# gener1al preprocessing includes removing HF variables from the data frame and filtering:
data = preprocess(data, pids= pids)
data= data %>% filter(PheCode %in% phecodes)
Phe_dic = data %>% distinct(PheCode, Phenotype, category)
# second preprocessing removes disease with lower count than 10
#data= preprocess2(data, pids = pids, filter_abs = 10)

# 1) use cross sectional comorbidity profiles ---------------------------------------------------------------

#### full
mod_df= table_to_model_frame(data)
mod_df= add_response_variable(mod_df)

dim(mod_df)
rownames(mod_df)

set.seed(20)

# do a universal cv_split
cv_splits <- vfold_cv(mod_df, strata = hf, pool = 0.1)

cross.all = wrap_ml(mod_df, cv_splits )
#cross.all.dim= get_dim_reductions(df, pids.list)

saveRDS(cross.all, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/hfpef_hfref_fit_newID_noicd2022.rds")
cross.all= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/classifier_output/hfpef_hfref_fit_newID_noicd2022.rds")

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
  ggplot(., aes(x= model, y= value, color = model))+
  geom_point(size= 3)+
  scale_color_manual(values=c("red", "black"))+
  theme_minimal()+
  scale_y_continuous(limits = c(0.5, 1))+
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

# test for correlation
cor.test(comp_feat$importance, abs(comp_feat$estimate))
# plot absolute comparison
ggplot(data= comp_feat, aes(x= abs(estimate), y= importance))+
  geom_point()


hfref= comp_feat %>% filter( importance >0.01) %>% arrange(desc(estimate)) %>% pull(PheCode)
hfpef= comp_feat  %>% filter( importance >0.01)  %>% arrange((estimate)) %>% pull(PheCode)

comp_feat = comp_feat %>%
  mutate(label = ifelse(importance>2 & abs(estimate)>0.1, Phenotype, ""))
comp_feat = comp_feat %>%
  mutate(label = ifelse(PheCode %in% c(hfpef[1:20], hfref[1:20]), Phenotype, "")
         )

comp_feat= comp_feat %>%
  mutate(category = ifelse(category =="NULL","injuries & poisonings", category ),
         category = ifelse(category =="injuries and poisonings","injuries & poisonings", category )
         )


p1= ggplot(data= comp_feat, aes(x= estimate, y= importance, color = category))+
  geom_point(size= 2)+
  scale_color_manual(values= col_vector)+
  # geom_label_repel(aes(label= label, color= category),
  #                 box.padding = 1,
  #                 max.overlaps = 80,
  #                 show.legend = FALSE,force = 0.5)+
  # #scale_y_log10()+
  labs(color = "Disease category")+
  theme_minimal()+
  #theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(x= "Elastic net parameter",
       y= "RF importance")


p1

p2= ggplot(data= comp_feat, aes(x= estimate, y= importance, color = category))+
  geom_point(size= 2)+
  scale_color_manual(values= col_vector)+
 geom_label_repel(aes(label= label, color= category),
                 box.padding = 1,
                 size=3,
                 max.overlaps = 70,
                 show.legend = FALSE,force = 0.5)+
  # #scale_y_log10()+
  labs(color = "Disease category")+
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(x= "Elastic net parameter",
       y= "RF importance")

unify_axis(p2)

p.df= comp_feat%>%
  #filter(PheCode %in% c(hfpef[1:25], hfref[1:25]))%>%
  filter(PheCode %in%feat.vector[1:50])%>%
  mutate(hf= ifelse(estimate<0, "HFpEF", "HFrEF"))

p3.1= p.df %>%
  filter(hf== "HFpEF") %>%
  ggplot(., aes(y= reorder(Phenotype, abs(estimate)), x= abs(estimate)))+
  geom_point()+
  geom_col(width=0.1)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(x= "Elatsic net estimate")

p3= p.df %>%
  ggplot(., aes(y= reorder(Phenotype, abs(estimate)), x= abs(estimate),
                color= category))+
  geom_point(aes(color= category), size= 3)+
  scale_color_manual(values= col_vector)+
  geom_col(width=0.1)+
  facet_grid(rows= vars(hf),scales= "free", space= "free_y" )+
  theme_minimal()+
  theme(axis.text.x = element_text(angle= 90, hjust= 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")+
  labs(x= "absolute parameter",
       y  = "Phenotype")

unify_axis(p3)
unify_axis(p1)

hist(comp_feat$estimate, breaks= 50)
hist(comp_feat$importance, breaks= 50)


saveRDS(comp_feat, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")
comp_feat= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")

# calculate the top 25% of most importante features
mod = rf_full_cross$variables
cut_off= quantile(mod$importance)[4]
mod %>% filter(importance> cut_off)
median(rf_full_cross$variables$importance)

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/classifier_features.pdf",
    width= 4,
    height= 8
    )
unify_axis(p2)
unify_axis(p1)+theme(axis.text.x = element_text(angle= 45, hjust= 1))
dev.off()

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/classifier_feature_numb50.pdf",
    width= 7,
    height= 8
)
unify_axis(p3)

dev.off()

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/classifier_feature_hyperparas.pdf",
    width= 5,
    height= 5
)
unify_axis(cross.all$lr$p.hyper)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

unify_axis(cross.all$rf$p.hyper)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()



p3.1= ggMarginal(p1+
             theme(legend.position = "bottom")+
             labs(col= "Disease\ncategory"),
             type="density", size=3, groupFill= F,margins = "both")

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/classifier_feature_numb.pdf",
    width= 8,
    height= 6
)
p3.1

dev.off()


# assess dis.sig diversity ------------------------------------------------

feat.vector= comp_feat%>% filter(estimate!= 0) %>% arrange(desc(abs(estimate))) %>%
  pull(PheCode)

df= comp_feat%>% filter(PheCode  %in% feat.vector[1:100])
table(df$estimate>0)

#x= "infectious diseases"

map(unique(comp_feat$category), function(x){

  comp_feat %>% filter(category == x) %>% pull(estimate)%>% hist(., breaks= 25)

  #lm(as.formula(paste0(x, " ~ estimate")), data= comp_feat)
})

p1= comp_feat%>% filter(PheCode  %in% feat.vector[1:100])%>%
  arrange(estimate)%>%
  mutate(label2= ifelse(estimate>0, "HFrEF", "HFpEF"))%>%
  ggplot(., aes(fill = category,y= abs(estimate),
                         x = label2))+
  geom_bar(position = "fill", stat = "identity")+
  scale_fill_manual(values= col_vector)+
  theme_bw()+
  #theme(panel.border = element_rect(size= 1, fill= NA))+
  labs(x= "",y= "%", fill  = "")

unify_axis(p1)

x= comp_feat%>% filter(estimate != 0)%>% arrange(estimate)%>%
  mutate(label2= ifelse(estimate>0, "HFrEF", "HFpEF"))

#check for relative contribution to estimats
df1= comp_feat%>% filter(PheCode  %in% feat.vector[1:100])%>%
  mutate(label2= ifelse(estimate>0, "HFrEF", "HFpEF"))

df1 %>% group_by(
                 label2)%>%
  summarise(s= sum(abs(estimate)) )

df1 %>% group_by(category,
                 label2)%>%
  summarise(s= sum(abs(estimate)) )%>%
  mutate(s.= ifelse(label2== "HFpEF", s/25.4, s/11))



x2= table(x$label2, x$category)
p2= t(x2)%>%as.data.frame()%>%
  ggplot(., aes(x= Var2, y= Freq, fill = Var1))+
  geom_bar(position = "stack", stat = "identity")+
  scale_fill_manual(values= col_vector)+
  theme_minimal()+
  theme(panel.border = element_rect(size= 1, fill= NA))+
  labs(x= "",y= "count", fill  = "")
x2= prop.table(table(x$label2, x$category),1)
p3= t(x2)%>%as.data.frame()%>%
  ggplot(., aes(x= Var2, y= Freq, fill = Var1))+
  geom_bar(position = "stack", stat = "identity")+
  scale_fill_manual(values= col_vector)+
  theme_minimal()+
  theme(panel.border = element_rect(size= 1, fill= NA))+
  labs(x= "",y= "%", fill  = "")
pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/feature_distribution.pdf",
    width= 3.5,
    height= 5
)
unify_axis(p1)+theme(axis.text.x= element_text(angle= 45, hjust=1 ))
unify_axis(p2)
unify_axis(p3)



# Check classifier for dc 1 &6 --------------------------------------------


mod_df= table_to_model_frame(data)
mod_df= add_response_variable(mod_df)

library(igraph)

cl.= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
se= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/full_clinic_info.rds")
sum.t= readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/patient_metadata_2022.rds")

net= igraph::as_data_frame(cl. , "vertices")%>% as_tibble()
dc.6= net %>% filter(group_louv== 6)%>% pull(name)
dc.1= net %>% filter(group_louv== 1)%>% pull(name)

fit_LR_subset= function(dc, mod_df){
  mod_df2= mod_df[,c("hf", dc)]

  f.p =sum.t %>% filter(sex== "f")%>% pull(pid)
  m.p= sum.t %>%filter(sex== "m")%>% pull(pid)

  mod_df2= mod_df2 %>% rownames_to_column("pid")  %>%
    mutate(sex = ifelse(pid %in% f.p,"f","u"),
           sex = ifelse(pid %in% m.p,"m",sex))%>%
    filter(sex %in% c("m", "f"))%>%
    mutate(sex= factor(sex, levels= c("f", "m"))) %>%
    mutate(hf= factor(hf, levels= c("hfpef", "hfref")))

  mod_df2= column_to_rownames(.data = mod_df2, var = "pid")

  colnames(mod_df2)= paste0("x", colnames(mod_df2))
  form  = paste(x =   unlist(colnames(mod_df2)[-1]), collapse= "+")

  fit = glm(formula= as.formula(paste0("xhf ~ ", form)),
            data = mod_df2,
            family= "binomial")
}

x= fit_LR_subset(dc.6, mod_df)
x2= fit_LR_subset(dc.1, mod_df)

diseases= dc.6

clusts= unique(net$group_louv)

p.estimates= map(clusts, function(y){
  #print(y)
  diseases= net %>% filter(group_louv== y)%>% pull(name)

  x= fit_LR_subset(diseases, mod_df)

  df= map(diseases, function(x){
    #print(x)
    f= fit_LR_subset(x, mod_df)
    coef(summary(f))
  })
  names(df)= diseases

  tag = map(df, function(x){
    x[2,4]
  })%>%
    unlist() %>%
    enframe(., value = "p_val")

  tag2 = map(df, function(x){
    x[2,1]
  })%>% unlist()%>% enframe(., value = "estimate")

  df= tag %>% left_join(tag2)%>%
    left_join(Phe_dic %>%
                rename(name= PheCode))%>%
    mutate(cluster= y)


})


comp_feat= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")
labels= comp_feat %>% arrange(desc(abs(estimate)))%>% slice(1:100)%>% pull(Phenotype)
labels= c("Neoplasm of uncertain behavior of breast",
          "Hypertensive heart disease" ,

          #"Atherosclerosis"                                            ,
           "Viral infection"                                            ,
           "Chronic pericarditis"                                       ,
           #"Postmenopausal bleeding",
          "Osteopenia",
          "Pulmonary heart disease",
          #"Valvular heart disease/ heart chambers",
          "Renal failure NOS",
          "Rheumatoid arthritis",
          #"Diseases of white blood cells" ,
          "Hypercholesterolemia",
          "Sleep apnea"                 ,
           "Cardiogenic shock"  ,
           "Mitral valve disease"   ,
          #"Secondary malignant neoplasm of digestive systems",
          "Inflammatory and toxic neuropathy" ,
              "Arrhythmia (cardiac) NOS" ,
          "Chronic renal failure [CKD]" ,
          "Type 2 diabetes" ,
          "Essential hypertension"   ,
          "Coronary atherosclerosis" ,
          #"Pericarditis" ,
          "Myocardial infarction"  ,
          #"Acquired hypothyroidism",
          "Ventricular fibrillation and flutter",
          #"Hyperparathyroidism"      ,
          "Tobacco use disorder"








)

plot.df= do.call(rbind,p.estimates) %>%
  mutate(p_val= p.adjust(p_val, "BH"))%>%
  mutate(sig= ifelse(p_val<0.001, "<0.001",
                     ifelse(p_val<0.01, "<0.01",
                            ifelse(p_val<0.05, "<0.05", "ns")
                            )
                     ),
         sig= factor(sig, levels=c("<0.001", "<0.01", "<0.05", "ns")),
         label = ifelse(Phenotype %in% labels, Phenotype, "")
         )


library(ggrepel)

p.distro=
plot.df%>%mutate(dummy= 1)%>%
  ggplot(., aes(y=  estimate,
                x= as.factor(cluster),
                fill=  as.factor(cluster),
                #
                label= label))+

  geom_hline(yintercept = 0, lty= 2)+
  geom_jitter(aes(col = sig), size= 3)+
  geom_boxplot(alpha= 0.9, width= 0.3, outlier.colour  = NA)+
  #scale_color_manual(values = rev(c("darkgrey", "#AA1430", "black", "#D62747")))+
  scale_color_manual(values = c("black",cols.nice[1:3]))+
  scale_fill_manual(values = col.set)+
  #facet_grid(cols= vars(cluster))+
  theme_bw()
  #scale_fill_gradient(low= "black" , high =  "yellow")+
  geom_label_repel(aes(label= label),
                   alpha= 0.8 ,
                   size= 3,
                   max.overlaps = 100)
  coord_flip()

  # theme(axis.text.x  = element_blank(),
  #       axis.title.x = element_blank())

  map(clusts, function(x){
    pos <- position_jitter(width = 0.3, seed = 2)
    plot.df%>%mutate(dummy= 1)%>%
      filter(cluster== x)%>%
      ggplot(., aes(y=  estimate,
                    x= dummy,
                    fill=  as.factor(cluster),
                    #
                    label= label))+

      geom_hline(yintercept = 0, lty= 2)+
      geom_jitter(aes(col = sig), size= 3, position = pos)+
      geom_boxplot(alpha= 0.9, width= 0.3, outlier.colour  = NA)+
      #scale_color_manual(values = rev(c("darkgrey", "#AA1430", "black", "#D62747")))+
      scale_color_manual(values = c("black",cols.nice[1:3]))+
      scale_fill_manual(values = col.set)+
      #facet_grid(cols= vars(cluster))+
      theme_bw()+
      #scale_fill_gradient(low= "black" , high =  "yellow")+
      geom_label_repel(aes(label= label),
                       alpha= 0.8 ,
                       size= 3,
                       max.overlaps = 100, position = pos)


  })


  plot.df%>%mutate(dummy= 1)%>%
    ggplot(., aes(y=  estimate,
                  x= dummy,
                  fill=  as.factor(cluster),
                  #
                  label= label))+

    geom_hline(yintercept = 0, lty= 2)+
    geom_jitter(aes(col = sig), size= 3, position = pos)+
    geom_boxplot(alpha= 0.9, width= 0.3, outlier.colour  = NA)+
    #scale_color_manual(values = rev(c("darkgrey", "#AA1430", "black", "#D62747")))+
    scale_color_manual(values = c("black",cols.nice[1:3]))+
    scale_fill_manual(values = col.set)+
    facet_grid(cols= vars(cluster))+
    theme_bw()+
  #scale_fill_gradient(low= "black" , high =  "yellow")+
  geom_label_repel(aes(label= label),
                   alpha= 0.8 ,
                   size= 3,
                   max.overlaps = 100, position = pos)
  coord_flip()

  # theme(axis.text.x  = element_blank(),
  #       axis.title.x = element_blank())







  plot.df%>%
    ggplot(., aes(x= reorder(Phenotype,estimate),
                  y= estimate,
                  col= sig,
                  label= label))+
    geom_point()+
    scale_color_manual(values = rev(c("darkgrey", "#AA1430", "#E45570", "#D62747")))+
    scale_fill_gradient(low= "black" , high =  "yellow")+
    facet_grid(cols= vars(cluster))+
    theme_minimal()+
    theme(axis.text.x  = element_blank(),
          axis.title.x = element_blank())+
    geom_hline(yintercept = 0)+
    geom_label_repel(aes(label= label))
