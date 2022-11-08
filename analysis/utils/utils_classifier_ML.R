## wrapper functions to perform unsupervised  patient stratification
# 1. Visualization techniques
#     a)  PCA
#     b) tsne
#     c) UMAP
#2. build HF classifier (baseline models)
#     a) LR
#     b) RF
#3. tf idf transformation

## dependencies:

library(Rtsne)
#library(tidyverse)
#library(rlang)
library(tidymodels)
#library(readr)
library(factoextra)
library(umap)
library(ggExtra)


col.set = c("#fb8500",
            "#0099cc",
            "#f60000",
            "#fde800",
            "#00b600",
            "#00896e",
            "#0053c8",
            "#ff3399",
            "#b374e0",
            "#00cc99",
            "#fbc9e2",
            "#8a7437",
            "#2d3b61",
            "#612d3f",
            "#c9386b",
            "#357785",
            "#91b5bd",
            "#564e73")

# ## 1.visualizations --------------------------------------------------------------

# 1.a PCA
#' @param df, data frame with rownames as patient id (pids), and columns are features.
#' @param pids.list, list of pids here very specific for hfpef vs hfref, needs generalization
#' @return list object with pca df and scree plot

do.pca = function(df, pids.list){
  pca.res= princomp(x = df, cor = F )

  pca_meta= summary(pca.res)

  pca.plot = tibble(pid= rownames(df)) %>%
    mutate(hf= ifelse(pid %in% pids.list$hfpef,
                      "hfpef",
                      "hfref")
                      # ifelse(pid %in% pids.list$hfref,
                      #        "hfref",
                      #        "unlabeled"))
  )

  print(paste0("checkpoint: pid order is correct:",
               all(as.character(pca.plot$pid) == names(pca_meta$scores[,1]))))

  pca.plot = pca.plot %>% mutate(PCA1 = pca_meta$scores[,1],
                                 PCA2 = pca_meta$scores[,2]) %>%
    left_join(hf %>% select(pid, sex) %>% mutate(pid = as.character(pid)))

  PoV <- pca.res$sdev^2/sum(pca.res$sdev^2)

  p.pca = ggplot(pca.plot, aes(x= PCA1, y= PCA2, color = hf))+
    geom_point()+
    scale_color_manual(values= c("red", "black","grey"))+
    labs(x= paste0("PC1 ",round(PoV[1]*100, 2), "%"),
         y= paste0("PC2 ",round(PoV[2]*100, 2), "%"))
  p.pca= ggMarginal(p.pca, type="boxplot", size=10, groupFill= T)

  return(list(pca= p.pca,
              scree= fviz_eig(pca.res, ncp = 50),
              pca.obj= pca.res))

}

## 1.b tsne

#' @param df, data frame with rowsnames as patient id (pids), and columns are features.
#' @param pids.list, list of pids here very specific for hfpef vs hfref, needs generalization
#' @param perplexity, can be a single value or vector of perplexity values to run tsne with


do.tsne = function(df, #pids.list,
                   perplexity){

  #remove duplicates,
  df= df[!duplicated(df), ]

  pls= map(perplexity, function(x){

    tsne_res= Rtsne(df, check_duplicates = F, perplexity = x)

    # patient.df = tibble(pid= rownames(df)) %>% mutate(hf= ifelse(pid %in% pids.list$hfpef,
    #                                                            "hfpef",
    #                                                            ifelse(pid %in% pids.list$hfref,
    #                                                                   "hfref",
    #                                                                   "unlabeled")
    #                                                            )
    #                                                   )
    #
    # tsne.plot = as_tibble(cbind(as_tibble(tsne_res$Y), pid=  rownames(df))) %>%
    #   #mutate(hf= ifelse(pid %in% pids.list$hfpef,"hfpef", ifelse(pid %in% pids.list$hfref, "hfref", "hfmref")))
    #   mutate(hf= ifelse(pid %in% pids.list$hfpef,"hfpef", "hfref"))
    # p.tsne= ggplot(tsne.plot, aes(x= V1, y= V2, color = hf))+
    #   geom_point()+
    #   scale_color_manual(values= c("red", "black","grey"))+
    #   ggtitle("tsne")
    # p.tsne= ggMarginal(p.tsne,type="boxplot", size=10, groupFill= T)
    #
    #res= list(data=tsne_res,
    #         plot= p.tsne)
  })

  return(list(tsne_runs= pls, pids= rownames(df))
  )

}



## 1.c umap

#' @param df, data frame with rowsnames as patient id (pids), and columns are features.
#' @param pids.list, list of pids here very specific for hfpef vs hfref, needs generalization

do.umap= function(df, pids.list){
  library(ggExtra)
  umap_res= umap(df)

  umap.plot = as_tibble(cbind(as_tibble(umap_res$layout), pid=  rownames(df))) %>%
    mutate(hf= ifelse(pid %in% pids.list$hfpef,"hfpef", "hfref"),
           hf= factor(hf, levels= c("hfpef", "hfref") ))
    # mutate(hf= ifelse(pid %in% pids.list$hfpef,"hfpef",
    #                   ifelse(pid %in% pids.list$hfref, "hfref", "unlabeled")))
    #
  p.umap =ggplot(umap.plot, aes(x= V1, y= V2, color = hf))+
    geom_point()+
    scale_color_manual(values= c("red", "black","grey"))+
    ggtitle("umap")

  p.umap= ggMarginal(p.umap, type="boxplot", size=10, groupFill= T)

  return(list(data= umap_res,
         plot= p.umap))
  }


# ## 2.classifier - base line models -------------------------------------------

#. LR
# performs test-train split, does simple fit on the train model and reports results on test
#' @param df, data frame with rowsnames as patient id (pids), and columns are features.
#' @param pids.list, list of pids here very specific for hfpef vs hfref, needs generalization
#' @param ratio, elastic net ratio between L1 and L2 penalty
#for testing:
#seed= 20
#ratio=0.5

do.elasticnet1 = function(model_df,penalty, ratio, seed= 20 ){

  set.seed(seed)

  lr_mod <-
    logistic_reg(penalty = penalty,  mixture = ratio) %>%
    set_engine("glmnet")

  lr_recipe <-
    recipe(hf ~ ., data = model_df)
  #%>%
   # step_rm(c("428.2", "428.1"))

  lr_workflow <-
    workflow() %>%
    add_model(lr_mod) %>%
    add_recipe(lr_recipe)

  #validation_data <- mc_cv(hf_train, prop = 0.9, times = 10)

  #hf_trained <- prep(lr_recipe, training = hf_train, verbose = TRUE)

  #hf_trained %>% juice()

  lr_fit= fit(lr_workflow, model_df)


  lr_training_pred <-
    predict(lr_fit, model_df) %>%
    bind_cols(predict(lr_fit, model_df, type = "prob")) %>%
    # Add the true outcome data back in
    bind_cols(model_df %>%
                mutate(hf= as.factor(hf))%>%
                select(hf))

  roc_res= lr_training_pred %>%                # training set predictions
    roc_auc(truth = hf, .pred_hfpef)


  p1= autoplot(conf_mat(lr_training_pred, truth = hf, estimate = .pred_class), type ="heatmap")
  p2= autoplot(roc_curve(lr_training_pred,truth = hf, .pred_hfpef))

  coefs= tidy(lr_fit) %>% arrange(desc(estimate)) %>% rename(PheCode= term)

  return(list("conf.matrix"= p1,"roc_plot"= p2,
              "roc_num"= roc_res,
              "fit"= lr_fit,
              "variables"= coefs)
         )

}


## LR
# performs test-train split
# does 10fold CV on train split for each hyperparamter combination
# selects best combination,
# fits this model on whole train
# and reports on test
#
#' @param df, data frame with rowsnames as patient id (pids), and columns are features.
#' @param pids.list, list of pids here very specific for hfpef vs hfref, needs generalization

cv.elasticnet=function(df,pids.list, seed= 20){

  set.seed(seed)

  #### 1. Preprocessing data and split in train / test

  # create df with target column as classifier result:
  model_df= df %>% rownames_to_column("pid")  %>%
    mutate(hf = ifelse(pid %in% pids.list$hfpef,"hfpef","hf"),
           hf = ifelse(pid %in% pids.list$hfref,"hfref",hf),
           hf = ifelse(pid %in% pids.list$hfmref,"hfref",hf)) %>% ## important here we treat hfmref as hfref
    select(-pid) %>%
    filter(hf %in% c("hfpef", "hfref"))%>%
    mutate(hf= factor(hf, levels= c("hfpef", "hfref")))


  splits = initial_split(model_df, strata = hf)

  hf_train <- training(splits)
  hf_test  <- testing(splits)

  hf_train %>%
    count(hf) %>%
    mutate(prop = n/sum(n))

  hf_test %>%
    count(hf) %>%
    mutate(prop = n/sum(n))

  # val_set <- validation_split(hf_train,
  #                             strata = hf,
  # )
  #


  #### 2. Define general modelling approach

  # Model spec
  lr_mod <-
    logistic_reg(penalty = tune(),  mixture = tune()) %>%
    set_engine("glmnet")

  # Recipe (data transformation and formula)
  lr_recipe <-
    recipe(hf ~ ., data = hf_train)


  #### 3. Tuning hyperparameters (penalty and mixture)

  # create tuning grid:
  # v1 by hand:
  #lr_reg_grid <- tibble(penalty = 10^seq(-3, -1, length.out = 30))

  # v2 by using grid_regular:
  mixture_param <- parameters(penalty(), mixture())
  lr_reg_grid= grid_regular(mixture_param, levels= c(20, 10))

  # create cv splits to be tested for each hyperparameter combination
  train_folds <- vfold_cv(hf_train)

  # create workflow for tuning
  tune_wf <-
    workflow() %>%
    add_model(lr_mod) %>%
    add_recipe(lr_recipe)

  # fit models for cross validated tuning grid
  tune_res= tune_wf %>%
    tune_grid(
      resamples= train_folds,
      grid= lr_reg_grid
    )

  # define best models:
  best_mods= tune_res %>%
    show_best("roc_auc")

  # define best model:
  best_mod= tune_res %>%
    select_best("roc_auc")

  point=best_mods %>% filter(penalty== best_mod$penalty &
                               mixture== best_mod$mixture)

  # plot results of tuning:
  p.hyper = collect_metrics(tune_res) %>%
    mutate(mixture= factor(round(mixture,1)))%>%
    ggplot(., aes(x= penalty, y= mean))+
    geom_line(mapping = aes(color= mixture), size= 1.5, alpha = 0.4)+
    scale_color_viridis_d(option = "plasma", begin = .9, end = 0)+
    facet_grid(rows=vars(.metric))+
    geom_point(alpha = 0.2)+
    scale_x_log10(labels = scales::label_number())+
    geom_vline(xintercept= best_mods$penalty)+
    geom_point(data= point,aes( x= penalty, y= mean),
               size= 3,
               color = "red")+
    theme_minimal()


  #### 4. Use cross tuning parameters to test performance


  ## create a new workflow using the best hyperparameters
  final_wf<-
    tune_wf %>%
    finalize_workflow(best_mod)

  # get final TRAIN performance
  lr_res <-
    final_wf %>%
    fit(hf_train)

  lr_res %>%
    pull_workflow_fit() %>% tidy() %>%
    ggplot(aes(x= estimate))+
    geom_histogram()

  p.coeff= lr_res %>%
    pull_workflow_fit() %>% tidy() %>%
    arrange(desc(abs(estimate))) %>%
    #arrange(desc(estimate)) %>%
    slice(1:50) %>%
    left_join(data %>% rename(term= PheCode) %>% distinct(term, Phenotype), by= "term")%>%
    ggplot(aes(x =reorder(Phenotype, estimate), y = estimate)) +
    # geom_vline(xintercept = lambda_opt, lty = 3) +
    geom_col()+
    coord_flip()+
    theme_minimal()


  # get final TEST performance on initial split.
  final_res =
    final_wf %>%
    last_fit(splits)

  final_metrics= final_res %>%
    collect_metrics()

  p1= final_res %>%
    collect_predictions() %>%
    roc_curve(hf, .pred_hfpef) %>%
    autoplot()+
    geom_text(x=0.9, y=0, label=paste0("AUC: ",round(final_metrics[2,3], 3) ))

  p2= final_res %>%
    collect_predictions() %>%
    conf_mat(., truth= hf, estimate= .pred_class) %>%
    autoplot(type = "heatmap")

  tidy_coefs <- final_res

  final_res %>% pull_workflow_fit()



    broom::tidy() %>%
    filter(term != "(Intercept)") %>%
    select(-step, -dev.ratio)


  results= list(roc= p1,
                conf= p2,
                hyper= p.hyper,
                coeff= p.coeff,
                metrics= final_metrics)

  return()

  }

## LR
# use this function for model comparison.
# takes a cv split (for the use on whole data )
# does hyperparameter tuning
# reports mean performance (from 10cv) for best combination, as well as thos parameters

# Cross validate LR
#' @param df, data frame with rowsnames as patient id (pids), and columns are features.
#' @param pids.list, list of pids here very specific for hfpef vs hfref, needs generalization


do.elasticnet=function(model_df,cv_splits, seed= 20){


  set.seed(seed)
  #### 2. Define general modelling approach

  # Model spec
  lr_mod <-
    logistic_reg(penalty = tune(),  mixture = tune()) %>%
    set_engine("glmnet")

  # Recipe (data transformation and formula)
  lr_recipe <-
    recipe(hf ~ ., data = model_df)
  summary(lr_recipe)

  #### 3. Tuning hyperparameters (penalty and mixture)

  # create tuning grid:
  # v1 by hand:
  #lr_reg_grid <- tibble(penalty = 10^seq(-3, -1, length.out = 30))

  # v2 by using grid_regular:
  mixture_param <- parameters(penalty(), mixture())
  lr_reg_grid= grid_regular(mixture_param, levels= c(10, 10))


  # create workflow for tuning
  tune_wf <-
    workflow() %>%
    add_model(lr_mod) %>%
    add_recipe(lr_recipe)

  # fit models for cross validated tuning grid
  tune_res= tune_wf %>%
    tune_grid(
      resamples= cv_splits,
      grid= lr_reg_grid
    )

  # define best models:
  best_mods= tune_res %>%
    show_best("roc_auc")

  # define best model:
  best_mod= tune_res %>%
    select_best("roc_auc")

  point=best_mods %>% filter(penalty== best_mod$penalty &
                               mixture== best_mod$mixture)

  # plot results of tuning:
  p.hyper = collect_metrics(tune_res) %>%
    mutate(mixture= factor(round(mixture,1)))%>%
    ggplot(., aes(x= penalty, y= mean, color = mixture))+
    geom_line(size= 1.5, alpha = 0.4)+
    scale_color_viridis_d(option = "plasma", begin = 1, end = 0)+
    facet_grid(rows=vars(.metric))+
    geom_point(alpha = 0.8)+
    scale_x_log10(labels = scales::label_number())+
    geom_vline(xintercept= best_mods$penalty)+
    geom_point(data= point,aes( x= penalty, y= mean),
               size= 3,
               color = "red")+
    theme_minimal()


  #### 4. Use cross tuning parameters to test performance

    best_mods= best_mods %>% mutate(best= ifelse(penalty== best_mod$penalty &
                                      mixture== best_mod$mixture,
                                    "yes",
                                    "no"))

    results= list(res= tune_res,
                best_mods= best_mods,
                p.hyper= p.hyper)

  return(results)

}


do.LR=function(model_df,cv_splits, penalty, mix, seed= 20){


  set.seed(seed)
  #### 2. Define general modelling approach

  # Model spec
  lr_mod <-
    logistic_reg(penalty = penalty,  mixture = penalty) %>%
    set_engine("glmnet")

  # Recipe (data transformation and formula)
  lr_recipe <-
    recipe(hf ~ ., data = model_df)
  summary(lr_recipe)

  #### 3. Tuning hyperparameters (penalty and mixture)

  # create tuning grid:
  # v1 by hand:
  #lr_reg_grid <- tibble(penalty = 10^seq(-3, -1, length.out = 30))
  #
  # # v2 by using grid_regular:
  # mixture_param <- parameters(penalty(), mixture())
  # lr_reg_grid= grid_regular(mixture_param, levels= c(10, 10))
  #

  # create workflow for tuning
  tune_wf <-
    workflow() %>%
    add_model(lr_mod) %>%
    add_recipe(lr_recipe)

  # fit models for cross validated tuning grid
  tune_res= tune_wf %>%
    tune_grid(
      resamples= cv_splits
    )

  # define best models:
  best_mods= tune_res %>%
    show_best("roc_auc")


  return(best_mods)

}


############### RF

#. RF
#' @param df, data frame with rownames as patient id (pids), and columns are features.
#' @param cv_splits, split object from tidy models
#model_df= mod_df.pre[1:100, 7900:length(colnames(mod_df.pre))]
#model_df= mod_df
do.randomf = function(model_df, cv_splits, seed= 20 ){

  set.seed(seed)

  rf_mod <-
    rand_forest(mtry = tune(), trees = tune(), min_n= 1) %>%
    set_engine("ranger") %>%
    set_mode("classification")

  rf_recipe <-
    recipe(hf ~ ., data = model_df)
  #dim(model_df)
  rf_workflow <-
    workflow() %>%
    add_model(rf_mod) %>%
    add_recipe(rf_recipe)

  # use grid_regular:
  #mixture_param <- parameters(trees(),min_n() )

  #rf_reg_grid= grid_regular(mixture_param, levels= c(10, 10))

  # for mtry use 10-100% of n(parameters) and also add rule of thumb sqrt(n.param)
  n.param = dim(model_df)[2]-1 # subtract 1 for response variable

  #define tuning grid:
  rf_reg_grid= expand.grid(trees = c(10, 50, 100, 500),
                mtry= c(seq(n.param/10,n.param,n.param/10), sqrt(n.param), 5,15)) %>%
     mutate(mtry= as.integer(round(mtry))) %>%
     arrange(trees, mtry)

  rf_tune <-
    rf_workflow %>%
    tune_grid(resamples = cv_splits,
              grid = rf_reg_grid
              )

  # define best models:
  best_mods= rf_tune %>%
    show_best("roc_auc")

  # define best model:
  best_mod= rf_tune %>%
    select_best("roc_auc")
  #
  point=best_mods %>% filter(trees== best_mod$trees &
                               mtry== best_mod$mtry)
  #
  # # plot results of tuning:

  p.hyper = collect_metrics(rf_tune) %>%
    mutate(trees = as.factor(trees))%>%
    ggplot(., aes(x=mtry , y= mean, color = trees))+
    geom_line(size= 1.5, alpha = 0.4)+
    scale_color_viridis_d(option = "plasma", begin = .9, end = 0)+
    facet_grid(rows=vars(.metric))+
    geom_point(alpha = 0.8)+
    geom_vline(xintercept= best_mods$mtry)+
    geom_point(data= point, aes( x= mtry, y= mean),
               size= 3,
               color = "red")+
    theme_minimal()


  results= list(res= rf_tune,
                best_mods= best_mods,
                p.hyper= p.hyper)

 return(results)
}



#. LR
# performs test-train split, does simple fit on the train model and reports results on test
#' @param df, data frame with rowsnames as patient id (pids), and columns are features.
#' @param pids.list, list of pids here very specific for hfpef vs hfref, needs generalization
#' @param ratio, elastic net ratio between L1 and L2 penalty
#for testing:
#seed= 20
#ratio=0.5

do.randomf1 = function(model_df,pids.list,mtry,trees, seed= 20 ){

  set.seed(seed)

  rf_mod <-
    rand_forest(mtry = mtry, trees = trees, min_n= 1) %>%
    set_engine("ranger", importance = "impurity") %>%
    set_mode("classification")

  rf_recipe <-
    recipe(hf ~ ., data = model_df)

  rf_workflow <-
    workflow() %>%
    add_model(rf_mod) %>%
    add_recipe(rf_recipe)


  rf_fit= fit(rf_workflow, model_df)


  rf_training_pred <-
    predict(rf_fit, model_df) %>%
    bind_cols(predict(rf_fit, model_df, type = "prob")) %>%
    # Add the true outcome data back in
    bind_cols(model_df %>%
                mutate(hf= as.factor(hf))%>%
                select(hf))

  roc_res= rf_training_pred %>%                # training set predictions
    roc_auc(truth = hf, .pred_hfpef)


  p1= autoplot(conf_mat(rf_training_pred, truth = hf, estimate = .pred_class), type ="heatmap")
  p2= autoplot(roc_curve(rf_training_pred,truth = hf, .pred_hfpef))

  cofs= pull_workflow_fit(rf_fit)

  v.importance= cofs$fit$variable.importance

  v.importance= enframe(v.importance, name = "PheCode", value = "importance")

  return(list("conf.matrix"= p1,
              "roc_plot"= p2,
              "roc_num"= roc_res,
              "fit"= rf_fit,
              "variables"= v.importance)
  )

}



# ## 3.data transformation -----------------------------------------------------

## transform to tfidf space:
#' @data , data frame with patients (pid), and diseases (PheCode) per row
#' @return matrix (patients x diseases) with tf idf values

table_to_tfidf_frame= function(data){

  #1) keep per patient only diagnosis with different entry dates

  #df= data %>% distinct(pid, entry_date, PheCode) %>% drop_na

  #remove NAs


    ## 1 TF
  # count disease occurences per patient (tf)
  dfx= df %>%
    group_by(pid, PheCode) %>%
    count(PheCode) %>%
    ungroup() %>%
    pivot_wider( values_from= n, names_from= PheCode) %>%
    column_to_rownames("pid")

  # change NAs into 0
  dfx[is.na(dfx)] = 0

  #get total number of patients:
  n.pat= dim(dfx)[1]

  # normalize tf by max tf per patient
  tf.max= apply(dfx, 1, max)

  print(paste0("maximal term frequency range: ", range(tf.max)[1]," - ", range(tf.max)[2]))

  tf.norm= dfx / tf.max
  print(range(tf.norm))
  ## 2 idf

  # calculate number of patients that have disease (n) and the idf which
  # is log(n.pat / n)

  dis.occ= df %>%
    distinct(PheCode, pid) %>%
    count(PheCode) %>%
    mutate(idf= log(n.pat/n))

  # make sure tf.norm and idf data frames are in the same order
  tf.norm= tf.norm[, dis.occ$PheCode]

  message(paste("dims agree for tfidf calc = ", dim(tf.norm)[2]== length(dis.occ$idf)))

  # calculate tfidf by multiplying each COLUM of Tf.norm with the idf value for each disease:
  tfidf = sweep(tf.norm, MARGIN=2, dis.occ$idf, `*`)

  return(tfidf)

}


#transform to count matrix with disease as columns, patients as rows and 0 and 1 defining the disease profile

#' @param df, data frame with aptient info (must have pid and PheCode column)
#' @param char, boolean , if T transformed dataframe will be type character (for running MCA downstream)

table_to_model_frame = function(df, char= F){
  #this model frame counts disease occurence only once :
  df= df %>% distinct(PheCode, pid)  %>% drop_na()
  # data transformation into n x p Matrix
  df= split(as.character(df$PheCode), df$pid) %>%
    lapply(function(x) strsplit(x, "-")) %>%
    mtabulate()%>%
    qdapTools::matrix2df("pid") %>%
    column_to_rownames("pid")
  if(char== T){
    df2= data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
    rownames(df2) = rownames(df)

  }else{df2= df}

  return(df2)
}



# Helper functions (data transformation, filterting etc.) -----------------

# function to filter samples and features.
# samples= only samples in pids.list ( hfpef and hfref labeled patients)
# features= filter icd codes to not contain those that were used to define those patients


preprocess= function(data,
                     HF_allcause= c("I11.0", "I13.0", "I13.2", "I25.5", "I42.0", "I42.5", "I42.8", "I42.9", "I50.0", "I50.1", "I50.9"),
                     pids= pids.list[2:4]){

  data_red= data %>%
    filter(pid %in% pids) %>% # decide here which patient cohort to use
    filter(!is.na(PheCode)) %>%
    filter(!icd4 %in% HF_allcause) %>%
    filter(!PheCode %in% c("428.0",
                           "428.1",
                           "428.2"))



}

preprocess2= function(data,
                      filter_abs= 0,
                      filter_rel= 0,
                      pids){

  diseases= disease_frequencies(pids, data, feature= "PheCode") %>%
    filter(freq> filter_abs,
           rel_freq > filter_rel) %>% pull(PheCode)

  data2= data %>% filter(pid %in% pids,
                         PheCode %in% diseases
                         )

}


add_response_variable= function(df){

  model_df= df %>% rownames_to_column("pid")  %>%
    mutate(hf = ifelse(pid %in% pids.list$hfpef,"hfpef","hf"),
           hf = ifelse(pid %in% pids.list$hfref,"hfref",hf),
           hf = ifelse(pid %in% pids.list$hfmref,"hfref",hf)) %>% ## important here we treat hfmref as hfref
      filter(hf %in% c("hfpef", "hfref"))%>%
    mutate(hf= factor(hf, levels= c("hfpef", "hfref")))# %>%  drop_na()

  model_df= column_to_rownames(.data = model_df, var = "pid")

  return(model_df)
}

get_dim_reductions= function(df,
                             ...){

  #pca= do.pca(df, pids.list)
  #df = df[rownames(df)%in% unlist(pids.list[2:4]),]
  tsne= do.tsne(df, pids.list, perplexity = c(100))

  umap= do.umap(df, pids.list)

  results= list(#"pca"= pca,
    "tsne"= tsne,
    "umap"= umap)

  return(results)
}

wrap_ml = function(model_df, cv_splits){
  res_lr= do.elasticnet(model_df, cv_splits)
  message("lr done")
  res_rf= do.randomf(model_df, cv_splits)
  message("rf done")
  list("lr"= res_lr,
       "rf"= res_rf)
}


get_aucs = function(all_res){

  auc_res= sapply(names(all_res), function(x){
    map(models, function(y){
      all_res[[x]][[y]]$best_mods %>% arrange(desc(mean)) %>% slice(1) %>%
        pull(mean)

    })
  })
  as.data.frame(auc_res) %>% mutate(model= models) %>% pivot_longer(-model)%>%
    mutate(value= unlist(value))
}


# perform full fits  ------------------------------------------------------

#function useful to extract optimal hyperparameters fromlisted results from the hyperparameter training
#' @param results listed results basically, featurespace[[subset]][[model]]

get_best_hyper= function(results, model){

  hyperparam= results[[model]]$best_mods %>% arrange(desc(mean))

  if(model== "lr"){
    return(list("penalty"= hyperparam%>% top_n(1)%>% pull(penalty),
                "mixture"= hyperparam%>% top_n(1) %>% pull(mixture)))
  }else if(model == "rf"){
    hyperparam = hyperparam %>% arrange(desc(mean), trees) # filter for best model with less trees
    return(list("mtry"= hyperparam%>% top_n(1)%>% pull(mtry),
                "trees"= hyperparam%>% top_n(1) %>% pull(trees)))

  }
}

#function to fit either elastic net or random forest on whole training data (model functions are in srced)
# using optimal hyperparameters
#' @param model, either "rf" or "lr"
#' @param data, full data, preprocessing done in the function
#' @param tfidf, boolean if tfidf space should be used
#' @param hyp.param, out put from get_best_hyper, a list with model specific best hyperparametrs
#'
perform_full_fit= function(model= "lr", df,  tfidf= F, hyp.param) {


  if(tfidf == T){
    df2 = table_to_tfidf_frame(df)
  } else {
    df2= table_to_model_frame(df)
  }

  df2= add_response_variable(df2)
  df2= df2[,colSums(df2 %>% select(-hf))>0]

  #model_df= mod_df[1:100, 1000:1153]
  if(model == "lr"){
    fit_all= do.elasticnet1(model_df= df2,
                            #pids.list,
                            penalty = hyp.param$penalty,
                            ratio=hyp.param$mixture, seed= 20 )
  }else if(model =="rf"){
    fit_all= do.randomf1(model_df= df2,
                         pids.list,
                         mtry= hyp.param$mtry,
                         trees= hyp.param$trees
    )

  }


}


#function to query two modelfit objects to return a combined data frame for auroc comparison

get_roc_df= function(mod_fit1, mod_fit2){
  mod_pred1 <-
    predict(mod_fit1, mod_fit1$pre$mold$predictors) %>%
    bind_cols(predict(mod_fit1, mod_fit1$pre$mold$predictors, type = "prob")) %>%
    # Add the true outcome data back in
    bind_cols(mod_fit1$pre$mold$outcomes)

  mod_pred2 <-
    predict(mod_fit2, mod_fit2$pre$mold$predictors) %>%
    bind_cols(predict(mod_fit2, mod_fit2$pre$mold$predictors, type = "prob")) %>%
    # Add the true outcome data back in
    bind_cols(mod_fit2$pre$mold$outcomes)

  df.plot= rbind(roc_curve(mod_pred1,truth = hf, .pred_hfpef)%>% mutate(model="rf"),
                 roc_curve(mod_pred2,truth = hf, .pred_hfpef)%>% mutate(model="lr")) %>%
    mutate(model= factor(model))

}
