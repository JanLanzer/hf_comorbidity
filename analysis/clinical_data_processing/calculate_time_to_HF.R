## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2023-02-03
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2023
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
##  We will calculate for each patient the time distance to his first HF diagnosis
## ---------------------------

library(tidyverse)
library(lubridate)
library(qdapTools)

source("~/GitHub/RWH_analysis/scripts/utils.R")
# read tables.

directory= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/"

icd = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")

pid.df= read.csv( "T:/fsa04/MED2-HF-Comorbidities/data/RWH_September2022/raw/levinson_comorbidities_pids_2022-10-06.csv", sep = ";")
#add age calculation
pid.df= pid.df  %>%
  mutate(birthday= lubridate::as_date(birthday),
         entry_date = lubridate::as_date(diag_date_min ),
         ICDint = lubridate::interval(birthday, entry_date),
         age.at.icd = ICDint/dyears(1)) %>%
  select(-ICDint)%>%
  as_tibble()

#load features
phecodes= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/topfreq_disease.rds")

pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")
map(pids.list,length)

#load comorbidity profiles_
comp_feat= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")
feat.vector= comp_feat%>% filter(estimate!= 0) %>% arrange(desc(abs(estimate))) %>%
  pull(PheCode)

comp_feat= comp_feat%>% filter(PheCode  %in% feat.vector[1:100])

ML.class= comp_feat %>% mutate(hf = ifelse( estimate< 0, "hfpef",
                                            ifelse(estimate>0, "hfref", "ns")
)
)


#pheno.data= readRDS(file= "T:/fsa04/MED2-HF-Comorbidities/data/processed_data/full_clinic_info.rds")
#pheno.data2= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/patient_metadata_2022.rds")



# calculate for each diagnosis in our cohort and comorbidities tim --------

#get df with earliest HF diagnosis
dfX= pid.df %>%
  select(pid,diag_date_min)%>%
  filter(pid %in% pids.list$hf_all)

# filter pid and codes, join with early hf diagnosis and calculate interval of each
# diagnosis in months
time.calc= icd%>%
  filter(pid %in% pids.list$hf_all,
         PheCode %in% phecodes)%>%
  left_join(dfX, by= "pid")%>%  #join with earlies HF
  mutate(time_to_HF= lubridate::interval(diag_date_min , entry_date)/dmonths(1))%>%
  select(time_to_HF, everything())%>%
  mutate(hf.type = ifelse(pid %in% pids.list$hfref,
                          "hfref",
                          ifelse(pid %in% pids.list$hfpef,
                                 "hfpef",
                                 ifelse(pid %in% pids.list$hfmref,
                                        "hfmref",
                                        "unlabeled")
                                 )
                          )
         )


## plot results for all for full data

#by cohort
time.calc %>%
  filter(hf.type!= "unlabeled")%>%
  group_by(pid, PheCode)%>%
  mutate(min_time_to_HF = min(time_to_HF))%>%
  ggplot(.,aes(x= min_time_to_HF, fill= hf.type))+
  facet_grid(rows= vars(hf.type), scales= "free")+
  geom_histogram(bins = 100)+
  labs(x= "time to HF (months)")
  #scale_y_log10()+
  geom_vline(xintercept = 0)



time.calc%>%
    filter( hf.type!= "unlabeled",)%>%
    ggplot(aes(x= hf.type, y= time_to_HF, fill= hf.type))+
    geom_hline(yintercept = c(-12,12,-24, 24,0))+
    geom_violin()+
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA)+
    stat_compare_means(method = "wilcox.test",comparisons =  list(c("hfpef", "hfref"),
                                                                c("hfpef", "hfmref"),
                                                                c("hfmref", "hfref")
                                                                )
                     )

  #by disease cat
time.calc %>%
  filter(hf.type!= "unlabeled")%>%
  ggplot(.,aes(x= time_to_HF, fill= hf.type))+
  facet_grid(cols= vars(category), scales= "free_y")+
  geom_histogram(bins = 100)+
  #scale_y_log10()+
  geom_vline(xintercept = 0)


# filter for comorbidity profiles and test wheter time to hf is different
time.calc %>%
  filter(hf.type!= "unlabeled",
         PheCode %in% (ML.class%>% filter(hf== "hfref")%>% pull(PheCode)))%>%
  ggplot(.,aes(x= time_to_HF, fill= hf.type))+
  facet_grid(rows= vars(hf.type), scales= "free")+
  geom_histogram(bins = 100)+
  labs(x= "time to HF (months)")+


library(ggpubr)

#boxxplots
time.calc%>%
  left_join(ML.class)%>%
  filter(!is.na(hf) &  hf.type!= "unlabeled",)%>%
  ggplot(aes(x= hf.type, y= time_to_HF, fill= hf))+
  geom_boxplot()
  stat_compare_means(method = "wilcox.test",comparisons =  list(c("hfpef", "hfref"),
                                                                c("hfpef", "hfmref"),
                                                                c("hfmref", "hfref")
  )
  )



# plot patient visits over time -------------------------------------------

df= icd%>%  filter(pid %in% pids.list$hf_all,
               PheCode %in% phecodes) %>%
    mutate(hf.type = ifelse(pid %in% pids.list$hfref,
                            "hfref",
                            ifelse(pid %in% pids.list$hfpef,
                                   "hfpef",
                                   ifelse(pid %in% pids.list$hfmref,
                                          "hfmref",
                                          "unlabeled")
                            )
    )
    )%>%
    distinct(pid, entry_date, hf.type)


df  %>%
  filter(hf.type != "unlabeled")%>%
  ggplot(aes(x= as.Date(entry_date),fill= hf.type))+
  facet_grid(rows= vars(hf.type))+
  geom_histogram(aes(y=..count../sum(..count..)))+
  theme(axis.text.x = element_text(angle= 90, hjust= 1))

df  %>%
  filter(hf.type != "unlabeled")%>%
  ggplot(aes(x= as.Date(entry_date),fill= hf.type))+
  facet_grid(rows= vars(hf.type))+
  geom_histogram(bins= 100)+
  theme(axis.text.x = element_text(angle= 90, hjust= 1))




df  %>%
  filter(hf.type != "unlabeled")%>%
  ggplot(aes(y= as.Date(entry_date),x= hf.type))+
  geom_boxplot()+
  stat_compare_means(method = "wilcox.test",comparisons =  list(c("hfpef", "hfref"),
                                                                c("hfpef", "hfmref"),
                                                                c("hfmref", "hfref")
  )
  )

  facet_grid(rows= vars(hf.type))+
  theme(axis.text.x = element_text(angle= 90, hjust= 1))

# test whether the time to HF as a variable to in a log regression  ----------------
source("analysis/utils/utils_classifier_ML.R")
source("analysis/utils/utils.R")

cl.= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
se= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/full_clinic_info.rds")
sum.t= readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/patient_metadata_2022.rds")


## prepare the modelling data with all comorbidities and hfsubtype as response V:

prep_data= function(time.calc){
  mod_df= table_to_model_frame(time.calc)
  mod_df= add_response_variable(mod_df)

  ## additional data will be sex, time to diagnosis and age
  #ML.class$PheCode

  dis = colnames(mod_df)[colnames(mod_df)  %in% ML.class$PheCode]
  mod_df2= mod_df[,c("hf", dis)]


  ## add sex
  f.p =sum.t %>% filter(sex== "f")%>% pull(pid)
  m.p= sum.t %>%filter(sex== "m")%>% pull(pid)

  mod_df2= mod_df2 %>% rownames_to_column("pid")  %>%
    mutate(sex = ifelse(pid %in% f.p,"f","u"),
           sex = ifelse(pid %in% m.p,"m",sex))%>%
    filter(sex %in% c("m", "f"))%>%
    mutate(sex= factor(sex, levels= c("f", "m"))) %>%
    mutate(hf= factor(hf, levels= c("hfpef", "hfref")))

  ## add age
  #sum.t
  mod_df2= mod_df2 %>% inner_join(sum.t %>% select(pid, age.at.icd)%>% mutate(pid= as.character(pid)))

}

mod_df2= prep_data(time.calc)
## fit each model now for a single disease and include the time to HF
df= map(ML.class$PheCode, function(dis){

  #add time do hf


  t.vec= time.calc %>%
    filter(PheCode== dis) %>%
    select(pid, time_to_HF) %>%
    group_by(pid)%>%
    summarise(time= min(abs(time_to_HF)))%>%
    mutate(pid= as.character(pid))

  # mod_df3= mod_df2[, c("pid", "hf", dis,"sex", "age.at.icd" ) ]%>%
  #   left_join(t.vec)%>%
  #   mutate(time= ifelse(is.na(time), 0, time))
  mod_df3= mod_df2[, c("pid", "hf", dis) ]%>%
    left_join(t.vec)%>%
    mutate(time= ifelse(is.na(time), 0, time))


  mod_df3= column_to_rownames(.data = mod_df3, var = "pid")

  colnames(mod_df3)[2] = paste0("x", dis)
  form  = paste(x =   unlist(colnames(mod_df3)[-1]), collapse= "+")

  fit = glm(formula= as.formula(paste0("hf ~ ", form)),
          data = mod_df3,
          family= "binomial")
})


df2= map(df, function(x){
  coef(summary(x))%>% as.data.frame()%>%
    rownames_to_column("coefficient")%>%
    as_tibble() %>%
    mutate("model"=  names(x$coefficients[2]))
})

bind_rows(df2)%>%
  filter(grepl("x", coefficient) | coefficient== "time")%>%
  mutate(coefficient= ifelse(grepl("x", coefficient),"comorbidity", coefficient))%>%
  select(model,coefficient, Estimate)%>%
  pivot_wider(names_from = coefficient, values_from = "Estimate")%>%
  ggplot(aes(x= comorbidity, y=time
             ))+
  geom_point()

bind_rows(df2)%>%


bind_rows(df2)%>%
  filter(coefficient %in% c("sexm", "age.at.icd", "time"))%>%
  ggplot(aes(x= coefficient, y= model, fill = Estimate))+
  geom_tile()

bind_rows(df2)%>%
  filter(coefficient== "time")%>%
  arrange(`Pr(>|z|)`)

bind_rows(df2)%>%
  filter(model== "x401.21")%>%
  arrange(`Pr(>|z|)`)



##### next we will split the original data into diagnosis made before and after HF

preHF= time.calc%>% filter(time_to_HF< -6)
postHF= time.calc%>% filter(time_to_HF> 6)


preHF= prep_data(preHF)
postHF = prep_data(postHF)
allHF=  prep_data(time.calc)

# map different subset of data
res= lapply(list("pre"=preHF,
                 "post"= postHF,
                 "all"= allHF), function(mod_df2){

  df= map(ML.class$PheCode, function(dis){



  mod_df3= mod_df2[, c("pid", "hf", dis,"sex", "age.at.icd" ) ]

  mod_df3= column_to_rownames(.data = mod_df3, var = "pid")

  colnames(mod_df3)[2] = paste0("x", dis)
  form  = paste(x =   unlist(colnames(mod_df3)[-1]), collapse= "+")

  fit = glm(formula= as.formula(paste0("hf ~ ", form)),
            data = mod_df3,
            family= "binomial")
  })

})
#names(res) = c("pre", "post")
res2= lapply(names(res), function(df){
  df2= map(res[[df]], function(x){
    coef(summary(x))%>% as.data.frame()%>%
      rownames_to_column("coefficient")%>%
      as_tibble() %>%
      mutate("model"=  names(x$coefficients[2]))
  })
  bind_rows(df2)%>% mutate(data= df)

})
names(res2) = c("pre", "post", "all")

bind_rows(res2) %>%
  filter(grepl("x", coefficient))%>%
  mutate(sig= ifelse(`Pr(>|z|)`< 0.01, "**",
                     ifelse(`Pr(>|z|)`<0.05, "*", "")))%>%
  ggplot(aes(x= data, y= reorder(coefficient,Estimate),  fill = Estimate,label=sig))+
  geom_tile()+
  scale_fill_gradient2(low= "blue", mid= "white", high="red")+
  geom_text()

bind_rows(res2)%>%
  filter(grepl("x1", coefficient))

## check neoplasms:


neoplasms= ML.class%>% filter(category == "neoplasms")%>%
  mutate(PheCode= paste0("x", PheCode))%>%
  pull(PheCode)

x2= ML.class %>%
  filter(category == "neoplasms")%>%
  mutate(PheCode= paste0("x", PheCode))%>%
  rename(coefficient= PheCode,
         mult_LR= estimate)


bind_rows(res2)%>%
  left_join(x2)%>%
  filter(coefficient %in% neoplasms)%>%
  mutate(sig= ifelse(`Pr(>|z|)`< 0.01, "**",
                     ifelse(`Pr(>|z|)`<0.05, "*", "")))%>%
  ggplot(aes(x= data, y= reorder(Phenotype,Estimate),  fill = Estimate,label=sig))+
  geom_tile()+
  scale_fill_gradient2(low= "blue", mid= "white", high="red")+
  geom_text()+
  labs(x= "", y= "")



# test time blocks of patient visits to assess possible confoundin --------


#define time blocks:

time.sub= time.calc %>% filter(hf.type != "unlabeled")%>%filter(entry_date> as_date("2008-01-01"))
t.r= range(time.sub$entry_date)


number_of_breaks= 3
# we translate the start and endpoint into a time interval in months:
total_y= interval(start = as_date(t.r[1]), as_date(t.r[2]))/dyears()

one_int = total_y/number_of_breaks


time.blocks= lapply(seq(1,number_of_breaks,1),  function(x){
  print(x)
  #compute intervals
  end.t= one_int*x
  start.t= one_int*(x-1)

  #translate intervals into dates
  min.date= as_date(t.r[1]) %m+% years(round(start.t, 0))
  max.date= as_date(t.r[1]) %m+% years(round(end.t, 0))

  #subset date to time interval
  df= time.calc %>% filter(entry_date> min.date,
                           entry_date< max.date)

  #prepare model dataframe
  mod_df2= prep_data(df)


  #map every disease of interest to get a single comorbidity log regression
  #estimate
  df2= map(ML.class$PheCode, function(dis){
    print(dis)
    #if disease is not in df,  abort)
    if(! dis %in% colnames(mod_df2)){
      return(NULL)
    }

    #if the disease has less than three patients in that period than we cannot model it
    if(sum(mod_df2[,dis])<10){
      return(NULL)
    }

    mod_df3= mod_df2[, c("pid", "hf", dis,"sex", "age.at.icd" ) ]

    mod_df3= column_to_rownames(.data = mod_df3, var = "pid")

    colnames(mod_df3)[2] = paste0("x", dis)
    form  = paste(x =   unlist(colnames(mod_df3)[-1]), collapse= "+")

    fit = glm(formula= as.formula(paste0("hf ~ ", form)),
              data = mod_df3,
              family= "binomial")
  })


})

names(time.blocks)= paste0("time_intervall", 1:length(time.blocks))

res2= lapply(names(time.blocks), function(df){

  df2= map(time.blocks[[df]], function(x){

    if(is.null(x)){return(NULL)}

    y= coef(summary(x))%>% as.data.frame()%>%
      rownames_to_column("coefficient")%>%
      as_tibble() %>%
      mutate("model"=  names(x$coefficients[2]))
    return(y)
  })
  bind_rows(df2)%>% mutate(data= df)

})


names(res2) = names(time.blocks)


bind_rows(res2) %>%
  filter(grepl("x", coefficient))%>%
  mutate(sig= ifelse(`Pr(>|z|)`< 0.01, "**",
                     ifelse(`Pr(>|z|)`<0.05, "*", "")))%>%
  ggplot(aes(x= data, y= reorder(coefficient,Estimate),  fill = Estimate,label=sig))+
  geom_tile()+
  scale_fill_gradient2(low= "blue", mid= "white", high="red")+
  geom_text()


#use data frame form above to check possible outliers:

df= icd%>%  filter(pid %in% pids.list$hf_all,
                   PheCode %in% phecodes) %>%
  mutate(hf.type = ifelse(pid %in% pids.list$hfref,
                          "hfref",
                          ifelse(pid %in% pids.list$hfpef,
                                 "hfpef",
                                 ifelse(pid %in% pids.list$hfmref,
                                        "hfmref",
                                        "unlabeled")
                          )
  )
  )%>%
  distinct(pid, entry_date, hf.type, PheCode)


outliers= c("709.7", "800.1", "348.2")

map(outliers, function(x){
  df %>%
    filter(hf.type != "unlabeled",
           PheCode ==x)%>%
    ggplot(aes(x= as.Date(entry_date),fill= hf.type))+
    facet_grid(rows= vars(hf.type))+
    geom_histogram(bins= 50)+
    theme(axis.text.x = element_text(angle= 90, hjust= 1))+
    ggtitle(x)

})






