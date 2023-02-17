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
library(ggpubr)
library(ComplexHeatmap)

#source("~/GitHub/RWH_analysis/scripts/utils.R")
source("analysis/utils/utils_classifier_ML.R")
source("analysis/utils/utils.R")
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

# Add patient visit frequency ---------------------------------------------

visit.freq= time.calc %>% distinct(pid,hf.type, entry_date) %>% #group_by(pid)%>%
  count(pid)
time.calc= time.calc %>% left_join(visit.freq)

p.histo_visits= time.calc %>%
  filter(hf.type!= "unlabeled")%>%
  ggplot(aes(x= n))+
  geom_histogram(bins= 100)+
  facet_grid(rows = vars(hf.type))

p.histo_visits

p.box_visits= time.calc %>%
  filter(hf.type!= "unlabeled")%>%
  ggplot(aes(x= hf.type, y= n))+
  geom_violin()+
  geom_boxplot(width= 0.5)+
  scale_y_log10()+
  stat_compare_means(method = "wilcox.test",
                    comparisons = list(c("hfpef", "hfref"),
                                       c("hfmref", "hfpef"),
                                       c("hfref", "hfmref"))
                    )  +
  labs(x= "", y= "number of unique entry dates per patient")+
  theme_bw()+
  theme(axis.text = element_text(color="black", size= 11))

p.box_visits

# plot time_to_HF ---------------------------------------------------------


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



p.time_to_HF= time.calc%>%
    filter( hf.type!= "unlabeled",)%>%
    ggplot(aes(x= hf.type, y= time_to_HF))+
    geom_hline(yintercept = c(-12,12,-24, 24,0), color= "grey", lty= 2)+
    geom_violin()+
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA)+
    stat_compare_means(method = "wilcox.test",comparisons =  list(c("hfpef", "hfref"),
                                                                c("hfpef", "hfmref"),
                                                                c("hfmref", "hfref")
                                                                ))+
  labs(x= "", y= "time to first HF diagnosis (months)")+
  theme_bw()+
  theme(axis.text = element_text(color="black", size= 11))

p.time_to_HF


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
  geom_violin()
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
  filter(entry_date> as_date("2008-01-01"))%>%
  filter(hf.type != "unlabeled")%>%
  ggplot(aes(x= as.Date(entry_date),fill= hf.type))+
  facet_grid(rows= vars(hf.type))+
  geom_histogram(aes(y=..count../sum(..count..)))+
  theme(axis.text.x = element_text(angle= 90, hjust= 1))

p.visitsHIST= df  %>%
  filter(entry_date> as_date("2008-01-01")) %>%
  filter(hf.type != "unlabeled")%>%
  ggplot(aes(x= as.Date(entry_date),fill= hf.type))+
  facet_grid(rows= vars(hf.type))+
  geom_histogram(bins= 100)+
  theme(axis.text.x = element_text(angle= 90, hjust= 1))+
  theme_bw()+
  theme(axis.text = element_text(color="black", size= 11))+
  labs(x="date of comorbidity recording")

p.visitsHIST

p.visits.BOX= df  %>%
  filter(hf.type != "unlabeled")%>%
  filter(entry_date> as_date("2008-01-01"))%>%
  ggplot(aes(y= as.Date(entry_date),x= hf.type))+
  geom_violin()+
  geom_boxplot(width= 0.5)+
  stat_compare_means(method = "wilcox.test",comparisons =  list(c("hfpef", "hfref"),
                                                                c("hfpef", "hfmref"),
                                                                c("hfmref", "hfref")
  )
  )+
  labs(x= "", y= "date of comorbidity recording")+
  theme_bw()+
  theme(axis.text = element_text(color="black", size= 11))
p.visits.BOX

p.visits= cowplot::plot_grid(p.visitsHIST, p.visits.BOX, labels="AUTO",
                   rel_widths = c(1,0.5))


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/patient_visits_overtime.pdf",
    height= 6, width= 9)
p.visits
dev.off()


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/boxplots_time_confounders.pdf",
    height= 8, width= 4)

cowplot::plot_grid(#p.box_visits,
                   p.time_to_HF,
                   p.visits.BOX,
                   ncol = 1)

dev.off()

# test whether the time to HF as a variable to in a log regression  ----------------


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

  ## add n visits

  mod_df2 = mod_df2 %>%
    left_join(visit.freq %>% mutate(pid= as.character(pid)), by= "pid")


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



# test whether coefficents change if assesed before or after HF -----------


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

    if(! dis %in% colnames(mod_df2)){
      return(NULL)
    }

    #if the disease has less than three patients in that period than we cannot model it
    if(sum(mod_df2[,dis])<15){
      return(NULL)
    }




  mod_df3= mod_df2[, c("pid", "hf", dis,"sex", "age.at.icd") ]

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
    if(is.null(x)){return(NULL)}
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
  filter(coefficient =="n")%>%
  pull(Estimate)%>%
  hist()

bind_rows(res2)%>%
  filter(coefficient =="n")%>% print(n=2000)

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

time.sub= time.calc %>% filter(hf.type != "unlabeled")%>%
  filter(entry_date> as_date("2008-01-01"))

t.r= range(time.sub$entry_date)

# set number of equal time intervals:
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
    #print(dis)
    #if disease is not in df,  abort)
    if(! dis %in% colnames(mod_df2)){
      return(NULL)
    }

    #if the disease has less than three patients in that period than we cannot model it
    if(sum(mod_df2[,dis])<15){
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

names(time.blocks)=
  lapply(seq(c(1:number_of_breaks)), function(x){
    print(x)
    paste(as_date(t.r[1]) %m+% years(round(one_int*(x-1), 0)),
          as_date(t.r[1]) %m+% years(round(one_int*x, 0)),
          sep=" -> "
    )

  })%>% unlist()

res3= lapply(names(time.blocks), function(df){

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


names(res3) = names(time.blocks)

### plot coefficent of the diseases:

H.timeblocks= bind_rows(res3) %>%
  filter(grepl("^x", coefficient))%>%
  mutate(coefficient= str_replace_all(coefficient, "x", ""))%>%
  left_join(ML.class %>%
              rename(coefficient= PheCode)%>%
              select(Phenotype, coefficient, category))%>%
  mutate(sig= ifelse(`Pr(>|z|)`< 0.01, "**",
                     ifelse(`Pr(>|z|)`<0.05, "*", "")))%>%
  ggplot(aes(x= data, y= reorder(Phenotype,Estimate),  fill = Estimate,label=sig))+
  geom_tile(color="black")+
  scale_fill_gradient2(low= "blue", mid= "white", high="red", na.value = "darkgrey")+
  #scale_fill_continuous(na.value = 'yellow')+
  geom_text()

H.timeblocks %>%
  arrange(Estimate)%>%
  filter(Phenotype=="Malignant neoplasm of female breast")

H.timeblocks= unify_axis(H.timeblocks+theme_classic())+
  theme(axis.text.x = element_text(angle= 60, hjust= 1))+
  labs(x= "data subset to time interval",
       y= "")

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/timeblocks_coeff3.pdf",
    height= 15, width= 7)
H.timeblocks
  dev.off()






# plot large HMAP together: -----------------------------------------------

big.HMAP= rbind(bind_rows(res2), bind_rows(res3))
saveRDS(big.HMAP, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/single_disease_LR_models.rds")
hmap.df= big.HMAP %>%
  filter(grepl("^x", coefficient))%>%
    mutate(coefficient= str_replace_all(coefficient, "x", ""))%>%
    left_join(ML.class %>%
                rename(coefficient= PheCode)%>%
                select(Phenotype, coefficient, category)
              )%>%
  mutate(sig= ifelse(`Pr(>|z|)`< 0.01, "**",
                     ifelse(`Pr(>|z|)`<0.05, "*", "")),
         data = str_replace_all(data, "pre", "pre HF"),
         data = str_replace_all(data, "post", "post HF"),
         data = str_replace_all(data, "all", "full data"))

df1= hmap.df%>% select(coefficient, Estimate,`Pr(>|z|)`, model , data, Phenotype, category, sig )
df2= ML.class%>% rename(Estimate= estimate, coefficient= PheCode)%>%
  mutate(data= "elastic.net",
         `Pr(>|z|)`= NA,
         model = "PheCode",
         sig= NA)%>%
  select(coefficient, Estimate,`Pr(>|z|)`, model , data, Phenotype, category, sig )

hmap.df2= rbind(df1, df2)

### Quantify agreement of direction:

dr= hmap.df2 %>%
  group_by(data)%>%
  mutate(s.Est= sign(Estimate))%>%
  group_by(coefficient)%>%
  mutate(sum. = sum(s.Est))

dr %>%
  filter(data %in% c("pre HF", "post HF"))%>%
  mutate(s.Est= sign(Estimate))

dr%>%
  ggplot(aes(y= sum., x= reorder(Phenotype,sum.)))+
  geom_point()
  scale_fill_gradient2(low= "blue", mid= "white", high="red", na.value = "darkgrey")+
  #scale_fill_continuous(na.value = 'yellow')+
  geom_text()+
  theme(axis.text.x = element_text(angle= 40, hjust=1))


hmap.df2%>%
    ggplot(aes(x= data, y= reorder(Phenotype,Estimate),  fill = Estimate,label=sig))+
    geom_tile(color="black")+
    scale_fill_gradient2(low= "blue", mid= "white", high="red", na.value = "darkgrey")+
    #scale_fill_continuous(na.value = 'yellow')+
    geom_text()+
  theme(axis.text.x = element_text(angle= 40, hjust=1))

## Complex Heatmap with blocks:
hmap.df%>%
  filter(category=="neoplasms")%>%
  print(n=100)
hmap.df2= hmap.df2 %>%
  select(data, Estimate, Phenotype)%>%
  pivot_wider(names_from= data, values_from = Estimate)%>%
  as.data.frame()%>%
  column_to_rownames("Phenotype")

hmap.df2= hmap.df2[, c("full data", "pre HF", "post HF",
       "2008-01-02 -> 2013-01-02",
       "2013-01-02 -> 2018-01-02" ,
       "2018-01-02 -> 2022-01-02" )]

df.P= hmap.df%>%
  select(data, `Pr(>|z|)`, Phenotype)%>%
  pivot_wider(names_from= data, values_from = `Pr(>|z|)`)%>%
    as.data.frame()%>%
  column_to_rownames("Phenotype")


df.P= df.P[, c("full data", "pre HF", "post HF",
               "2008-01-02 -> 2013-01-02",
               "2013-01-02 -> 2018-01-02" ,
               "2018-01-02 -> 2022-01-02" )]


big.hmap = Heatmap(hmap.df2[, c("full data", "pre HF", "post HF",
               "2008-01-02 -> 2013-01-02",
               "2013-01-02 -> 2018-01-02" ,
               "2018-01-02 -> 2022-01-02" )],
        cluster_columns = F, row_names_side = "left",
        row_dend_side = "right",
        name= "Logistic\nRegression\nEstimate",
        column_names_rot = 60,
        row_names_max_width = unit(15, "cm"),
        column_split = c("full",
                         "time to HF","time to HF",
                         "subset to years","subset to years","subset to years"),
        # cell_fun = function(j, i, x, y, w, h, fill) {
        #   if(df.P[i, j] < 0.01) {
        #     grid.text("**", x, y)
        #   } else if(df.P[i, j] < 0.05) {
        #     grid.text("*", x, y)
        #   }
        # }

        )
big.hmap

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/combined_time_hmap.pdf",
    height= 14.6, width= 9)
draw(big.hmap)
dev.off()

big.hmap_slim= Heatmap(hmap.df2[, c("full data", "pre HF", "post HF",
                          "2008-01-02 -> 2013-01-02",
                          "2013-01-02 -> 2018-01-02" ,
                          "2018-01-02 -> 2022-01-02" )],


                   cluster_columns = F,
                   show_row_names = F,
                   row_dend_side = "right",
                   name= "Logistic\nRegression\nEstimate",
                   column_names_rot = 90,
                   #row_names_max_width = unit(15, "cm"),
                   column_split = c("full",
                                    "time to HF","time to HF",
                                    "years","years","years"),

                   # cell_fun = function(j, i, x, y, w, h, fill) {
                   #   if(df.P[i, j] < 0.01) {
                   #     grid.text("**", x, y)
                   #   } else if(df.P[i, j] < 0.05) {
                   #     grid.text("*", x, y)
                   #   }
                   # }

)
big.hmap_slim
pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/combined_time_hmap_slim2.pdf",
    height= 14.6, width= 4)
draw(big.hmap_slim)
dev.off()



# sandbox -----------------------------------------------------------------



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










