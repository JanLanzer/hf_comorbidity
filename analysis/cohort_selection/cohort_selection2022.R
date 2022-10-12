## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-12-07
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##   This script implements a diagnostic algorithm to select HF cohorts from the RWH
# in the first part an inclusive data frame of nyha , bnp and echo data is generated
# then this data and icd10 data is used to generate different patient cohorts
# i.e. all_hf , hfpef, hfref, hfmref, finally those cohorts are used to plot different characteristics
##
## ---------------------------###


# libs and data -----------------------------------------------------------
library(ggpubr)
library(cowplot)
library(ggVennDiagram)
library(tidyverse)
library(lubridate)



source("analysis/utils/utils.R")

directory= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/"

# full data frame:
data= readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe.rds")


#new cohort 2022
pid.df= read.csv( "T:/fsa04/MED2-HF-Comorbidities/data/RWH_September2022/raw/levinson_comorbidities_pids_2022-10-06.csv", sep = ";")
#add age calculation
pid.df= pid.df  %>%
  mutate(birthday= lubridate::as_date(birthday),
         entry_date = lubridate::as_date(diag_date_min ),
         ICDint = lubridate::interval(birthday, entry_date),
         age.at.icd = ICDint/dyears(1)) %>%
  select(-ICDint)

range(pid.df$age.at.icd, na.rm= T)

pids.age= pid.df %>% filter(age.at.icd<110 & age.at.icd>39)%>% pull(pid)

#full pheno data

full_df= readRDS( file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/full_clinic_info.rds")
pharma= readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/data/processed_data/sept2022/Medication_df.rds")
echo_df= readRDS( "T:/fsa04/MED2-HF-Comorbidities/data/processed_data/sept2022/echo_measures_most_extreme.rds")
EF= readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/data/processed_data/sept2022/EF_df.rds")
NY= readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/data/processed_data/sept2022/NYHA.rds")
bnp= readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/data/processed_data/sept2022/labs_bnp.rds")

bnp = bnp %>% group_by(patient_id) %>%
  mutate(max.bnp= max(entry_value),
         med.bnp= median(entry_value))

EF= EF%>% group_by(patient_id) %>%
  mutate(min.ef= max(ef),
         med.ef= median(ef))

NY= NY %>% group_by(patient_id) %>%
  mutate(max.nyha= max(entry_value))

# get pids for different cohorts  ---------------------------------------------------
### create different cohorts for different features/filters and merge them at a final level

## 1) filter for ICD 10
# get pids for everybody that has HF

pids.HF = pid.df %>%
  distinct(pid) %>%
  pull(pid)

length(pids.HF)

pids.I42 = data %>%
  #filter(pid %in% patients)%>%
  group_by(pid) %>%
  #filter(any(icd3 == "I42")) %>%
  filter(any(icd4 %in% c(#"I42.0", DCM?
    "I42.1",
    "I42.2",
    #"I42.3", #eosinophil
    "I42.4",
    "I42.5",
    #"I42.6", #alcohol
    #"I42.7",# iatrogen. pharmaceutical
    "I42.8",#arvc
    "I42.9")# other
    ))%>%
  arrange(pid)%>%
  distinct(pid) %>%
  pull(pid)


# visualize
unique(data$pid)
ggVennDiagram(list("all_cause"= pids.HF,
                   "CM"= pids.I42,
                   "all_data" = unique(data$pid)))


# 443 patients in our cohort received an I42 diagnossis that is unwanted (example HCM)

# we remove those patients:
pids.HF.noCM= pids.HF[!pids.HF %in% pids.I42]

## 2) filter for BNP
pids.maxbnp = bnp %>%
  group_by(patient_id) %>%
  filter(med.bnp>120)%>%
  distinct(patient_id) %>%
  pull(patient_id) %>% as.integer()

### 3) filter for EF
pids.pEF = EF %>% filter(min.ef >= 50) %>% distinct(patient_id) %>%pull(patient_id)
pids.rEF = EF %>% filter(min.ef <= 40) %>% distinct(patient_id) %>%pull(patient_id)
pids.mrEF = EF %>% filter(min.ef >40 & min.ef <50) %>% distinct(patient_id) %>%pull(patient_id)

### 4) filter for NYHA
pids.nyha= NY %>% filter(!is.na(max.nyha)) %>%distinct(patient_id)%>%  pull(patient_id)

### Other echo data
pids.LVDD= echo_df %>% filter(edd.max>55) %>% pull(patient_id)
pids.LvEE= echo_df%>% filter(ee.max >15) %>%pull(patient_id)
pids.LVEA.min= echo_df %>% filter(ea.min<1 )%>% pull(patient_id)
pids.LVEA.max = echo_df %>% filter(ea.max>3 )%>% pull(patient_id)

pids.icd10 = as.integer(pids.HF.noCM)
pids.rEF= as.integer(pids.rEF)
pids.pEF= as.integer(pids.pEF)
pids.mrEF = as.integer(pids.mrEF)
pids.maxbnp= as.integer(pids.maxbnp)
pids.LvEE= as.integer(pids.LvEE)

#calculate number of hf diagnosis per patient ( cannot be the same day)
HF_allcause= c("I11.0", "I13.0", "I13.2", "I25.5", "I42.0","I42.3", "I42.6", "I42.7","I50.0", "I50.1", "I50.9")
data_count_HF= data %>%
  filter(pid %in% pids.icd10) %>%
  distinct(pid, entry_date, entry_value, icd3, icd4) %>%
  mutate(x= ifelse(icd3 %in% c("I50") | icd4 %in% HF_allcause, "yes", "no")) %>%
  filter(x== "yes") %>%
  arrange(pid) %>%
  distinct(pid, entry_date, x)%>%
  group_by(pid)%>%
  count(x)

hist(data_count_HF$n, breaks= 50)

#now pull
pid.icd10.2x= data_count_HF %>% filter(n>1) %>% pull(pid)

# add loop diuretic:
loop_diuretic= c("lasix", "torem","torasemid", "furosemid")
pids.loopd= pharma %>% filter(therapie.medikament %in% loop_diuretic) %>% distinct(patient_id) %>% pull(patient_id)



# hf= ICD (mandatory)+ either (ef, ee, or bnp) or 3x codes or pharma
pid.htx =readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/pids_endpoints.rds")

pid.HF.multi= data %>%
  distinct(pid) %>%
  filter(pid %in% pids.icd10,
         pid %in% pids.age,
         !pid %in% pid.htx$htx,
         pid %in% pids.maxbnp | pid %in% pids.nyha | pid %in% pids.LvEE | pid %in% c(pids.rEF,pids.mrEF) | pid %in% pid.icd10.2x | pid %in% pids.loopd)%>%
  pull(pid)

ggVennDiagram(list("loopd"= as.numeric(pids.loopd),
                   "bnp"= as.numeric(pids.maxbnp),
                   "2x"= as.numeric(pid.icd10.2x),
                   "ee"= as.numeric(pids.LvEE),
                   "icd10"= as.numeric(pids.icd10)))


# cohort definitions Jul 2022 --------------------------------------------------
# remove HTX patients (EF after HTX is wrongfully normal)

pid.HFpEF= intersect(pids.pEF, pid.HF.multi)
pid.HFpEF= pid.HFpEF[!pid.HFpEF %in% pid.htx$htx]

pid.HFrEF= intersect(pids.rEF, pid.HF.multi)
pid.HFrEF= pid.HFrEF[!pid.HFrEF %in% pid.htx$htx]

pid.HFmrEF= intersect(pids.mrEF, pid.HF.multi)
pid.HFmrEF= pid.HFmrEF[!pid.HFmrEF %in% pid.htx$htx]
#remove htx

hf.list= list("hf_all"= pid.HF.multi,
              "hfpef"= pid.HFpEF,
              "hfref"= pid.HFrEF,
              "hfmref"= pid.HFmrEF)

map(hf.list, length)
map(hf.list, function(x) length(intersect( x,unique(data$pid))))

saveRDS(hf.list, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds" )

x= data%>% group_by(pid) %>%
  distinct(pid, PheCode)%>% count()
x %>% filter(pid %in% hf.list$hfpef)%>% pull(n)%>% hist()
x %>% filter(pid %in% hf.list$hfref)%>% pull(n)%>% hist

map(hf.list, function(x){
  df= echo_df %>% filter(patient_id  %in% x)
  apply(df, 2,function(x) (median(x,na.rm = T)))

})

echo_df




# redo the cohort for a the less frequent visit cohort  ---------------------------------------------------

if(!exists("hf_meta")){
  hf_meta = read.csv("T:fsa04/MED2-HF-Comorbidities/data/RWH_March2020/levinson_comorbidities_in_hf_patients_2020-03-25.csv",
                     sep = ";",
                     na.strings=c("","NA")) %>% as_tibble
}

df = data %>% left_join(hf_meta, by = "pid")

HF_allcause= c("I11.0", "I13.0", "I13.2", "I25.5", "I42.0","I42.3", "I42.6", "I42.7","I50.0", "I50.1", "I50.9")

range(data$entry_date)

length(unique(data$pid))

data%>%
  group_by(pid)%>%
  mutate(hf = ifelse(any(icd3 %in% c("I50") | icd4 %in% HF_allcause), "HF", "noHF"))%>%
  distinct(pid, hf)%>%
  group_by(hf)%>%
  count()

               pids.HF = data %>%
  filter(icd3 %in% c("I50") | icd4 %in% HF_allcause) %>%
    distinct(pid) %>%
  pull(pid)%>% unique()

#calculate age at the first HF diagnosis:
age= df  %>%
  mutate(birthday= as_date(birthday),
         entry_date = as_date(entry_date),
         ICDint = interval(birthday, entry_date),
         age.at.icd = ICDint/dyears(1)) %>%
  select(-ICDint,-birthday) %>%
  filter(icd4 %in% HF_allcause) %>%
  group_by(pid) %>%
  summarise("age_at_HF"= min(age.at.icd))

df= df %>% left_join(age)

pids.age = df %>% filter(age_at_HF>39)%>% pull(pid)


pids.icd10= df %>%
  filter(#age_at_HF>40,
         pid %in% pids.HF
  ) %>%
  pull(pid)%>%
  unique()

unique(full_df$patient_id)

x= list("fullcode"= as.numeric(pids.icd10),
     "pheno"= as.numeric(unique(full_df$patient_id)),
     "patients"= as.numeric(patients))
ggVennDiagram(x)

pids.request= x$fullcode[!x$fullcode %in% x$patients]
pids.request= pids.request[!pids.request %in% pids.I42]
pids.request %>% write.csv("T:fsa04/MED2-HF-Comorbidities/data/PID_request_august2022.csv")
