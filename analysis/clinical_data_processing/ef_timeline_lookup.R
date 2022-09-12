## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-06-30
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## check for heart transplant
## ---------------------------
# -------------------------------------------------------------------------


htx= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/pids_endpoints.rds")
hftype = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids.rds")
icd10= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe.rds")

icd10 = icd10 %>% filter(pid %in% htx$htx)



ef= readRDS("T:/fsa04/MED2-HF-Comorbidities/data/processed_data/EF_df.rds")

ef = ef %>% rename(entry_date_ef= entry_date)

comb.df= ef %>% full_join(icd10%>% rename(patient_id= pid), by= "patient_id")

comb.df2 = comb.df %>%
  filter(patient_id %in% htx$htx) %>% filter(icd3 == "Z94")%>%
  distinct(entry_date_ef, ef, entry_date, icd3)

comb.df2 = comb.df2 %>% group_by(patient_id)%>%
  mutate(min.htx = min(entry_date),
         min.ef= min(entry_date_ef))    %>%
  mutate(delta.t = (min.htx     %--% min.ef )/dmonths(1) )

df= comb.df2 %>% distinct(delta.t, patient_id)

hist(df%>% filter(patient_id %in% hftype$hfpef)%>% pull(delta.t))
hist(df%>% filter(patient_id %in% hftype$hfref)%>% pull(delta.t))
df= comb.df2 %>% group_by(patient_id) %>% mutate(min.t= min(delta.t))
hist(df$min.t)

df2= df%>%distinct(min.t, patient_id)%>% filter(patient_id %in% hftype$hfpef)
hist(df2$min.t)

### get patients

pre.htx.ef= comb.df2 %>% distinct(patient_id, entry_date_ef, ef, min.htx  , delta.t)%>%
  group_by(patient_id)%>%
  mutate(min.ef = min(ef)) %>%
  filter(delta.t<0)%>%
  mutate(hf= ifelse(min.ef>=50 , "hfpef", "hfref"))%>% distinct(patient_id, hf)

table(pre.htx.ef$hf)

## remove data from HTX patients after to HTX

icd10%>% filter(pid %in% htx$htx)


