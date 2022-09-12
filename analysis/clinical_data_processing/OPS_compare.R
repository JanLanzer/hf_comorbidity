## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-09-12
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## Analyze OPS data
## ---------------------------


library(tidyverse)

ops= readRDS(file= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ops_data.rds")

icd = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe.rds")
#pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids.rds")

icd = icd %>% left_join(phedic %>% distinct(PheCode, Phenotype))

phedic= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/icd10_phewas_dictionary.rds")

phecodes= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/top300_disease.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids.rds")
#pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids_filhtx.rds")
dis.l= map(pids.list,length)

ops= ops%>%   mutate(hf.type = ifelse(pid %in% pids.list$hfref,
                                            "hfref",
                                            ifelse(pid %in% pids.list$hfpef,
                                                   "hfpef",
                                                   ifelse(pid %in% pids.list$hfmref,
                                                          "hfmref",
                                                          "unlabeled"))))


# calc freqs by Klassentitel
df2= ops %>%
  distinct(pid,entry_value , Klassentitel, hf.type)%>%
  group_by(hf.type,Klassentitel)%>%
  count()%>%
  ungroup()

df2 = df2%>% mutate(prop= ifelse(hf.type == "hfmref", n/length(pids.list$hfmref),
                           ifelse(hf.type=="hfpef", n/length(pids.list$hfpef),
                                  ifelse(hf.type=="hfref", n/length(pids.list$hfref),
                                         n/length(pids.list$hf_all)))))


df3= df2 %>%
  pivot_wider(names_from = hf.type, values_from = prop,-n)%>%
  arrange(desc(hfref))#%>% left_join(ops%>% distinct(entry_value, Klassentitel))

df3 %>% mutate(prop.diff = hfpef-hfref)%>%
  arrange(desc(prop.diff))%>%
  print(n=100)
df3 %>% mutate(prop.diff = hfpef-hfref)%>%
  arrange((prop.diff))


df3


df1 %>%
  group_by(hf.type)%>%
  arrange(desc(n))


df1= ops %>%
  distinct(pid, title.4, hf.type)%>%
  group_by(hf.type, title.4)%>%
  count()


df1 %>% arrange(desc(n)) %>%
  ungroup()%>%
  group_by(hf.type)%>%
  count()

