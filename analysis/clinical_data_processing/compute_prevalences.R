## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-11-04
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## compute comorbidity prevalences
## ---------------------------


source("analysis/utils/utils.R")

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")
phecodes= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/topfreq_disease.rds")



data= data %>%
  mutate(hf = ifelse(pid   %in% pids.list$hfref,
                     "HFrEF",
                     ifelse(pid   %in% pids.list$hfpef,
                            "HFpEF",
                            ifelse(pid   %in% pids.list$hfmref,
                                   "HFmrEF",
                                   "unlabeled"))))%>%
  mutate(hf= factor(hf, levels= c("HFpEF","HFmrEF", "HFrEF", "unlabeled")))


data %>% filter(hf %in% c("HFpEF", "HFrEF"))

pef= disease_frequencies(pids.list$hfpef, data)
ref= disease_frequencies(pids.list$hfref, data)

df.= ref %>% rename(ref_freq  = "rel_freq")%>%
  left_join(pef %>% rename(pef_freq  = "rel_freq"),
            by= c("PheCode", "Phenotype"))%>% mutate(diff= ref_freq-pef_freq)
df. %>% print(n=100)

