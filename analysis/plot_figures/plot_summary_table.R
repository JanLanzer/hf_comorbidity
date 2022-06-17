## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-05-10
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ------------https://www.danieldsjoberg.com/gtsummary/index.html---------------
##
## Notes:
##
## plot summary stats for hfpef, hfref cohort
## ---------------------------
library(rlang)
library(gtsummary)
library(gtable)
library(tidyverse)

source("analysis/utils/utils.R")
data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe.rds")

pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/pidslist_oct2021.rds")
phecodes= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/top300_disease.rds")

data.r= data %>% filter(pid %in% c(pids.list$hfpef, pids.list$hfref),
                PheCode %in% phecodes)

df= get_summary_table(data.r, pids.list[2:3])
table.df= df %>%
  filter(patient_cohort != "none") %>%
  distinct(pid,
           #age_at_HF,
           median.BMI,
           mean.sys,
           mean.dias,
           median.LDL,
           median.HDL,
           median.Chol,
           median.Trigs,
           intu,
           htx,
           defi,
           pci,
           PheCode_count,
           icd10gm_count,
           charlson_score,
           elixhauser_wscore,
           nyha.max,
           median.bnp,
           patient_cohort,
           edd.max,
           ea.min,
           ea.max,
           ee.max,
           ef.min,
           median.gfrcg,
           median.hba1c,
           age.at.mean,
           sex)

gt.tab= table.df %>%
  select(-pid) %>%
  mutate(patient_cohort= ifelse(patient_cohort=="hfpef", "HFpEF", "HFrEF"))%>%
  select(patient_cohort,
         sex,
         age.at.mean,
         median.BMI,
         mean.sys,
         mean.dias,
         median.LDL,
         median.HDL,
         median.Trigs,
         median.Chol,
         #median.gfrcg,
         #median.hba1c,
         #median.bnp,
         PheCode_count,
         #icd10gm_count,
         #charlson_score,
         elixhauser_wscore,
         intu,
         htx,
         defi,
         pci,
         ef.min,
         #edd.max,
         #ea.min,
         #ea.max,
         #ee.max

  )%>%
  dplyr::rename(Sex= sex,
                Age = age.at.mean,
                BMI= median.BMI,
                Systolic_RR= mean.sys,
                Diastolic_RR= mean.dias,
                "LDL(md/dl)"= median.LDL,
                "HDL(md/dl)"= median.HDL,
                "TC(md/dl)"= median.Trigs,
                "Chol(mg/dl)"= median.Chol,
                "n(PheCodes)"= PheCode_count,
                Intubated= intu,
                Elixhauser= elixhauser_wscore,
                HTX= htx,
                PCI= pci,
                ICD_implant= defi,
                "EF(%)"= ef.min)  %>%
  tbl_summary(by = patient_cohort) %>%
  add_p()%>%
  add_overall()%>%
  add_n() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Main Cohort**")


gt::gtsave(as_gt(gt.tab), file ="output/cohor_summary.png")
