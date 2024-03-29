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
## --------------------------
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

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")
phecodes= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/topfreq_disease.rds")

map(map(pids.list, unique), length)
#8062+6585+3018
#pids.list$tavi= pids.oi

# data.r= data %>% filter(pid %in% c(pids.list$hf_all),#,pids.list$hfpef, pids.list$hfref),
#                 PheCode %in% phecodes)
pids.list$hf_unl= pids.list$hf_all[!pids.list$hf_all %in% unlist(pids.list[2:4])]

data.r= data %>% filter(pid %in% c(pids.list$hfmref,pids.list$hfpef, pids.list$hfref,pids.list$hf_unl ),
                        PheCode %in% phecodes)

df= get_summary_table(data = data.r, pids.list = pids.list[2:5])
#df= get_summary_table(data = data.r, pids.list = pids.list[c(2,3,4,5)])

length((table.df$patient_cohort))


table.df %>% distinct(pid, patient_cohort) %>% group_by(patient_cohort)%>% count()

table.df= df %>%
  filter(patient_cohort != "none") %>%
  distinct(pid,
           age.at.icd,
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
           #icd10gm_count,
           charlson_score,
           charlson_wscore,
           elixhauser_wscore,
           elixhauser_score,
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
           #age.at.mean,
           sex)

table.df$patient_cohort= str_replace_all(table.df$patient_cohort, "hfpef", "HFpEF")
table.df$patient_cohort= str_replace_all(table.df$patient_cohort, "hfref", "HFrEF")
table.df$patient_cohort= str_replace_all(table.df$patient_cohort, "hfmref", "HFmrEF")

table(table.df$patient_cohort)
table(df$patient_cohort)

table.df %>%
  saveRDS(., file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/patient_metadata_2022.rds")

table.df= readRDS( file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/patient_metadata_2022.rds")



gt.tab= table.df %>%
  filter(patient_cohort %in% c("HFpEF", "HFrEF"))%>%
  #filter(sex.y != "u",
  #       patient_cohort != "hf_all")%>%
  #mutate(patient_cohort= factor(patient_cohort, levels= c("HFrEF", "HFmrEF", "HFpEF")))%>%
  select(patient_cohort,
         sex,
         age.at.icd,
         median.BMI,
         mean.sys,
         mean.dias,
         median.LDL,
         median.HDL,
         median.Trigs,
         median.Chol,
         #median.gfrcg,
         #median.hba1c,
         median.bnp,
         PheCode_count,
         #icd10gm_count,
         #charlson_score,
         #charlson_wscore,
         #elixhauser_wscore,
         elixhauser_score,
         intu,
         #htx,
         defi,
         pci,
         #ef.min,
         #edd.max,
         #ea.min,
         #ea.max,
         #ee.max

  )%>%
  mutate("logBNP"= log10(median.bnp))%>%
  select(-median.bnp)%>%
  dplyr::rename(Sex= sex,
                "Age (y)" = age.at.icd,
                BMI = median.BMI,
                #"EF(%)"= ef.min,
                "Systolic RR (mmHg)"= mean.sys,
                "Diastolic RR (mmHg)"= mean.dias,
                "LDL (mg/dl)"= median.LDL,
                "HDL (mg/dl)"= median.HDL,
                "TC (mg/dl)"= median.Trigs,
                "Chol (mg/dl)"= median.Chol,
                "log(NT-BNP)"= logBNP,
                "n(PheCodes)"= PheCode_count,
                "Intubated"= intu,
                "Elixhauser"= elixhauser_score,
                #Charlson Index= charlson_score,
                #HTX= htx,
                "PCI"= pci,
                "ICD implantation"= defi
                )  %>%
  #mutate(Elixhauser = Elixhauser +0000000.1)%>%
  tbl_summary(by = patient_cohort ,
              type= list(Elixhauser ~ "continuous"),
              statistic = list(Elixhauser ~ " {mean} ({sd})")) %>%
  #add_p()%>%
  add_p(test = is_numeric() ~ "t.test")%>%
  #add_difference(by= patient_cohort)%>%
  add_overall()%>%
  add_n() %>%
  modify_header(label ~ "**Variable**") %>%
  #modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**HF subtypes**")
  modify_spanning_header(c("stat_1", "stat_2") ~ "**HF subtypes**")
gt.tab

table(gt.tab$patient_cohort)
gt::gtsave(as_gt(gt.tab), file ="output/cohort_summary.ltx")
gt::gtsave(as_gt(gt.tab), file ="output/cohort_summary_2groups.png")




df.test= table.df%>% filter(patient_cohort %in% c("HFpEF", "HFrEF"))

rstatix::t_test(formula= elixhauser_score  ~ patient_cohort,
       data = df.test)


df.test%>% group_by(patient_cohort)%>% summarise(mean(median.BMI, na.rm= T))

t_test(formula= median.BMI ~ patient_cohort,
       df.test)

# full cohort -------------------------------------------------------------

data.r= data %>% filter(pid %in% c(pids.list$hf_all),
                        PheCode %in% phecodes)

df= get_summary_table(data.r, pids.list)
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
  #mutate(patient_cohort= ifelse(patient_cohort=="hfpef", "HFpEF", "HFrEF"))%>%
  select(patient_cohort,
         sex,
         age.at.icd,
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
         nyha.max,
         #edd.max,
         #ea.min,
         #ea.max,
         #ee.max

  )%>%
  dplyr::rename(Sex= sex,
                Age = age.at.icd,
                BMI= median.BMI,
                "EF(%)"= ef.min,
                "Systolic RR"= mean.sys,
                "Diastolic RR"= mean.dias,
                "LDL(mg/dl)"= median.LDL,
                "HDL(mg/dl)"= median.HDL,
                "TC(mg/dl)"= median.Trigs,
                "Chol(mg/dl)"= median.Chol,
                "n(PheCodes)"= PheCode_count,
                Intubated= intu,
                Elixhauser= elixhauser_score,
                #HTX= htx,
                PCI= pci,
                ICD_implant= defi,
                "NYHA"= nyha.max,
  )  %>%
  tbl_summary(by = patient_cohort) %>%
  add_p()%>%
  add_overall()%>%
  add_n() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Main Cohort**")


gt.tab

gt::gtsave(data= gt.tab, filename= ".html", path= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/cohorttable.html")


# save for manual process -------------------------------------------------

df.= table.df %>%
  #select(-pid) %>%
  #mutate(patient_cohort= ifelse(patient_cohort=="hfpef", "HFpEF", "HFrEF"))%>%
  #filter(sex.y != "u",
  #       patient_cohort != "hf_all")%>%
  mutate(patient_cohort= factor(patient_cohort, levels= c("HFrEF", "HFmrEF", "HFpEF")))%>%
  select(pid,
         patient_cohort,
         sex,
         age.at.icd,
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
         charlson_score,
         charlson_wscore,
         #elixhauser_wscore,
         elixhauser_score,
         intu,
         #htx,
         defi,
         pci,
         #ef.min,
         #edd.max,
         #ea.min,
         #ea.max,
         #ee.max

  )%>%
  dplyr::rename(Sex= sex,
                "Age (y)" = age.at.icd,
                BMI = median.BMI,
                #"EF(%)"= ef.min,
                "Systolic RR (mmHg)"= mean.sys,
                "Diastolic RR (mmHg)"= mean.dias,
                "LDL (mg/dl)"= median.LDL,
                "HDL (mg/dl)"= median.HDL,
                "TC (mg/dl)"= median.Trigs,
                "Chol (mg/dl)"= median.Chol,
                "n(PheCodes)"= PheCode_count,
                "Intubated"= intu,
                "Elixhauser Score"= elixhauser_score,
                #Charlson Index= charlson_score,
                #HTX= htx,
                "PCI"= pci,
                "ICD implantation"= defi
  )

df. %>%
