### This script implements a diagnostic algorithm to select HF cohorts from the RWH
# in the first part an inclusive data frame of nyha , bnp and echo data is generated
# then this data and icd10 data is used to generate different patient cohorts
# i.e. all_hf , hfpef, hfref, hfmref
# finally those cohorts are used to plot different characteristics


# libs and data -----------------------------------------------------------
library(ggpubr)
library(cowplot)
library(ggVennDiagram)
library(tidyverse)
library(lubridate)

directory= "T:/fsa04/MED2-HF-Comorbidities/data/processed_data/"
# full data frame:
icd10= readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/ICD10_labeled_phe.rds")

hf = read.csv("T:/fsa04/MED2-HF-Comorbidities/data/RWH_March2020/levinson_comorbidities_in_hf_patients_2020-03-25.csv",
              sep = ";",
              na.strings=c("","NA")) %>% as_tibble

# full cohort (28k with at least 3 visits in more than 2 years)
patients =readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/time_range_patientIDs.rds")

#endpoints:
endpoints= readRDS(file= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/pids_endpoints.rds")

## get single data values
medis= readRDS(file.path(directory, "Medication_df.rds"))

ef= readRDS( file = file.path(directory, "EF_df.rds"))

nyha = readRDS(file = file.path(directory, "NYHA.rds") )

pharma= readRDS(file = file.path(directory, "Medication_df.rds"))

nyha.max = nyha %>%
  group_by(patient_id) %>%
  summarise(nyha.max= max(entry_value))

echo_df= readRDS( file = file.path(directory, "echo_measures_most_extreme.rds") )

labs.bnp = readRDS(file = file.path(directory, "labs_bnp.rds") )

labs.bnp.median = labs.bnp %>%
  group_by(patient_id) %>%
  summarise(median.bnp= median(entry_value),
            max.bnp= max(entry_value)
  )

bp.mean = readRDS(file.path(directory, "diag_bp.rds"))


# merge additional data to full matrix ------------------------------------


#### FULL MERGE
full_df = echo_df %>% full_join(nyha.max) %>% filter(ea.min<25)%>%
  full_join(labs.bnp.median) %>%
  full_join(bp.mean %>% select(-count)) %>%
  filter(patient_id %in% patients) # control check

range(full_df$mean.sys[!is.na(full_df$mean.sys)])

#add sex info to this data

full_df= full_df %>% #mutate(patient_id = as.integer(patient_id)) %>%
  left_join(hf %>% rename(patient_id = pid) %>%
                        select(patient_id, sex),by = "patient_id")

saveRDS(full_df, file = file.path(directory, "full_clinic_info.rds"))



# get pids for different cohorts  ---------------------------------------------------
### create different cohorts for different features/filters and merge them at a final level

## 1) filter for ICD 10
# get pids for everybody that has HF
HF_allcause= c("I11.0", "I13.0", "I13.2", "I25.5", "I42.0", "I42.5", "I42.8", "I42.9",
               "I50.0", "I50.1", "I50.9")

pids.HF = icd10 %>%
  filter(pid %in% patients) %>%
  filter(icd3 %in% c("I50") | icd4 %in% HF_allcause) %>%
  #filter(icd3 %in% c("I50") | icd4 %in% c("I11.0","I13.0", "I13.2","I25.5", "I42.")) %>% ### RE EVALUATE IF I25 should be included!! (wihtout hf.pids. shrink to 9.500)
  distinct(pid) %>%
  pull(pid)

length(pids.HF)

pids.I42 = icd10 %>%
  filter(pid %in% patients)%>%
  group_by(pid) %>%
  #filter(any(icd3 == "I42")) %>%
  filter(any(icd4 %in% c(#"I42.0",
                         "I42.1",
                         "I42.2",
                         #"I42.3",
                         "I42.4",
                         "I42.5",
                         #"I42.6",
                         #"I42.7",
                         "I42.8",
                         "I42.9"
                         )
             )) %>%
  arrange(pid)%>%
  distinct(pid) %>%
  pull(pid)

length(pids.I42)

# visualize
ggVennDiagram(list("all_cause"= pids.HF,
                   "CM"= pids.I42))
# 443 patients in our cohort received an I42 diagnossis that is unwanted (example HCM)

# we remove those patients:
pids.HF.noCM= pids.HF[!pids.HF %in% pids.I42]

saveRDS(pids.HF.noCM, file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/pids_HF_noCM.rds")

## 2) filter for BNP
pids.maxbnp = full_df %>%
  group_by(patient_id) %>%
  filter(max.bnp>120)%>%
  distinct(patient_id) %>%
  pull(patient_id) %>% as.integer()

### 3) filter for EF
pids.pEF = full_df %>% filter(ef.min >= 50) %>% distinct(patient_id) %>%pull(patient_id)
pids.rEF = full_df %>% filter(ef.min <= 40) %>% distinct(patient_id) %>%pull(patient_id)
pids.mrEF = full_df %>% filter(ef.min >40 & ef.min <50) %>% distinct(patient_id) %>%pull(patient_id)

### 4) filter for NYHA
pids.nyha= full_df %>% filter(!is.na(nyha.max)) %>%distinct(patient_id)%>%  pull(patient_id)

### Other echo data
pids.LVDD= full_df %>% filter(edd.max>55) %>% pull(patient_id)
pids.LvEE= full_df%>% filter(ee.max >15) %>%pull(patient_id)
pids.LVEA.min= full_df %>% filter(ea.min<1 )%>% pull(patient_id)
pids.LVEA.max = full_df %>% filter(ea.max>3 )%>% pull(patient_id)


# plot venn diagrams for different cohorts --------------------------------------------------------------

# hf general =
hf.all= list("Diagnosis_code_HF"= pids.HF.noCM,
         "BNP>125 "=pids.maxbnp,
         "EF>50"= pids.pEF,
         "EF<40"= pids.rEF)
ggVennDiagram::ggVennDiagram(x = hf.all, color= "black")
map(hf.all, class)
# compare with i42 paients.
hf.i42= list(Icd= pids.HF.noCM,
             "I42"=pids.I42,
             "EF>50"= pids.pEF,
             "EF<40"= pids.rEF)
ggVennDiagram(hf.i42, color= "black")


# COHORT DEFINITION HFPEF, HFREF, HFMREF ------------------------------------------
pids.icd10= pids.HF.noCM

ggVennDiagram(list(nyha= as.integer(pids.nyha), icd10 =pids.icd10, i42= pids.I42))
# all patients with a nyha record are captured by the icd10 definition, therefore
# nyha is not needed

pids.icd10 = as.integer(pids.icd10)
pids.rEF= as.integer(pids.rEF)
pids.pEF= as.integer(pids.pEF)
pids.mrEF = as.integer(pids.mrEF)
pids.maxbnp= as.integer(pids.maxbnp)

### HFrEF :
# 1) ICD10 code for HF +
# 2) lowest EF below 40 recorded
pidsHFREF = pids.icd10[pids.icd10 %in% pids.rEF]

### HFmrEF
# 1) ICD10 code for HF +
# 2) lowest EF below 40 recorded
pidsHFMREF= pids.icd10[pids.icd10 %in% pids.mrEF]

### HFpEF
# 1) ICD10 code for HF +
# 2) lowest EF above or equal 50 recorded
pidsHFPEF = pids.icd10[pids.icd10 %in% pids.pEF]
pidsHFPEF = pidsHFPEF[pidsHFPEF %in% pids.maxbnp]

pids.list= list(hf_all= pids.icd10, hfref= pidsHFREF, hfpef= pidsHFPEF, hfmref= pidsHFMREF)

map(pids.list, length)

saveRDS(pids.list,
        file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids_fil.rds" )

#remove htx patients:

pids.list2= map(pids.list, function(x){x[!x %in% htx$htx]})

map(pids.list, length)
map(pids.list2, length)


saveRDS(pids.list2,
        file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids_filhtx.rds" )

# COHORT DEFINITION HF vs CT ----------------------------------------------
patients =readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/time_range_patientIDs.rds")
HF_allcause= c("I11.0", "I13.0", "I13.2", "I25.5", "I42.0", "I42.5", "I42.8", "I42.9", "I50.0", "I50.1", "I50.9")
hf.pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/pidslist_v2021.rds")


ct.pid= patients[!patients %in% hf.pids.list$hf_all]
hf.pid= hf.pids.list$hf_all

#refine the ct cohort by removing patients with hf code and nyha

# double check hf diagnosis
if(!exists("pheno_df")){
  pheno_df= readRDS("T://fsa04/MED2-HF-Comorbidities/data/processed_data/full_clinic_info.rds")
}

data2= data %>% left_join(pheno_df %>% dplyr::rename(pid= patient_id),by= "pid")

pid.ct.2= data2 %>% filter(pid %in% ct.pid,
                           !icd4 %in% HF_allcause & !icd3 %in% c("I50", "I42"),
                           is.na(nyha.max)) %>%
  pull(pid) %>% unique()

pids.hf_ct= list("ct"= pid.ct.2, "hf"= hf.pid)
map(pids.hf_ct, length)

saveRDS(pids.hf_ct, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/pidslist_hf_ct.rds")



# plot summary tables of features for cohorts -----------------------------

df.ct_hf= get_summary_table(data = data,pids.list =  pids.hf_ct)
df.hfpef_hfref= get_summary_table(data = data,pids.list =  hf.pids.list[2:3])

df= df.ct_hf

table.df= df %>%
  filter(patient_cohort != "none") %>%
  distinct(pid,
           patient_cohort,
           sex,
           age.at.mean,
           median.BMI,
           mean.sys,
           mean.dias,
           median.LDL,
           median.HDL,
           median.Trigs,
           median.Chol,
           median.gfrcg,
           median.hba1c,
           median.bnp,
           PheCode_count,
           icd10gm_count,
           charlson_score,
           elixhauser_wscore,
           intu,
           htx,
           defi,
           ef.min,
           edd.max,
           ea.min,
           ea.max,
           ee.max)

tab.hf= table.df %>% select(-pid) %>%
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
         median.gfrcg,
         median.hba1c,
         median.bnp,
         PheCode_count,
         icd10gm_count,
         charlson_score,
         elixhauser_wscore,
         intu,
         htx,
         defi

  )%>%
  tbl_summary(by = patient_cohort) %>%
  add_p()%>%
  add_overall()%>%
  add_n() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Main Cohort**")

tab.hf%>%    # build gtsummary table
  as_gt()%>%
  gt::gtsave(             # save table as image
    filename = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/manuscript/main/cohort_description_hf.png"
  )
## no hfpef vs hfref
df= df.hfpef_hfref

table.df= df %>%
  filter(patient_cohort != "none") %>%
  distinct(pid,
           patient_cohort,
           sex,
           age.at.mean,
           nyha.max,
           median.BMI,
           mean.sys,
           mean.dias,
           median.LDL,
           median.HDL,
           median.Trigs,
           median.Chol,
           median.gfrcg,
           median.hba1c,
           median.bnp,
           PheCode_count,
           icd10gm_count,
           charlson_score,
           elixhauser_wscore,
           intu,
           htx,
           defi,
           ef.min,
           edd.max,
           ea.min,
           ea.max,
           ee.max)

tab.hfpef= table.df %>% select(-pid) %>%
  select(patient_cohort,
         sex,
         age.at.mean,
         nyha.max,
         median.BMI,
         mean.sys,
         mean.dias,
         median.LDL,
         median.HDL,
         median.Trigs,
         median.Chol,
         median.gfrcg,
         median.hba1c,
         median.bnp,
         PheCode_count,
         icd10gm_count,
         charlson_score,
         elixhauser_wscore,
         intu,
         htx,
         defi,ef.min,
         edd.max,
         ea.min,
         ea.max,
         ee.max

  )%>%
  tbl_summary(by = patient_cohort) %>%
  add_p()%>%
  add_overall()%>%
  add_n() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**HF sub cohorts**")


class(tab.hfpef)
tab.hfpef=
  as_gt(tab.hfpef)

tab.hfpef%>%    # build gtsummary table
  gt::gtsave(             # save table as image
    filename = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/manuscript/main/cohort_description_hfpef.png"
  )
# Tips
# pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/manuscript/main/cohort_description_hfpef.pdf")
#
# tab.hfpef
# dev.off()


dis.freq= map(pids.hf_ct, function(x){
  d = disease_frequencies(pids=x , icd= data)
  d= d %>% left_join(Phe_dic)
})

dis.freq$hf

# use hfpef,hfref, hfnocm to plot clinical data ---------------------------

###  1) add age to the data frame (age at HF diagnosis)

age= icd10  %>%
  inner_join(hf %>% select(pid, birthday),by= "pid") %>%
  mutate(birthday= as_date(birthday),
         entry_date = as_date(entry_date),
         ICDint = interval(birthday, entry_date),
         age.at.icd = ICDint/dyears(1)) %>%
  select(-ICDint,-birthday)

age= age %>% #filter(pid %in% pids.HF.noCM) %>%
  filter(icd3 %in% c("I50", "I25") | icd4 %in% c("I11.0","I13.0", "I13.2","I25.5")) %>%
  group_by(pid) %>%
  summarise(min(age.at.icd))

# create a full dataframe that inlcudes the cohort labels:

ef_full_df = full_df %>%
  filter(patient_id %in% pids.list$hf_all) %>%
  mutate(HF= "hf_no_ef")%>%
    mutate(HF= ifelse(patient_id %in% pids.list$hfpef, "hfpef", HF),
           HF= ifelse(patient_id %in% pids.list$hfref, "hfref", HF),
           HF= ifelse(patient_id %in% pids.list$hfmref, "hfmref", HF))%>%
  mutate(HF= factor(HF, levels= c("hf_no_ef","hfpef", "hfmref", "hfref"))) %>%
     filter(sex %in% c("f", "m"))

ef_full_df= ef_full_df %>% left_join(age %>% rename(patient_id= pid,
                                        min.age.at.icd10 = `min(age.at.icd)`), by= "patient_id") %>%
  filter(min.age.at.icd10 >18)

### 2) Plotting:

cols= c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51")

gender= ggplot(ef_full_df, aes(x= HF, fill= sex))+
  geom_histogram(stat="count")+
  theme_minimal()+
  scale_fill_manual(values= cols[4:5])+
  labs(x= "")

gender

nyha= ef_full_df %>% filter(!is.na(nyha.max)) %>%  count(HF, nyha.max) %>%
  ggplot(., aes(x= HF, y= n,fill = nyha.max))+
  geom_col(aes(x = HF, y = n, fill = nyha.max), position = "fill", na.rm= T)


bnp= ef_full_df %>%
  ggplot(., aes(x= HF, y= log(labs.bnp.max), fill = HF))+
  geom_jitter(alpha= 0.1)+
  geom_violin(alpha= 0.9)+
  theme_minimal()+
  scale_fill_manual(values= cols)+
  geom_hline(yintercept = log(125))+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "hfref") +
  labs(x= "")
bnp

bp= ef_full_df %>% pivot_longer(cols = c(mean.sys, mean.dias)) %>%
  ggplot(., aes(x= HF,y= value, fill = name))+
  scale_fill_manual(values= cols[4:5])+
  geom_jitter(alpha= 0.1)+
  geom_violin(alpha= 0.9)+
  labs(x= "", y= "mean_blood_pressure", fill ="")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  labs(x= "")

bp

edd= ef_full_df %>%
  ggplot(., aes(x= HF, y= edd.max, fill = HF))+
  geom_jitter(alpha= 0.1)+
  scale_fill_manual(values= cols)+
  geom_violin(alpha= 0.9)+
  theme_minimal()+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "hfref") +
  labs(x= "")
edd

ages= ef_full_df %>%
  ggplot(., aes(x= HF, y= min.age.at.icd10, fill = HF))+
  geom_jitter(alpha= 0.1)+
  scale_fill_manual(values= cols)+
  geom_violin(alpha= 0.9)+
  theme_minimal()+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "hfref") +
  labs(x= "")


ages

p1 = plot_grid(gender,
               ages+theme(legend.position = "none"),
               nyha+theme(legend.position = "none"),
               bnp+theme(legend.position = "none"),
               edd+theme(legend.position = "none"),
               bp,
               labels = "AUTO",
               ncol= 3)
p1

pdf(file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/clinical_characterization_HFcohorts.pdf",
    width= 12,
    height= 9)
p1
dev.off()


# explore pharma data -----------------------------------------------------

top_medis = pharma %>% filter(patient_id %in% pids.list$hf_all) %>%
  distinct(patient_id, entry_date, therapie.medikament ) %>%
  count(therapie.medikament) %>%
  arrange(desc(n)) %>% print(n= 100)

top_medis$therapie.medikament[grepl("unat", top_medis$therapie.medikament)]

loop_diuretic= c("lasix", "torem","torasemid", "furosemid")

pharma2 = pharma %>% filter(patient_id %in% pids.list$hf_all) %>% filter(therapie.medikament %in% loop_diuretic)


#  distinct(patient_id, entry_date, therapie.medikament )

unique(pharma$therapie.medikament)

pids.hfpef_pot= pids.icd10[pids.icd10 %in% pids.pEF]

# get patients with pEF and diagnosis, but those that did not have elevated BNP
x= pharma %>% filter(patient_id %in% pids.hfpef_pot & !patient_id %in% pids.list$hfpef) %>%
  filter(therapie.medikament %in% loop_diuretic) %>% pull(patient_id)

map(pids.list, function(y){x %in% y})

pids.loopd= pharma %>% filter(therapie.medikament %in% loop_diuretic) %>% distinct(patient_id) %>% pull(patient_id)


hf.all= list("Diagnosis_code_HF"= pids.HF.noCM,
             "EF>50 "=pids.pEF,
             "loopD"= pids.loopd,
             "EF<40"= pids.rEF)
ggVennDiagram(hf.all, color= "black")

