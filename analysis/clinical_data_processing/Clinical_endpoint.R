# load OPS data and annotate with metafile from :https://www.dimdi.de/dynamic/de/klassifikationen/downloads/?dir=ops/
# start query OPS and identify potential clinical endpoints and number of patients behind those features.

library(tidyverse)
library(purrr)
library(broom)
library(mosaic)
library(reshape2)
library(stringr)


# Load data and create a map of OPS codes from the RDH  ---------------------------------------------
ops = read.csv("T:/fsa04/MED2-HF-Comorbidities/data/RWH_March2020/levinson_comorbidities_in_hf_ops_2020-03-25.csv",
               sep = ";",
               na.strings = c("",NA)) %>% as_tibble


directory= "T:/fsa04/MED2-HF-Comorbidities/"



icd10 = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds") %>%
  filter(pid %in% patients)
# get pids from hf cohort:
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")

map(pids.list, length)

patients= pids.list$hf_all

ops_dictionary = read.table(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/ops2020syst_kodes.txt", sep = ";") %>%
  as_tibble() %>%
  rename(entry_value = V7,
         Klassentitel =V9,
         title.4= V10,
         title.5= V11,
         title.6= V12,
         chapter = V4 ,
         code3 = V5,
         code4= V6)



ops_dictionary2=ops_dictionary %>%
    select(chapter, code3, code4,entry_value, Klassentitel, title.4, title.5, title.6)

#1.
ops_red = ops %>%
  distinct(entry_value) %>%
  mutate(entry_value_3 = substr(entry_value, 1, 5) )

ops_red = ops_red %>% left_join(ops_dictionary2)



ops_red2 = ops_red %>%
  dplyr::filter(is.na(chapter )) %>%
  distinct(entry_value, entry_value_3) %>%
  left_join(ops_dictionary2 %>%
              rename(entry_value_3= entry_value) ,
            by= "entry_value_3") %>%
  distinct()

## OPS_dic can now be joined via entry_value. some codes remain NA.
OPS_dic= ops_red %>% drop_na %>% rbind(ops_red2)


# calculate frequencies of ops -------------------------------

# add descriptions
ops= ops %>% left_join(OPS_dic)


# calculate frequencies of ops in patient cohort

  ops_pids = ops %>%
    filter(pid %in% patients) %>%
    select(pid, entry_value)%>%
    unique()%>% # important, icd3 codes are duplicated among single patients(this is only a summary of all following codes)
    drop_na

  # calculate frequencies with table-function (this is only based on the icd3 code)
  ops_freq= as_tibble(as.data.frame(sort(table(ops_pids %>% pull(entry_value)))))
  colnames(ops_freq) = c("entry_value","freq")

  #     add relative frequency and description
ops_freq= ops_freq %>%
    mutate(rel_freq = freq / length(unique(patients))) %>%
    arrange(desc(freq)) %>%
    left_join(OPS_dic)

ops_freq  %>% filter(code3 == "5-35" & code4 == "5-37") %>% print(n=100)

ops_freq  %>% filter(entry_value_3 == "5-375") %>% print(n=100) ## 89 patients with heart transplant


# identify htx patients ---------------------------------------------------

#via ops
ops_htx = ops %>%  filter(entry_value_3 == "5-375")

#via icd
icd_htx= icd10 %>% filter(icd4 == "Z94.1")
length(unique(icd_htx$pid))

#via phecode
phecode_htx= icd10 %>% filter(PheCode == "429.1")
length(unique(phecode_htx$pid))

htx = unique(c(ops_htx$pid, icd_htx$pid, phecode_htx$pid))
table(htx %in% unlist(pids.list[2:4]))
length(htx)

table(htx %in% patients)
#440 htx patients..


# identify defibrillator placement patients -------------------------------

ops_defi = ops %>%  filter(entry_value_3%in% c("5-378", "5-377") )%>% print(n=200)
defi = unique(ops_defi$pid)
#remove patients with htx from general OPS code
defi =  defi[!defi %in% htx]


table(defi %in% unlist(pids.list[2:4]))

# identify coronary intevention patients ----------------------------------
ops_PCI = ops %>%  filter(code4 == "8-83") %>% print(n=200)
PCI = unique(ops_PCI$pid)

table(PCI %in% unlist(pids.list[2:4]))


table(ops_PCI %>% filter(pid %in% pids.list$hfref)%>% distinct(pid, entry_value_3)%>% pull(entry_value_3))
# intubation -----------------------------------------------------

ops_intu= ops %>% filter(entry_value_3 %in% c("8-701","8-704","8-706", "8-852") )
unique(ops_mort$title.4)

intu= ops_intu %>%
  pull(unique(pid))

table(intu %in% unlist(pids.list[2:4]))

clinical_endpoint= list("intubation"= intu,
     "defi"= defi,
     "htx" = htx,
     "pci"= PCI)

saveRDS(clinical_endpoint, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/pids_endpoints.rds")



