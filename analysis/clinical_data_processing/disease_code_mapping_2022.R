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
##  This script maps the icd10 codes from the RDWH to PheWas codes and creates a dictionary
## for this purpose
##
## ---------------------------


library(tidyverse)
library(lubridate)
library(ICD10gm)

source("analysis/utils/utils.R")

directory= "T:/fsa04/MED2-HF-Comorbidities/data/RWH_September2022/raw/"

#load phewas maps

# Standard WHO version, used here
phewas2= read_csv(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/databases/PheWas/phecode_icd10.csv") %>% filter(!is.na(PheCode))

# Load Phecode definitions:
Phecode_definitions= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/databases/PheWas/phecode_definitions1.2.csv",
                              colClasses = rep("character", 7)) %>% as_tibble

## icd10 new
raw= read.csv(file= "T:/fsa04/MED2-HF-Comorbidities/data/RWH_September2022/raw/levinson_comorbidities_data_2022-10-06.csv",sep = ";")

## demo_new
data= read.csv(paste0(directory,"levinson_comorbidities_pids_2022-10-06.csv"), sep = ";")

# Load ICD10-dataset#old df
icd10_lab_old = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled.rds")

data = as_tibble(raw) %>% filter(entity %in% c("diagnostik.nebendiagnose.icd.kode","diagnostik.hauptdiagnose.icd.kode"))

## assess which disease codes are new=
table(unique(data$entry_value) %in% unique(icd10_lab_old$entry_value))

old_phe_dic = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/icd10_phewas_dictionary.rds")
## assess which disease codes are new

# add ICD10 labels --------------------------------------------------------
#disease Concept = PheWas codes
# explore how many icd10 codes are in the phewas bank
table(unique(data$entry_value) %in% phewas2$ICD10)


#add another icd10 code version XXX.X (icd4)
data2= data %>%
  mutate(icd4 =  str_replace(entry_value, "(\\..).","\\1"),
         icd3 = substr(entry_value,1,3 )) %>%
  distinct(entry_value, pid, entity, icd3, icd4)

table(unique(data2$entry_value) %in% unique(old_phe_dic$entry_value))
table(unique(data2$icd4) %in% unique(old_phe_dic$icd4))
table(unique(data2$icd3) %in% unique(old_phe_dic$icd3))

icd10_lab = data2
#phewas2= phewas2.mod

#now perform a stepwise mapping, using different versions of the icd10 code to map to phewas
#1. map XXX.XX
icd10_1= icd10_lab %>%
  left_join(phewas2 %>%
              rename(entry_value = ICD10),
            by= "entry_value")

#2. of those XXX.XX that couldnt be mapped, take the XXX.X version and map again
icd10_2 = icd10_1 %>%
  dplyr::filter(is.na(PheCode)) %>%
  select(1:5) %>%
  left_join(phewas2 %>%
              rename(icd4 = ICD10 ),
            by = "icd4")

#3. of those XXX.X that couldnt be mapped, take the XXX version and map again
icd10_3 = icd10_2 %>%
  filter(is.na(PheCode)) %>%
  select(1:5) %>%
  left_join(phewas2 %>% rename(icd3 = ICD10 ), by = "icd3")

#4. combine all three mappings for greater coverage:
icd10_phewas = rbind(icd10_1 %>% drop_na(PheCode),
      icd10_2 %>% drop_na(PheCode),
      icd10_3 %>% drop_na(PheCode) )


#still ~1000 XXX.XX codes couldnt be mapped to a PheWAS code
length(unique(icd10_lab$entry_value))
table(unique(icd10_lab$entry_value) %in% icd10_phewas$entry_value) #all



### add disease categories:

icd10_phewas= icd10_phewas %>%
  left_join(Phecode_definitions %>%
              rename(PheCode = phecode)%>%
              select(-phenotype, -phecode_exclude_range))


#add entry_date back

icd10_phewas= icd10_phewas%>% full_join(data%>% select(-entry_time),
                          by= c("pid", "entity", "entry_value")
                          )
saveRDS(icd10_phewas, file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")




# check mapping stats -----------------------------------------------------

unmapped_codes= sort(unique(data$entry_value[!data$entry_value %in% icd10_phewas$entry_value]))

x = data%>%
  filter(entry_value %in% unmapped_codes)%>%
  distinct(pid, entry_value, icd4, icd3)%>%
  group_by(entry_value)%>%
  count %>%
  arrange(desc(n))%>%
  filter(!grepl("Z", entry_value),
         !grepl("U", entry_value),
         !grepl("Y", entry_value),
         !grepl("Z", entry_value))

top_codes= x %>% filter(n>100)%>% pull(entry_value)

filter(ICD10 %in%  top_codes)

p1= x %>%
  filter(n>100)%>%
  ggplot(., aes(x= reorder(entry_value,n), y= n))+
  geom_col()+
  coord_flip()+
  scale_y_log10()+
  labs(x= "unmapped ICD10 codes")+
  theme_minimal()

p1 %>% unify_axis()

x %>% mutate(PheCode = ifelse(entry_value== "I27.28", "415.2", ""),
             PheCode = ifelse(entry_value== "F05.0", "290.2", PheCode),
             PheCode = ifelse(entry_value== "F32.9", "296.22", PheCode),
             PheCode = ifelse(entry_value== "F32.1", "296.22", PheCode),
             PheCode = ifelse(entry_value== "M51.2", "722.6", PheCode),
             PheCode = ifelse(entry_value== "L82", "702.2", PheCode))
#m51,
df= rbind(c("M51", "722.6"), # lumbar degen
      c("M81", "743.11"), #
      c("F32", "296" ), ##depressive mood
      c("F05", "290.2" ), #delir
      c("I27.28", "415.2"), #pulmonary dis
      c("M15.9","740.2") #polyarthorsis
      )
colnames(df)= c("ICD10", "PheCode")

maps= as_tibble(df)%>%
  left_join(phewas2 %>% select(-ICD10), by= "PheCode")%>%
  mutate(`ICD10 String` = NA)%>%
  distinct()

phewas2.mod= rbind(phewas2, maps)

saveRDS(phewas2.mod, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/databases/PheWas/phecode_icd10_mod.rds")

phe_dic = icd10_phewas%>% distinct(entry_value, icd4, icd3, PheCode, Phenotype, category)

# plot the # codes -------------------------------------------------------------------------

barplot(c("ICD10"= length(unique(icd10_phewas$entry_value)),
       "ICD10_4digit"= length(unique(icd10_phewas$icd4)),
       "ICD10_3digit"= length(unique(icd10_phewas$icd3)),
  "PheCode"= length(unique(icd10_phewas$PheCode))
),ylim = c(0,7000)
)
icd10_phewas %>% distinct(pid, entry_value, PheCode, icd3)





# analyze the number of mapped codes. step by step ------------------------
all_codes= unique(data$entry_value)
length(all_codes) #7517 unique codes

#1
codes= icd10_1 %>% filter(is.na(PheCode))%>% distinct(entry_value) %>% pull(entry_value)
length(unique(codes)) #not mapped: 3050 with xxx.xx
#dim(icd10_1)[1]-length(unique(codes))
xcodes= icd10_1 %>% filter(!is.na(PheCode)) %>% pull(entry_value)
length(unique(xcodes)) #mapped: 4787 with xxx.xx

#2
codes2= icd10_2 %>% filter(is.na(PheCode))%>% pull(entry_value)
length(unique(codes2)) #not mapped: 1170 with xxx.x
xcodes2= icd10_2 %>% filter(!is.na(PheCode)) %>% pull(entry_value)
length(unique(xcodes2)) #mapped: 1860 with xxx.xx

#3
codes3= icd10_3 %>% filter(is.na(PheCode))%>% pull(entry_value)
length(unique(codes3)) #not mapped: 783 with xxx
xcodes3= icd10_3 %>% filter(!is.na(PheCode)) %>% pull(entry_value)
length(unique(xcodes3)) #mapped: 410 with xxx.xx



#check if this provides more coverage
table(unique(icd10_lab$entry_value) %in% phewas2$ICD10)
table(unique(icd10_lab$icd4) %in% phewas2$ICD10)
table(unique(icd10_lab$icd3) %in% phewas2$ICD10)













# PROBLEM : muliple phecode for the same icd10 code -----------------------

## check for each level how many additional codes were mapped
dobled= phewas2%>%
  filter(ICD10 %in% phewas2$ICD10[duplicated(phewas2$ICD10)]) %>%
  arrange(ICD10)  %>%
  #distinct(ICD10,PheCode, Phenotype) %>%
  filter(ICD10 %in% all_codes)%>%
  filter(!grepl("Z", ICD10),
         !grepl("U", ICD10),
         !grepl("Y", ICD10),
         !grepl("Z", ICD10))
dobled

unique(dobled$ICD10)


icd10_phewas%>% filter(entry_value=="A15.9")

# Continue analyzing the mapping:  ----------------------------------------

table(icd10_lab$entry_value %in% icd10_phewas$entry_value)

unmapped_codes=
  unique(icd10_lab$entry_value[!icd10_lab$entry_value %in% icd10_phewas$entry_value])

#most codes are very infrequeutn, a few seem to be frequent
data%>% filter(entry_value %in% unmapped_codes)%>%
  ggplot(., aes(x= entry_value))+
  geom_histogram(stat= "count")

unmapped_codes_freq= data %>%
  distinct(pid, entry_value) %>%
  filter(entry_value %in% unmapped_codes)%>%
  group_by(entry_value)%>%
  count %>%
  filter(n>10)%>%
  pull(entry_value)

data%>% filter(entry_value %in% unmapped_codes_freq)%>%
  ggplot(., aes(x= entry_value))+
  geom_histogram(stat= "count")


#zoom in on those that are frequent.  J91 might be an importatn one, the others seem to bel less relevant
ggplot(data= x[1:50,] , aes(x= reorder(entry_value, freq), y= freq))+
  geom_col()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle =60, hjust= 1))

x %>% print
#save dictionary:
saveRDS(icd10_phewas, file = paste0(directory, "output_Jan/icd10_phewas_dictionary.rds"))


# how many unique features does every code contain? plot for overvie
df_feature= sapply(colnames(icd10_phewas), function(x){
  length(unique(icd10_phewas[[x]]))
})

data.frame( names= names(df_feature), nfeat= df_feature) %>%
  as_tibble() %>%
  ggplot(aes(x= reorder(names, nfeat), y= nfeat))+
  geom_col()


