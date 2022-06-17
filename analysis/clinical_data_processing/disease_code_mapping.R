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
directory= "T:/fsa04/MED2-HF-Comorbidities/data/RWH_March2020/"

#load phewas maps

# Standard WHO version, used here
phewas2= read_csv(file = paste0(directory,"data_Jan/phecode_icd10.csv")) %>% filter(!is.na(PheCode))

# Load Phecode definitions:
Phecode_definitions= read.csv(file.path(directory, "data_Jan/phecode_definitions1.2.csv"),
                              colClasses = rep("character", 7)) %>% as_tibble

# Load ICD10-dataset
icd10_labeled = readRDS("T:/fsa04/MED2-HF-Comorbidities/data/RWH_March2020/output_Jan/ICD10_labeled.rds")


#Load additional data for disease frequency calculation
#patient of interest
patients = readRDS(file.path(directory, "output_Jan/time_range_patientIDs.rds"))
#final_labeled data:
icd10 = readRDS("T:/fsa04/MED2-HF-Comorbidities/data/RWH_March2020/output_Jan/ICD10_labeled_phe.rds") %>%
  filter(pid %in% patients)
# calculate disease frequencies based on icd10 table, patient id, and ontology
source("~/GitHub/RWH_analysis/scripts/utils.R")

freqs = disease_frequencies(patients,icd10, feature = "entry_value") %>%
  select(-rel_freq) #%>%
  mutate(PheCode = as.character(PheCode))


# disease mapping of ICD10 codes to disease concepts ----------------------

# Disease Concept = PheWas codes

# explore how many icd10 codes are in the phewas bank
table(unique(icd10_labeled$entry_value) %in% phewas2$ICD10)

#add another icd10 code version XXX.X (icd4)
icd10_lab = icd10_labeled %>%
  mutate(icd4 =  str_replace(entry_value, "(\\..).","\\1"))  %>%
  distinct(entry_value, icd3, icd4) %>%
  arrange(entry_value)


#now perform a stepwise mapping, using different versions of the icd10 code to map to phewas
#1. map XXX.XX
icd10_1= icd10_lab %>%
  left_join(phewas2 %>%
              rename(entry_value = ICD10),
            by= "entry_value")

#2. of those XXX.XX that couldnt be mapped, take the XXX.X version and map again
icd10_2= icd10_1 %>%
  dplyr::filter(is.na(PheCode)) %>%
  select(1:3) %>%
  left_join(phewas2 %>%
              rename(icd4 = ICD10 ),
            by = "icd4")

#3. of those XXX.X that couldnt be mapped, take the XXX version and map again
icd10_3= icd10_2 %>%
  filter(is.na(PheCode)) %>%
  select(1:3) %>%
  left_join(phewas2 %>% rename(icd3 = ICD10 ), by = "icd3")

#4. combine all three mappings for greater coverage:
icd10_phewas= rbind(icd10_1 %>% drop_na(),
      icd10_2 %>% drop_na(),
      icd10_3 %>% drop_na() )

#still ~1000 XXX.XX codes couldnt be mapped to a PheWAS code
length(unique(icd10_labeled$entry_value))
table(unique(icd10_labeled$entry_value) %in% icd10_phewas$entry_value) #all


# analyze the number of mapped codes. step by step ------------------------
all_codes= unique(icd10_lab$entry_value)
length(all_codes) #7817 unique codes

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
table(unique(icd10_labeled$entry_value) %in% phewas2$ICD10)
table(unique(icd10_lab$icd4) %in% phewas2$ICD10)
table(unique(icd10_lab$icd3) %in% phewas2$ICD10)


# PROBLEM : muliple phecode for the same icd10 code -----------------------


## check for each level how many additional codes were mapped
dobled= phewas2%>%
  filter(ICD10 %in% phewas2$ICD10[duplicated(phewas2$ICD10)]) %>%
  arrange(ICD10)  %>%
  #distinct(ICD10,PheCode, Phenotype) %>%
  filter(ICD10 %in% all_codes)
dobled

length(unique(dobled$ICD10))
phewas2 %>% unique(ICD10)

phewas2[!unique(phewas2$ICD10),]%>% arrange(ICD10) %>% print(n=100)

table(all_codes %in% phewas2$ICD10)

dim(icd10_1 %>% filter(is.na(PheCode)))[1]

table(unique(icd10_lab$entry_value) %in% (icd10_1 %>% dplyr::filter(!is.na(PheCode)) %>% pull(entry_value)))# entry

level2y = dim(icd10_2%>% drop_na)[1]
level3y= dim(icd10_3 %>% drop_na)[1]
table(unique(icd10_labeled$entry_value) %in% icd10_3$entry_value)#icd4



# Continue analyzing the mapping:  ----------------------------------------


#find those code
x= icd10_labeled[!icd10_labeled$entry_value %in% icd10_phewas$entry_value,] %>%
  distinct(entry_value, icd3,label_icd3 ) %>%
  arrange(entry_value)

# find the frequency of those codes in our data to estimate whether those codes are important

x= x %>%
  left_join(freqs, by= "entry_value") %>%
  arrange(desc(freq)) %>%
  distinct(entry_value, freq, label_icd3)

#most codes are very infrequeutn, a few seem to be frequent
ggplot(data= x, aes(x= reorder(entry_value, freq), y= freq))+
  geom_col()

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


# Create actual mapping for icd10 table -----------------------------------

icd10 = readRDS(file = paste0(directory, "output_Jan/ICD10_labeled.rds"))

icd = icd %>%
  dplyr::select(pid,entry_date, entry_value) %>%
  left_join(icd_phewas, by= "entry_value")

saveRDS(icd, file =paste0(directory, "output_Jan/ICD10_labeled_phe.rds"))
icd_labeled_phe= readRDS(file =paste0(directory, "output_Jan/ICD10_labeled_phe.rds"))

# Add phenotype GROUPS ----------------------------------------------------
# Phenotypes can be grouped into higher level categories. This Information is taken
# from the Phecode_definition

table(unique(icd_labeled_phe$PheCode) %ni% Phecode_definitions$phecode)
# only 1 are not in the definition data. lets take a look at them
unique(icd_labeled_phe$PheCode)[unique(icd_labeled_phe$PheCode) %ni% Phecode_definitions$phecode] #only some NA..


icd_labeled_phe= icd_labeled_phe %>%
  left_join(Phecode_definitions %>%
              rename(PheCode = phecode)%>%
              select(-phenotype, -phecode_exclude_range))

saveRDS(icd_labeled_phe, file =paste0(directory, "output_Jan/ICD10_labeled_phe.rds"))

