
library(igraph)
library(tidyverse)
library(reshape2)
library(cowplot)

#process disgenet data: 
# library(GEOquery)
# gunzip(filename = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/all_gene_disease_associations.tsv.gz")
# gunzip(filename = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/disease_mappings.tsv.gz")

disgenet= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/all_gene_disease_associations.tsv", sep ='\t') %>% 
  as_tibble()

disdis= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/disease_mappings.tsv", sep ='\t') %>% 
  as_tibble()

unique(disdis$vocabulary)

disdis %>% filter(vocabulary == "ICD10CM")

disgenet %>% filter(diseaseId %in% c("C0235480", "C2585653", "C0694539")) %>% filter(score >0.3)
disgenet %>% filter(diseaseId %in% c("C0302358", "C0302360"))

#source("~/GitHub/RWH_analysis/scripts/network_scripts/create_networ/create_network_links.R")

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/ICD10_labeled_phe.rds") 
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/pidslist_v2021.rds") 
#identify the unique icd-codes to map :

#7000 values for whole data(30k pat)
unique(data$entry_value)


# compare coverage between disgenet and our icd10 data --------------------

# for 5k HF cohort only 4000 values
icd.hf= data %>% filter(pid %in% unlist(pids.list[1])) %>% pull(entry_value) %>% unique()
icd3.hf= data %>% filter(pid %in% unlist(pids.list[1])) %>% pull(icd3) %>% unique()
# filter for frequent diseases: 
icd3.hf= data %>% filter(pid %in% unlist(pids.list[1])) %>% distinct(pid, icd3) %>% count(icd3) %>% 
  #filter(n>50) %>%
  drop_na %>% 
  pull(icd3) %>% unique() 
dis_icd10 = disdis %>% filter(vocabulary == "ICD10CM") %>% pull(code)
dis_icd10cm= disdis %>% filter(vocabulary == "ICD10") %>% pull(code)

library("ggVennDiagram")
# 1
code_list= list("dis_icd10"= dis_icd10,
                "dis_icd10cm"= dis_icd10cm, 
                "data_hfcohort"= icd.hf)

ggVennDiagram(code_list,filename = NULL)

code_list= list("dis_icd10"= dis_icd10,
                "dis_icd10cm"= dis_icd10cm, 
                "data_hfcohort"= icd3.hf)

ggVennDiagram(code_list,filename = NULL)


## -> conclusion : 609 diseases are not mappable 
#lets take a look
problem_codes= icd3.hf[!icd3.hf %in% dis_icd10]

x= data%>% filter(icd3 %in% problem_codes) %>% distinct(icd3, Phenotype) %>% drop_na %>% arrange(icd3)
View(x)

########### step wise mapping 

# map the full length icd10 codes
hf_data= data %>% filter(pid %in% unlist(pids.list[1])) %>% distinct(entry_value, icd3, icd4)
disdis = disdis%>% filter(vocabulary == "ICD10") %>% select(code, diseaseId)

hf_data1= hf_data %>% left_join(., y = disdis  %>% rename(entry_value = code), by= "entry_value") %>%
  rename(diseaseId.full = diseaseId)

hf_data2= hf_data1 %>% left_join(., y = disdis %>% rename(icd4 = code), by ="icd4") %>% 
  rename(diseaseId.icd4 = diseaseId)

hf_data3= hf_data2 %>% left_join(., y = disdis %>% rename(icd3 = code), by ="icd3")%>% 
  rename(diseaseId.icd3 = diseaseId)

# combine mappings to merged column (checking each row for the first non NA value)
hf_data3 =hf_data3%>% mutate(disease_ID = coalesce(diseaseId.icd3, diseaseId.icd4, diseaseId.full))

hf_data3 %>% filter(is.na(disease_ID),
                    !is.na(icd3))

## i48 is for example still not mapped: 

# MESH:D001281

disdis %>% filter(vocabulary =="DO" ) %>% filter(code == "0060224")
disdis %>% filter(diseaseId == "C0004238")


#take a look at the more frequent diseases that fail to map
#use only icd3 based mappingS: 
hf_data %>% left_join(., y = disdis %>% rename(icd3 = code), by ="icd3")%>% 
  rename(diseaseId.icd3 = diseaseId)

mappable_icd3 = hf_data3 %>% distinct(icd3, disease_ID) %>% drop_na %>% pull(icd3) %>% unique

data %>%  filter(pid %in% unlist(pids.list[1])) %>% distinct(pid, icd3) %>%
  filter(!is.na(icd3)) %>% #group_by(pid)%>% 
  count(icd3) %>% arrange(desc(n)) %>% 
  mutate(mappable= ifelse(icd3 %in% mappable_icd3, "yes","no")) %>% 
  filter(n>50) %>%
  ggplot(.,aes(x= reorder(icd3, -n), y= n, alpha = mappable))+
  geom_col()


# map missing disease -----------------------------------------------------
# usually icd9 is covered for the problematic codes, map ICD10 to ICD9 to DisGeNet

library(touch)
?icd_map()
problems= hf_data3 %>% filter(is.na(disease_ID),
                              !is.na(icd3)) %>% 
  pull(entry_value) %>% unique


y= icd_map(problems, from = 10, to = 9 )
icd9 =as_tibble(cbind(icd10= problems, icd9 = y))

disdis %>% filter(vocabulary == "ICD9CM", 
                  code %in% icd9$icd9)
#nope that didnt help

problems

problem_codes



# another approach_ reduce mappings from disgenet side --------------------


#reduce mappings from disgenet side, cut all icd codes to three. 
##1 spell out icd10 ranges: 

#get the ranged mappings: 
test.df= disdis %>% 
  filter(vocabulary == "ICD10CM") %>% filter(grepl("-", code))
#split on minus:
x= strsplit(test.df$code, "-")
names(x)= test.df$diseaseId

# spell out the ranges by extracting the letter, start and endpoint
# then create single value for the whole sequence and map back to the disease id (UMLS-CUI)
spelled_out_range= map(x, function(y){
  # extract numbers:
  x1= str_extract(y[1], "([0-9]+).*$")
  x2= str_extract(y[2], "([0-9]+).*$")
  sequen = seq(x1,x2,1)
  letter= str_extract(y[1], "[A-Z]")
  sequen= map(sequen, function(xy){
    
    if (nchar(xy)==1){
      xy= paste0(0,xy)
      }
    return(xy)
    })
  paste0(letter, unlist(sequen))
})
#y= x$C0042721
# this will be added later to the dictionary
spelled_out_range= enframe(spelled_out_range, name = "diseaseId", value= "icd3") %>% unnest(icd3)


## create the main dictionary for icd3 codes (which have been collapsed from disgenet side)

dis3= disdis %>% 
  filter(vocabulary == "ICD10") %>% filter(!grepl("-", code)) %>%
  mutate(icd3= substr(code, 1, 3)) # str_replace(code, "(\\..).","")) 

# this mappping now contains for a single icd3 code multiple CUIs (including those that mapped to more
#specific daughter icd codes)

icd3_dic= hf_data3 %>% select(icd3, disease_ID) %>% drop_na %>% arrange(icd3) %>% rename(diseaseId= disease_ID)%>% 
  distinct(icd3, diseaseId)
# add those codes from the ranged map; 
icd3_dic= rbind(icd3_dic, spelled_out_range) %>% distinct(icd3, diseaseId)

#save dic 
saveRDS(icd3_dic,"T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/icd3_disgenet_dictionary.rds")

# create a relevant patient map 
patientmap_df= data %>% filter(pid %in% pids.list$hf_all) %>% distinct(pid, icd3, PheCode) %>% drop_na %>% left_join(icd3_dic)

disgenetmap_df= disgenet %>% filter(score>0.3) %>% select(geneSymbol, diseaseId, diseaseName, score)

pat_dis_gene = patientmap_df %>% left_join(disgenetmap_df, by= "diseaseId")

# most frequent disease where no gene is mapped: 
pat_dis_gene %>% filter(is.na(diseaseId)) %>% count(icd3) %>% arrange(desc(n)) %>% print(n= 100)

#maybe curate icd3 codes with over 100 patients manually? 

#~3746 genes are mapped
unique(pat_dis_gene$geneSymbol)

