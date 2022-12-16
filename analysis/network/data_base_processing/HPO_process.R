## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-12-15
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##  HPO ontology is used to to generat
##
## ---------------------------

library(corpustools)
library(tidyverse)
library(ontologyIndex)
library(ontologySimilarity)

hpo = read.delim("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/HPO/phenotype.hpoa", skip = 4)%>% as_tibble()

icd_pheno_map= read.delim("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/databases/Phenotype_mapping/semiautomatic_ICD-pheno.txt")%>%
  as_tibble()

#icdmaps= read.delim("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/HPO/HPO_to_ICD10.txt") %>% as_tibble()

phemaps= read.delim("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/HPO/HPO_to_phecode.txt")%>% as_tibble()

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
Phe_dic = data %>% distinct(entry_value  , icd3, icd4, PheCode)

link.data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")
# Subset HPO to disease net ----------------------------------------------------------

phecodes =unique(c(link.data$links$disease1, link.data$links$disease2))

relevant_diseases= Phe_dic %>% filter(PheCode %in% phecodes)

mapped_codes= relevant_diseases%>%
  left_join(icd_pheno_map %>% dplyr::rename(entry_value= ICD10.Code)) %>% print(n=100)

mapped_codes_filt= mapped_codes%>%
  distinct(PheCode, HPO.ID)%>%
  drop_na()#%>%
  filter(HPO.ID %in% hpo$id)

class(hpo)
information_content <- descendants_IC(hpo)
?descendants_IC
dis_phe_sets= split( mapped_codes_filt$HPO.ID, mapped_codes_filt$PheCode)

sim_mat <- get_sim_grid(ontology=hpo,
                        term_sets=dis_phe_sets,
                        information_content = information_content,
                        term_sim_method = "lin")

longformat= sim_mat %>%
  as.data.frame%>%
  rownames_to_column("nodeA") %>%
  pivot_longer(cols = -nodeA , names_to = "nodeB", values_to = "sim")%>%
  distinct(nodeA, nodeB, sim)

full_hpo_net = graph_from_data_frame(longformat, directed = F)

E(full_hpo_net)$weight= longformat$sim

full_hpo_net  = simplify(full_hpo_net, remove.multiple = T, remove.loops = T)

hpo_net= backbone_filter(full_hpo_net, alpha = 0.05)

saveRDS(hpo_net,"T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/hpo_net.rds" )


# FULL hpo net ------------------------------------------------------------


mapped_codes= data %>% distinct(PheCode, Phenotype, entry_value)%>% drop_na%>%
  left_join(icd_pheno_map %>% rename(entry_value= ICD10.Code)) %>%
  filter(HPO.ID %in% hpo$id)

information_content <- descendants_IC(hpo)

dis_phe_sets= split( mapped_codes$HPO.ID, mapped_codes$PheCode)

sim_mat <- get_sim_grid(ontology=hpo,
                        term_sets=dis_phe_sets,
                        information_content = information_content,
                        term_sim_method = "lin")

longformat= sim_mat %>%
  as.data.frame%>%
  rownames_to_column("nodeA") %>%
  pivot_longer(cols = -nodeA , names_to = "nodeB", values_to = "sim")%>%
  distinct(nodeA, nodeB, sim)

full_hpo_net = graph_from_data_frame(longformat, directed = F)
E(full_hpo_net)$weight= longformat$sim
full_hpo_net  = simplify(full_hpo_net, remove.multiple = T, remove.loops = T)

hpo_net= backbone_filter(full_hpo_net, alpha = 0.05)

saveRDS(hpo_net,"T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/hpo_net_full.rds" )
