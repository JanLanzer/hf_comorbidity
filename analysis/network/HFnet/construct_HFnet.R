## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-05-05
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## -create HF comorbidity network (HFnet)
## -run louvain cluster
## -save the network in various forms
## ---------------------------


source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_network.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_pairwise_disease.R")


library(tidyverse)
library(cowplot)
library(ggrepel)
library(igraph)
library(ComplexHeatmap)

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")

map(pids.list, length)

# PheCode dictionary:
Phe_dic= data %>%
  distinct(PheCode, Phenotype,category) %>%
  drop_na

# define which PheCodes will be used in the analysis ----------------------

#remove HF-codes
HF_allcause= c("I11.0", "I13.0", "I13.2", "I25.5", "I42.0", "I42.5", "I42.8", "I42.9", "I50.0", "I50.1", "I50.9")
#get all pids
pids= unlist(pids.list$hf_all)
#pids= unlist(pids.list[2:4])
pids= unlist(c(pids.list$hfref,pids.list$hfpef))
data= data%>% filter(!startsWith(entry_value, "Z"))

#calculate disease frequencies
phecode.frequency=
  disease_frequencies(pids, data %>% filter(!icd4 %in% HF_allcause | icd3 != "I50" | icd3!= "I42"), "PheCode") %>%
  arrange(desc(rel_freq)) %>%
  #top_n(topn_disease) %>%
  mutate(PheCode2 = as.character(PheCode)) %>%
  filter(!PheCode2 %in% c("428.1", "428.2"))%>%
  left_join(Phe_dic)

# at least 10 patients with that PheCode
cutoff= 50

phecodes= phecode.frequency%>%
  filter(freq>cutoff)%>%
  pull(PheCode2)

p.phe.counts=
  ggplot(phecode.frequency, aes(x= reorder(PheCode, -freq)  , y= log10(freq)))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept = log10(cutoff), col = "grey", size= 1, type= 2)+
  geom_vline(xintercept = length(phecodes),  col = "grey")+
  theme(axis.text.x = element_blank())+
  labs(x= "PheCodes",
       y= "log10 counts")


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/feature_cutoff.pdf",
    height = 3,
    width = 2)
p.phe.counts
dev.off()

saveRDS(phecodes, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/topfreq_disease.rds")


## function to add hfref and hfpef as diagnosis codes to the data

get_pid_df= function(pids, name2){
  enframe(pids, value= "pid") %>% mutate(PheCode= name2) %>% select(-name)

}

hfref= get_pid_df(pids.list$hfref, "hfref")
hfpef= get_pid_df(pids.list$hfpef, "hfpef")
#hfmref= get_pid_df(pids.list$hfmref, "hfmref")

#append data
data_grouped = data %>%
  distinct(pid,PheCode) %>%
  drop_na

data_comb= rbind(hfpef, hfref,# hfmref,
                 data_grouped)

#update phedic:
Phe_dic= rbind(Phe_dic,
               c("hfref", "hfref", "circulatory system"),
               c("hfpef", "hfpef", "circulatory system"))
               #c("hfmref", "hfmref", "circulatory system")


data_comb= data_comb %>% left_join(Phe_dic)
phecodes= c(phecodes, "hfref", "hfpef")#ä, "hfmref")

.links= create.links(pids= pids,
                     phecodes = phecodes,
                     data= data_comb)

##
saveRDS(.links, file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_table_hf_cohort_fil.rds" )
.links= readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_table_hf_cohort_fil.rds")

saveRDS(list("links"= .links,
             "data"= data_comb,
             "pid"= pids,
             "phecodes350"= phecodes,
             "phe_dic" = Phe_dic),
        "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")

# check hf nodes: ---------------------------------------------------------
#
.links %>%
  filter(rho_part_lower >0.00 & fisher.p.adj<1) %>%
  filter(grepl("hfpef", disease2)) %>% arrange(desc(pcorr.corpor)) %>%
  dplyr::select(disease1, dis1_phenotype, dis2_phenotype,pcorr.corpor, everything())%>%
  print(n=199)# %>% pull(disease1)


.links %>%
  filter( fisher.p.adj<0.05, odds.ratio>1) %>%
  filter(grepl("hfpef", disease2)) %>% arrange(desc(pcorr.corpor)) %>%
  dplyr::select(disease1, dis1_phenotype, dis2_phenotype, everything()) %>% print(n=199)# %>% pull(disease1)


.links %>%
  filter(rho_part_lower >0) %>%
  filter(grepl("hfref", disease2)) %>% arrange(desc(pcorr.corpor)) %>%
  dplyr::select(disease1, dis1_phenotype, dis2_phenotype, everything())%>% print(n=199) #%>% pull(disease1)

 .links %>%
  filter(fisher.p.adj<0.05, odds.ratio>1) %>%
  filter(grepl("hfref", disease2)) %>% arrange(desc(pcorr.corpor)) %>%
  dplyr::select(disease1, dis1_phenotype, dis2_phenotype, everything()) %>% print(n=199)#%>% pull(disease1)

# hfmref = .links %>%
#   filter(rho_part_lower >0, fisher.p.adj<0.05) %>%
#   filter(grepl("hfmref", disease2)) %>% arrange(desc(pcorr.corpor)) %>%
#   dplyr::select(disease1, dis1_phenotype, dis2_phenotype, everything())# %>% pull(disease1)

# cluster HFnet -----------------------------------------------------------------
hist(.links$corr.tet)
 hist(.links$pcorr.corpor)

ggplot(.links, aes(x= pcorr.corpor, y= corr.tet))+
  geom_point(size= 0.1)

ggplot(.links %>%
         mutate(edge= ifelse(rho_part_lower>0 , "y", "n")), aes(x= pcorr.corpor, y= corr.tet, col = edge))+
  geom_point(size= 0.1)



 ggplot(.links%>%
         mutate(edge= ifelse((corr.tet>0 & fisher.p.adj<0.05), "y", "n")), aes(x= pcorr.corpor, y= corr.tet,col = edge))+
  geom_point(size= 0.1)

 p1= ggplot(.links%>%
              mutate(edge= ifelse((corr.tet>0 & fisher.p.adj<0.05), "y", "n")), aes(x= pcorr.corpor, y= corr.tet))+
   geom_point(size= 0.1)
p1+ geom_density_2d_filled(alpha = 0.5)+
  geom_vline(xintercept = 0.01)+
  geom_hline(yintercept = 0.1)

ggplot(.links%>% filter(rho_part_lower>0 & fisher.p.adj<0.05), aes(x= pcorr.corpor, y= corr.tet))+
  geom_point(size= 0.1)


#define cut link cut off  and create HFnet

hist(.links$rho_part_lower)

net= quick_base_net(.links%>% filter(rho_part_lower>0.00, fisher.p.adj<1),#, fisher.p.adj<0.05),
                pids,
               data_comb,
               weight_col ="pcorr.corpor")

net= quick_base_net(.links%>% filter(corr.tet>0, fisher.p.adj<0.05),
                     pids,
                     data_comb,
                     weight_col ="corr.tet")


#clsuter HFnet with louvain

set.seed(2)
net.mod= modularize(net, method= "spinglass")
net.mod= modularize(net, method= "louvain")
net.mod= modularize(net, method= "leiden",resolution_parameter = 1)

Hfnet= net.mod$fullnet
net.mod$h.map
links= igraph::as_data_frame(Hfnet, "edges")
nodes= igraph::as_data_frame(Hfnet, "vertices")
nodes %>% filter(name %in% c("hfpef", "hfref", "hfmref"))


nodes%>% filter(group_louv ==1
              )%>% as_tibble()%>% arrange(desc(size))%>% print(n=100)

saveRDS(net.mod$fullnet, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
Hfnet= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
write_delim(links, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFnet_links.tsv", delim = "\t")
write_delim(nodes, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFnet_nodes.tsv",delim = "\t")

library(WriteXLS)
n.list= split(nodes, nodes$group_louv)
saveRDS(n.list, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFnet_nodes.rds")

WriteXLS(n.list, SheetNames = names(n.list), ExcelFileName = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFnet_nodes.xlsx")
