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
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_network.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_pairwise_disease.R")

library(tidyverse)
library(cowplot)
library(ggrepel)
library(igraph)
library(ComplexHeatmap)
library(ggExtra)

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

#remove process data
data= data%>% filter(!startsWith(entry_value, "Z"))%>% drop_na(PheCode)

#calculate disease frequencies
phecode.frequency=
  disease_frequencies(pids, data %>%
                        filter(!icd4 %in% HF_allcause | icd3 != "I50" | icd3!= "I42"), "PheCode") %>%
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
length(phecodes)

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

data_r= data %>% distinct(pid, PheCode, Phenotype)

.links= create.links(pids= pids,
                     phecodes = phecodes,
                     data= data_r)

##
saveRDS(.links, file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_table_hf_cohort_fil.rds" )
.links= readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_table_hf_cohort_fil.rds")

saveRDS(list("links"= .links,
             "data"= data_r,
             "pid"= pids,
             "phecodes"= phecodes,
             "phe_dic" = Phe_dic),
        "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")

link.data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")



# construct HFnet ---------------------------------------------------------

net= quick_base_net(.links%>% filter(corr.phi>0, fisher.p.adj<0.0001),
                     pids = pids.list$hf_all,
                     data,
                     weight_col ="corr.phi",
                    Phe_dic = link.data$phe_dic)

V(net)$group_cat[V(net)$group_cat== "NULL"]= "injuries & poisonings"


# get adjacecny matrix
A= as.matrix(as_adjacency_matrix(net,attr = "corr.phi", type = "both"))

# calculate mean correlation of all edges with >0 correlation
b= rowSums(as.matrix(A)) / rowSums(as.matrix(A)> 0)

table(rownames(A)== names(b))

# scale every row by mean correlation
B= ((A) / b)

#sanity check
B[2,] ==A[2,]/ b[2]
isSymmetric(as.matrix(B))

# create HFnet by selecting max value of disease i,j and j,i scaled correalation
HFnet2= graph_from_adjacency_matrix(B, "max", weighted = T)
x=igraph::as_data_frame(HFnet2, "edges")
x1= x%>% rename("disease1"= from , "disease2"= to)
x2= x%>% rename("disease2"= from , "disease1"= to)

#add weight value to original link.table
links2= .links%>% left_join(x1, by = c("disease1", "disease2") )%>%
  left_join(x2, by = c("disease1", "disease2") )

.links$weight= links2%>%rowwise() %>% summarise(weight= max(weight.x, weight.y , na.rm= T))%>%
  pull(weight)

net= quick_base_net(.links%>% filter(is.finite(weight)),
                    pids = pids.list$hf_all,
                    data,
                    weight_col ="weight",
                    Phe_dic = link.data$phe_dic)

net
V(net)$group_cat[V(net)$group_cat== "NULL"]= "injuries & poisonings"

##save HFnet
saveRDS(net, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
saveRDS(.links, file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_table_hf_cohort_fil.rds" )
.links= readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_table_hf_cohort_fil.rds")





# plot pairwise stats -----------------------------------------------------------------

ggplot(.links, aes(x= log10(odds.ratio), y= corr.phi))+
  geom_point(size= 0.1)


ggplot(.links, aes(x= corr.tet, y= corr.phi))+
  geom_point(size= 0.1)


p1= ggplot(.links%>%filter(log10(odds.ratio)<10)%>%
             mutate(edge= ifelse(corr.tet>0 & fisher.p.adj<0.0001, "yes", "no")),
           aes(x= corr.phi, y= log10(odds.ratio), col = edge))+
  geom_point(size= 0.1)+
  scale_color_manual(values= cols.nice)+
  theme_minimal()+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)
  )+
  guides(colour = guide_legend(override.aes = list(size=5)))

ggMarginal(p1+
             theme(legend.position = "bottom")+
             labs(col= "Edge"),
           type="density", size=5, groupFill= F,margins = "both" )

p1= unify_axis(p1)

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/pairwise.corr.pdf",
    width= 5,
    height= 5)
print(ggMarginal(p1+
                   theme(legend.position = "bottom")+
                   labs(col= "Edge"),
                 type="density", size=5, groupFill= F,margins = "both" )
)
dev.off()

p2= ggplot(.links%>%filter(log10(odds.ratio)<10) ,
           aes(x= weight, y=corr.phi))+
  geom_point(size= 0.1)+
  #scale_color_manual(values= cols.nice)+
  theme_minimal()+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

unify_axis(p2)

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/pairwise_disease/pairwise.corr.weight.pdf",
    width= 3,
    height=3)
unify_axis(p2)
dev.off()

x= .links%>%filter(log10(odds.ratio)<10)%>%
  mutate(edge= ifelse(corr.tet>0 & fisher.p.adj<0.0001, "y", "n"))

prop.table(table(x$edge))

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



