## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-05-13
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
##  plot supp figure with
## ---------------------------
library(tidyverse)
library(ComplexHeatmap)
link.data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")

.links = link.data$links
hfpef = .links %>%
  filter(fisher.p.adj<0.01, (pcorr.corpor)>0.0) %>%
  filter(grepl("hfpef", disease2)) %>% arrange(desc(pcorr.corpor)) %>%
  dplyr::select(disease1, dis1_phenotype, dis2_phenotype, everything()) %>% pull(dis1_phenotype)
hfref = .links %>%
  filter(fisher.p.adj<0.01, (pcorr.corpor)>0.0) %>%
  filter(grepl("hfref", disease2)) %>% arrange(desc(pcorr.corpor)) %>%
  dplyr::select(disease1, dis1_phenotype, dis2_phenotype, everything()) %>% pull(dis1_phenotype)

ggvenn::ggvenn(data = list("hfpef"= hfpef,"HFreF"=  hfref))

both= unlist(intersect(hfref,hfpef))
hfref= unique(hfref[!hfref %in% both])
hfpef= unique(hfpef[!hfpef %in% both])

length(both)/ length(c(both, hfpef, hfref))

vec1= rep("hfpef",length(hfpef))
names(vec1)= hfpef
vec2= rep("hfpef+hfref",length(both))
names(vec2)= both
vec3= rep("hfref",length(hfref))
names(vec3)= hfref


df.hflinks= unique(c(both, hfref,hfpef))

heatmap.df = .links %>% filter(dis1_phenotype %in% df.hflinks, grepl("hf", dis2_phenotype))%>%
  dplyr::select(dis1_phenotype, dis2_phenotype,pcorr.corpor)%>% ungroup()%>%
  pivot_wider(names_from = dis2_phenotype, values_from = pcorr.corpor)

mat= heatmap.df%>%  as.data.frame()%>% column_to_rownames("dis1_phenotype")
rows= c(vec1,vec2,vec3)
h.links= Heatmap(mat[names(rows), ], name = "partial rho", row_split =rows,cluster_columns = F,
                 row_names_side  = "right",
                 row_names_max_width =unit(50, "cm"))
h.links
pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/hf_neighbors.pdf",
    width= 8,
    height= 10)
h.links
dev.off()

