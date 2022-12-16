# ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-12-08
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## compare HF net with morbinet
## ---------------------------
library(tidyverse)
library(ComplexHeatmap)
library(igraph)
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_hetnet.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_network.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_delta_con.R")

.links= readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_table_hf_cohort_fil.rds")
link.data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")
hf_net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
edgelist= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds")
Phedic= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/icd10_phewas_dictionary.rds")
hpo_net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/hpo_net.rds" )
hpo_net_full = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/hpo_net_full.rds" )
data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
# we compare  -------------------------------------------------------------

# 1 HF_net (normal)
# 2 HF_net (OR cut off)
# 3 hpo net
# 4 hpo net subset to hfnet nodes
# 5 morbinet
# 6 morbinet subset to hfnet nodes

# 1 HF_net the way it is implemented
hf_real = hf_net


# 2 HF_net with OR
hf_or = table_to_links( link.data$links ,
                     p.val = 0.001,
                     weight_dd = 2,
                     weight_col = "odds.ratio")

hf_or= simplify_link_table(hf_or) %>%
  graph_from_data_frame(., directed = F)


#3) full hpo net
hpo_net_full
hpo_net_full =  simplify_link_table(igraph::as_data_frame(hpo_net_full))%>%
  graph_from_data_frame(., directed = F)

# 4) reduced hpo net
hpo_net= induced_subgraph(hpo_net_full, vids = V(hpo_net_full)$name %in%link.data$phecodes)


#5) full morbinet
morbinet= process_morbinet()
morbinet= simplify_link_table(morbinet)
morbinet= graph_from_data_frame(morbinet, directed = F)

#6) reduced morbinet
morbinet_r= induced_subgraph(
  morbinet,
  V(morbinet)$name %in% V(hf_real)$name)

#7) reduced hf_net to morbinet
hf_real_r= induced_subgraph(
  hf_real,
  V(hf_real)$name %in% V(morbinet_r)$name)

hf_real_r2= induced_subgraph(
  hf_real,
  V(hf_real)$name %in% V(hpo_net)$name)

#8) check high confidence versions
hfnet_hc= table_to_links( .links ,
                          p.val = 0.000001,
                          weight_dd = 2,
                          weight_col = "odds.ratio")


# 9 overlap of all three nodes

vec1= phecodes[phecodes %in% V(hpo_net_full)$name]
vec2= vec1[vec1 %in% V(morbinet)$name]


# compare sizes  ----------------------------------------------------------

net.list= list("HFnet"= hf_real,
         #"HFnet_"= hf_or,
         "HFnet_subset_M"= hf_real_r,
         "HFnet_subset_HPO"= hf_real_r2,
          "HPOnet"= hpo_net_full,
          "HPOnet_subset_HF"= hpo_net,
         "Morbinet"= morbinet,
         "Morbinet_subset_HF"= morbinet_r)

net.sizes= sapply(net.list, function(x){
  tibble("nodes"= length(V(x)),
         "edges"= length(E(x))
  )
})


p.network.sizes=
  net.sizes%>% as.data.frame() %>%rownames_to_column("feature") %>%
  pivot_longer(-feature)%>% mutate(value = unlist(value))%>%
  ggplot(., aes(x= name, y= value))+
  facet_grid(vars(cols= feature),  scales="free")+
  #scale_y_log10()+
  geom_col(fill = cols.nice[1])+
  theme_light()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle= 60, hjust= 1))


unify_axis(p.network.sizes)
# compare jaccard overlap: ------------------------------------------------

vertex.sim=
  sapply(net.list, function(x){
  map(net.list, function(y){
    jaccard_index(x, y)

  })%>% unlist()
})

edge.sim =
  sapply(net.list, function(x){
    map(net.list, function(y){
      jaccard_index(x, y, "edge")

    })%>% unlist()
  })


v.map = Heatmap(vertex.sim,
                name = "jaccard",
                col = cols.nice[c(3,2)],
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.1f", vertex.sim[i, j]), x, y, gp = gpar(fontsize = 8))
                  },
                show_row_dend = F,
                cluster_rows = F ,
                cluster_columns = F,
                show_column_dend = F)


e.map = Heatmap(edge.sim, name = "jaccard",
                col = cols.nice[c(4,1)],
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.2f", edge.sim[i, j]), x, y, gp = gpar(fontsize = 8))
                },
                show_row_dend = F,
                cluster_rows = F ,
                cluster_columns = F,
                show_column_dend = F)
e.map

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/disease_net_comparison_heatmaps.pdf",
    width= 5,
    height = 5)
print(v.map)
print(e.map)
print(d.map)
dev.off()

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/disease_net_comparison_size.pdf",
    width= 4,
    height = 5)
unify_axis(p.network.sizes)+theme(axis.title = element_blank())
dev.off()



# calculate delta con between netowrks ------------------------------------

df= sapply(net.list, function(x){
  sapply(net.list, function(y){
    deltacon_igraph(x,y)
  })
})

nets= c("HFnet_subset_M", "HFnet_subset_HPO","HPOnet_subset_HF", "Morbinet_subset_HF")
df2= df[nets,nets]
d.map = Heatmap(df2,
                name = "DetlaCon",
                col = cols.nice[c(5,1)],
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.2f", df2[i, j]), x, y, gp = gpar(fontsize = 8))
                },
                show_row_dend = F,
                cluster_rows = F ,
                cluster_columns = F,
                show_column_dend = F)
d.map

#deltacon_igraph(net.list$HFnet_subset_HPO,net.list$HPOnet_subset_HF)


#compare deltacon with randomized hf net

.gX = hf_net


# sequence of probabilities


run_randomization= function(.gX, probseq= seq(0,1, 0.1)){


dist.m= sapply(probseq, function(prob){
  sapply(c(1:5), function(seed){ # we repeat randomization 5x for each probability
    set.seed(seed)
    g.rand2= rewire(.gX, with = each_edge(prob= prob))
    E(g.rand2)$weight= sample(E(.gX)$weight)

    deltacon_igraph(.gX,g.rand2)
  })
})

colnames(dist.m)= probseq#paste0("Pr." , probseq)
df= dist.m%>% as_tibble()%>%
  pivot_longer(everything(), names_to = "Prob", values_to = "value")%>%
  mutate(Prob= as.numeric(Prob))
return(df)
}


hf_hpo= run_randomization(net.list$HFnet_subset_HPO)

hf_morb= run_randomization(net.list$HFnet_subset_M)
# Get mean of data values from data frame
hf_hpo = hf_hpo %>% mutate(net = "hpo.subset")
hf_morb = hf_morb %>% mutate(net = "M.subset")

mean1 <- hf_hpo %>%
  group_by(Prob) %>%
  summarize(average = median(value)) %>%
  ungroup()

mean2 <- hf_morb %>%
  group_by(Prob) %>%
  summarize(average = median(value)) %>%
  ungroup()

df= rbind(hf_hpo, hf_morb)

# Create Boxplot with a line plot using mean values
p.deltacon= df %>%
  #mutate(Prob= as.numeric(Prob))%>%
  ggplot(mapping = aes(x = Prob   , y = value)) +
  geom_jitter(color= "darkgrey") +
  geom_hline(yintercept = 0.46, color= "darkred", lty= 2)+
  geom_hline(yintercept = 0.39, color= "darkblue", lty= 2)+
  geom_line(data = mean1, #hpo
            mapping = aes(x = Prob, y = average, group=1),color="darkblue")+
  geom_line(data = mean2,# morbi
            mapping = aes(x = Prob, y = average, group=1),color="darkred")+
  theme_bw()+
  labs(y= "DeltaCon distance",
       x= "Rewire probability")+
  annotate("text", x=0.05, y=0.47, label="Morbinet", angle=0, size=4, color="darkred")+
  annotate("text", x=0.05, y=0.40, label="HPOnet", angle=0, size=4, color="darkblue")+
  scale_y_continuous(breaks = c(0.39, 0.40, 0.46, 0.6, 0.8, 1))+
  scale_x_continuous(breaks = seq(0, 1, 0.1))
p.deltacon

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/disease_net_detlcon.random.pdf",
    width= 4,
    height = 6)
#unify_axis(p.deltacon)
p.deltacon
dev.off()

edges =
  readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds")

xdd =randomize_links_by_layer(igraph::as_data_frame(edges$disease$heidelberg))

# Show edges conserved between both:  -------------------------------------

g1= net.list$hf_or
g2= net.list$morbinet_subset
g1_edg= igraph::as_data_frame(g1, what= "edges") %>% select(from,to) %>% as_tibble
g2_edg = igraph::as_data_frame(g2, what= "edges")%>% select(from,to)%>% as_tibble

# create an edgeID that is the same in both networks (sort the 2 nodes and than paste)
g1_IDs= map2(g1_edg$from, g1_edg$to, function(x,y){
  paste(sort(c(x,y)),collapse="_")
})
g1_edg= g1_edg %>% mutate(edgeID = unlist(g1_IDs))

g2_IDs= map2(g2_edg$from, g2_edg$to, function(x,y){
  paste(sort(c(x,y)),collapse="_")
})
g2_edg= g2_edg %>% mutate(edgeID = unlist(g2_IDs))

# now calculate jaccard based on this edge id union and intersect
u= union(g1_edg$edgeID, g2_edg$edgeID)
i= intersect(g1_edg$edgeID, g2_edg$edgeID)

length(g2_edg$edgeID)
length(i)
# hf net edges that intersect:
g1_edg %>% filter(edgeID %in% i) %>%
  left_join(Phedic%>% rename(from= PheCode,
                              Phenotype_from= Phenotype)%>% select(-category)) %>%
  left_join(Phedic%>% rename(to= PheCode,
                              Phenotype_to= Phenotype)%>% select(-category)) %>% print(n=500)

# hf net special edges:
g1_edg %>% filter(!edgeID %in% i) %>%
  left_join(Phedic%>% rename(from= PheCode,
                              Phenotype_from= Phenotype)%>% select(-category)) %>%
  left_join(Phedic%>% rename(to= PheCode,
                              Phenotype_to= Phenotype)%>% select(-category)) %>% print(n=500)
