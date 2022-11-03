## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-05-11
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
##  use disease cluster mapping to explore patient profiles
## ---------------------------
library(tidyverse)
library(igraph)
library(ComplexHeatmap)

dd.net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")
#link.data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data2.rds")
pids= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/pids_endpoints.rds")
c.info= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/full_clinic_info.rds")%>%
  distinct(patient_id, sex)%>% filter(sex != "u")%>% mutate(sex= ifelse(sex=="m", "male", ifelse(sex=="f", "female", sex)))


pids.interesting= list("hfpef"= pids.list$hfpef,
                       "hfref"= pids.list$hfref,
                       "hfmref"= pids.list$hfmref,
                       "hfmref"= pids.list$hfmref,
                       "female"= c.info%>% filter(sex=="female")%>% pull(patient_id)%>% unique(),
                       "male"= c.info%>% filter(sex=="male")%>% pull(patient_id)%>% unique()
                       )

pids.interesting= c(pids.interesting, pids)

node_df=igraph::as_data_frame(dd.net, "vertices")
diseasecluster= split(node_df$name, node_df$group_louv)


# transform patient profile to node sets ----------------------------------

# function to transform patient data into a list of diseases for each patient
create_nodesets= function(pids.list, data){

  #Prepare all patients
  all.pids = unique(unlist(pids.list))

  #prefilter to reduce comp time
  data2 = data %>%
    filter(pid %in% all.pids) %>%
    group_by(pid) %>%
    distinct(pid, PheCode) %>%
    drop_na() %>%
    arrange(pid)

  node_sets= split( data2$PheCode, data2$pid)

  return(node_sets)
}

# get sample and features
phecodes= V(dd.net)$name
targetpids= pids.list$hf_all

data.red =
  data %>%
  filter(PheCode %in% phecodes &
           pid %in% targetpids)

node_sets_all= create_nodesets(pids.list = targetpids, data.red)




# Get node centrality Df --------------------------------------------------

#Page ranks
pr= page_rank(dd.net, weights = E(dd.net)$weight)

#pr$vector= scale(pr$vector)

pr.frame= enframe(diseasecluster, name= "clust", value= "PheCode")%>% unnest(PheCode)%>%
  left_join(enframe(pr$vector, name= "PheCode", value= "PR"),  by= "PheCode")%>%
  group_by(clust)%>%
  mutate(clust.sum= sum(PR))

pr.frame

#Degree centrality:
str.frame= strength(dd.net, weights = E(dd.net)$weight)

str.frame= enframe(diseasecluster, name= "clust", value= "PheCode")%>% unnest(PheCode)%>%
  left_join(enframe(str.frame, name= "PheCode", value= "strength"),  by= "PheCode")%>%
  group_by(clust)%>%
  mutate(clust.sum= sum(strength))

str.frame


#network = enframe(node_sets_all,name= "pid", value= "PheCode")%>% unnest(PheCode)
wsum= sapply(node_sets_all[names(node_sets_all) %in% pids.list$hf_all], function(x){
  print(x)
  map(unique(str.frame$clust),function(y){
    print(y)
    pr.frame.filt= str.frame %>% filter(clust== y)
    prsum= pr.frame.filt %>% filter(PheCode %in% unlist(x))%>% pull(strength )%>% sum()
    if(prsum== 0){
      return(wsum= 0)}
    else{
      #total.sum= unique(pr.frame.filt$clust.sum)
      total2= str.frame%>% filter(PheCode %in% unlist(x))%>% pull(strength )%>% sum()

      wsum= prsum/(total2)
    }
  })%>% unlist()
})%>% unlist()



wsum.res= map(pids.list, function(x1){
  wsum= sapply(node_sets_all[names(node_sets_all) %in% x1], function(x){
    map(unique(str.frame$clust),function(y){
      pr.frame.filt= str.frame %>% filter(clust== y)
      prsum= pr.frame.filt %>% filter(PheCode %in% unlist(x))%>% pull(strength )%>% sum()
      if(prsum== 0){
        return(wsum= 0)}
      else{
        #total.sum= unique(pr.frame.filt$clust.sum)
        total2= str.frame%>% filter(PheCode %in% unlist(x))%>% pull(strength )%>% sum()

        wsum= prsum/(total2)
      }
    })%>% unlist()
  })%>% unlist()
})

wsum.df= cbind(wsum.res$hfpef, wsum.res$hfref)

wsum.s= apply(wsum.df, 1, scale)

rownames(wsum.s)= colnames(wsum.df)
colnames(wsum.s) = paste0("DC.", 1:dim(wsum.s)[2])
wsum.s


p.box= wsum.s %>%as.data.frame() %>% rownames_to_column("pid")%>%
  as_tibble() %>%
  mutate(pid= ifelse(pid %in% pids.list$hfpef,"hfpef", "hfref"))%>%
  pivot_longer(-pid, names_to= "cluster", values_to= "scaled_jacc")

p.box%>%
  ggplot(., aes(x= cluster, y= scaled_jacc, fill = pid))+
  geom_boxplot()

hfpef= apply(wsum.s[colnames(wsum.res$hfpef),], 2, median)
hfref= apply(wsum.s[colnames(wsum.res$hfref),], 2, median)

saveRDS(wsum.s, "output/scaled_wsum_network_profiles.rds")


Heatmap(rbind(hfpef,hfref), name = "mean scaled jaccard", rect_gp = gpar(col= "black"))



# jaccard! ----------------------------------------------------------------

targetpids= pids.list

jaccs.res= map(targetpids, function(x){
  jaccs.pef= sapply(node_sets_all[names(node_sets_all) %in% x], function(x){
    map(diseasecluster,function(y){
      length(intersect(x,y))/length(union(x,y))
    })%>% unlist()
  })

})

jaccs.res$hf_all[1:10, 1:10]
dim(jaccs.res$hf_all)

#scale the jaccard indices per CLUSTER
jaccs_sc= apply(jaccs.res$hf_all, 1, scale)
jaccs_sc[1:10, 1:10]
boxplot((jaccs_sc))
rownames(jaccs_sc)= colnames(jaccs.res$hf_all)

#tidy df

jacc.df = t(jaccs_sc) %>%
  as_tibble() %>%
  mutate(cluster= factor(names(diseasecluster),levels = names(diseasecluster)))%>%
  pivot_longer(-cluster)
# pids.interesting$all_patients = pids.list$hf_all
# pids.interesting= pids.interesting[names(pids.interesting) != "hfmref"]
# pids.interesting$hfmref = pids.list$hfmref
for (i in names(pids.interesting)){
  print(i)
   jacc.df= jacc.df %>%
     mutate({{i}} := ifelse(name %in% pids.interesting[[i]], i, "n"
                         )
            )
}

list.jaccs= lapply(names(pids.interesting), function(x){
  print(x)
  x= c(x, "cluster")
  df= jacc.df %>% group_by(across(all_of(x)))%>% summarise(mean= mean(value),
                                                       median= median(value))
  colnames(df)[1]= "group"
  df

})

plot.df = do.call(rbind, list.jaccs)

plot.df%>%
  ggplot(., aes(x= cluster, y=mean))+
  facet_grid(vars(rows= group))+
  geom_col(position = "dodge")+
  #scale_fill_manual(values= col.set)+
  theme_minimal()

p1 = plot.df%>%
  ggplot(., aes(x= cluster,y= group , fill =mean))+
  geom_tile()+
  #scale_fill_manual(values= col.set)+
  theme_minimal()

print(p1)
plot.hmap= plot.df %>% filter(group!= "n") %>% pivot_wider(-median, names_from = cluster, values_from = mean)
plot_h= column_to_rownames(plot.hmap, "group")
Heatmap(plot_h[c("hfpef", "hfref"),], cluster_columns = T,rect_gp = gpar(col= "black"))

Heatmap(plot_h[,], cluster_columns = T,rect_gp = gpar(col= "black"))

library(ComplexHeatmap)


p.mean_jacc= mean_jaccs %>%
  as_tibble() %>%
  mutate(cluster= factor(names(diseasecluster),levels = names(diseasecluster)))%>%
  pivot_longer(-cluster) %>%
  ggplot(., aes(x= cluster, y=value , fill = name))+
  facet_grid(vars(rows= name))+
  geom_col(position = "dodge")+
  scale_fill_manual(values= col.set)+
  theme_minimal()

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/jacc_patient_net_cluster.pdf",
    height = 5,
    width=4)
p.mean_jacc
dev.off()

# scale patient-wise jaccards

jacs.hf= cbind(jaccs.res$hfpef, jaccs.res$hfref, jaccs.res$hfmref)
jacs.hf[1:10, 1:10]

jacs.s= apply(jacs.hf, 1, scale)

rownames(jacs.s)= colnames(jacs.hf)
colnames(jacs.s) = paste0("DC.", colnames(jacs.s))


library(ComplexHeatmap)

ha = HeatmapAnnotation(bar = sample(letters[1:3], 10, replace = TRUE))
cohort= ifelse(rownames(jacs.s) %in% colnames(jaccs.res$hfpef), "hfpef", "hfref")

ra = rowAnnotation(cohort =cohort, col = list("cohort"= c("hfpef"= "black", "hfref"= "yellow")))

h.scaled_jac= Heatmap(jacs.s, show_row_names = F,#top_annotation = ha,
                      right_annotation = ra, cluster_columns = F, name = "scaled_jaccard")

p.box= jacs.s %>%as.data.frame() %>% rownames_to_column("pid")%>%
  as_tibble() %>%
  mutate(pid= ifelse(pid %in% pids.list$hfpef,"hfpef",
                     ifelse(pid %in% pids.list$hfref,"hfref","hfmref")
                     )
         )%>%
  pivot_longer(-pid, names_to= "cluster", values_to= "scaled_jacc")
unique(p.box$pid)
p.box%>%
  ggplot(., aes(x= cluster, y= scaled_jacc, fill = pid))+
  geom_boxplot()

hfpef= colMeans(jacs.s[colnames(jaccs.res$hfpef),])
hfref= colMeans(jacs.s[colnames(jaccs.res$hfref),])

hfpef= apply(jacs.s[colnames(jaccs.res$hfpef),], 2, median)
hfref= apply(jacs.s[colnames(jaccs.res$hfref),], 2, median)


hfpef= apply(jacs.s[colnames(wsum.res$hfpef),], 2, median)
hfref= apply(jacs.s[colnames(wsum.res$hfref),], 2, median)

saveRDS(jacs.s, "output/scaled_jaccard_network_profiles.rds")


Heatmap(rbind(hfpef,hfref), name = "mean scaled jaccard", rect_gp = gpar(col= "black"))

# run wilcox to test difference for each cluster between hfpef-hfref
x= 1
unique(jacc.df$cluster)
p.vals= map(unique(jacc.df$cluster), function(x){
  print(x)
  x1= jacc.df%>% filter(cluster== x, hfref=="hfref")%>% pull(value)
  x2= jacc.df%>% filter(cluster== x, hfpef=="hfpef")%>% pull(value)
  wilcox.test(x1, x2)$p.value

})
names(p.vals)= names(diseasecluster)

enframe(p.vals)%>% mutate(value = unlist(value))%>% mutate(logp= -log10(value), sig = (value<0.1))



p.vals= map(unique(jacc.df$cluster), function(x){
  wilcox.test(jaccs.res$hfpef[x, ],
              jaccs.res$hfref[x, ])$p.value

})
names(p.vals)= names(diseasecluster)

mat2= enframe(p.vals)%>%
  mutate(value = unlist(value))%>%
  mutate(logp= -log10(value),
         p.adj= p.adjust(value, method= "BH"),
         sig = (p.adj<0.05))



##run anovas when adding hfmref

jacc.df2= jacc.df %>% mutate(hf= ifelse(hfpef=="hfpef", "hfpef",
                              ifelse(hfref=="hfref", "hfref",
                                     ifelse(hfmref== "hfmref" , "hfmref", "n"
                                            )
                                     )
                              )

                   )%>%
  filter(hf !="n")




p.vals= map(seq(1:10), function(x){
  x1= jacc.df%>% filter(cluster== x, hfref=="hfref")%>% pull(value)
  x2= jacc.df%>% filter(cluster== x, hfpef=="hfpef")%>% pull(value)
  x3= jacc.df%>% filter(cluster== x, hfpef=="hfmref")%>% pull(value)
  anova()
  wilcox.test(x1, x2)$p.value

})

p.vals= map(seq(1:11), function(x){
  df= jacc.df2 %>% filter(cluster==x)
  p=kruskal.test(df$value~ df$hf)$p.value
  #-log10(p)
  #pairwise.wilcox.test(df$value, df$hf)$p.value
})

p.val.df= enframe(unlist(p.vals))%>%
  mutate(sig= ifelse(value<0.001, "<0.001",
                     ifelse(value < 0.01, "**",
                            ifelse(value< 0.05, "<0.05", "ns"))))


plot.hmap= plot.df %>% filter(group!= "n") %>% pivot_wider(-median, names_from = cluster, values_from = mean)

plot_h= column_to_rownames(plot.hmap, "group")
colnames(plot_h) = paste0("DC.", colnames(plot_h))

ha = HeatmapAnnotation(significance =p.val.df$sig,
                       col= list(significance= c("<0.001"= "darkred", "<0.05"= "red", "ns"= "grey")))

p_jaccard= Heatmap(plot_h[c("hfpef", "hfref"),],
        cluster_columns = F,
        cluster_rows = F,
        top_annotation = ha,
        rect_gp = gpar(col= "black"),
        name= "mean jaccard")

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/jacc_patient_net_cluster_sig.pdf",
    height = 2,
    width=5)
p_jaccard
dev.off()
