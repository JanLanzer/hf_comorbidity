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

source("analysis/utils/utils.R")

dd.net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet_clustered.rds")

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")
#link.data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data2.rds")
pids= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/pids_endpoints.rds")

c.info= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/full_clinic_info.rds")%>%
  distinct(patient_id, sex)%>% filter(sex != "u")%>% mutate(sex= ifelse(sex=="m", "male", ifelse(sex=="f", "female", sex)))


pids.interesting= list("hfpef"= pids.list$hfpef,
                       "hfref"= pids.list$hfref,
                       "hfmref"= pids.list$hfmref,
                       "female"= c.info%>% filter(sex=="female")%>% pull(patient_id)%>% unique(),
                       "male"= c.info%>% filter(sex=="male")%>% pull(patient_id)%>% unique()
                       )

pids.interesting= c(pids.interesting, pids[names(pids)!= "htx"])
names(pids.interesting)= c("HFpEF", "HFrEF", "HFmrEF", "Female", "Male", "Intubated", "Device.Implant", "PCI")

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



# calculate jaccard index of each patient with each cluster
jaccs.res= map(pids.list, function(x){
  jaccs.pef= sapply(node_sets_all[names(node_sets_all) %in% x], function(x){
    map(diseasecluster,function(y){
      length(intersect(x,y))/length(union(x,y))
    })%>% unlist()
  })

})

jaccs.res$hf_all[1:9, 1:9]
dim(jaccs.res$hf_all)

#scale the jaccard indices per CLUSTER
jaccs_sc= apply(jaccs.res$hf_all, 1, scale)
jaccs_sc[1:9, 1:9]

rownames(jaccs_sc)= colnames(jaccs.res$hf_all)

#tidy df

jacc.df = t(jaccs_sc) %>%
  as_tibble() %>%
  mutate(cluster= factor(names(diseasecluster),levels = names(diseasecluster)))%>%
  pivot_longer(-cluster)

for (i in names(pids.interesting)){
  print(i)
   jacc.df= jacc.df %>%
     mutate({{i}} := ifelse(name %in% pids.interesting[[i]], i, "n"
                         )
            )
}
library(ggpubr)
library(rstatix)
x= "HFpEF"
list.jaccs= lapply(names(pids.interesting), function(x){
  print(x)
  x= c(x, "cluster")
  df= jacc.df %>% group_by(across(all_of(x)))%>%
    summarise(mean= mean(value),
              median= median(value))

  colnames(df)[1]= "group"
  df

})

#calc p-values one vs all
list.jaccs.pvals= lapply(names(pids.interesting), function(x){
  print(x)

  df= jacc.df %>%
    group_by(cluster)%>%
    pairwise_wilcox_test(as.formula(paste(" value ~ ", x)), p.adjust.method = "BH")


})


plot.df = do.call(rbind, list.jaccs)%>%
  filter(!group %in% c("n","htx"))%>%
  mutate(cluster = paste0("DC.", cluster))

p1 = plot.df%>%
  mutate(group= factor(group, levels= rev(c("HFpEF", "HFrEF", "HFmrEF", "Female", "Male", "Intubated", "Device.Implant", "PCI"))))%>%
  ggplot(., aes(x= cluster,y= group , fill =mean))+
  geom_tile(col= "darkgrey")+
  scale_fill_gradient(low= "white", high= "darkred")+
  theme_minimal()+
  labs(fill= "scaled\nJaccard\nIndex")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle= 90, hjust= 1))+
  coord_equal()



pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/jacc_patient_net_cluster.pdf",
    height = 5,
    width=5)
print(unify_axis(p1))
dev.off()



plot.hmap= plot.df %>% filter(group!= "n") %>% pivot_wider(-median, names_from = cluster, values_from = mean)
plot_h= column_to_rownames(plot.hmap, "group")
Heatmap(plot_h[c("hfpef", "hfref"),], cluster_columns = T,rect_gp = gpar(col= "black"))

Heatmap(plot_h[,], cluster_columns = T,rect_gp = gpar(col= "black"))

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

mat2= enframe(p.vals)%>%
  mutate(value = unlist(value))%>%
  mutate( p.adj= p.adjust(value, method= "BH"),
          logp= -log10(p.adj),
         sig = (p.adj<0.01))



# ##run anovas when adding hfmref
#
# jacc.df2= jacc.df %>% mutate(hf= ifelse(hfpef=="hfpef", "hfpef",
#                               ifelse(hfref=="hfref", "hfref",
#                                      ifelse(hfmref== "hfmref" , "hfmref", "n"
#                                             )
#                                      )
#                               )
#
#                    )%>%
#   filter(hf !="n")
#
#
#
#
# p.vals= map(seq(1:10), function(x){
#   x1= jacc.df%>% filter(cluster== x, hfref=="hfref")%>% pull(value)
#   x2= jacc.df%>% filter(cluster== x, hfpef=="hfpef")%>% pull(value)
#   x3= jacc.df%>% filter(cluster== x, hfpef=="hfmref")%>% pull(value)
#   anova()
#   wilcox.test(x1, x2)$p.value
#
# })
#
# p.vals= map(seq(1:11), function(x){
#   df= jacc.df2 %>% filter(cluster==x)
#   p=kruskal.test(df$value~ df$hf)$p.value
#   #-log10(p)
#   #pairwise.wilcox.test(df$value, df$hf)$p.value
# })

p.val.df= enframe(unlist(p.vals))%>%
  mutate(sig= ifelse(value<0.0001, "*", "ns"))
p.val.df= p.val.df %>%
  mutate(sig= ifelse(value<0.0001, "<1e-04",
                     ifelse(value<0.01,"<1e-02",
                            "ns")))

plot.hmap= plot.df %>% filter(group!= "n") %>% pivot_wider(-median, names_from = cluster, values_from = mean)

plot_h= column_to_rownames(plot.hmap, "group")


ha = HeatmapAnnotation(HFpEF_v_HFrEF =p.val.df$sig,
                       col= list(HFpEF_v_HFrEF= c("<1e-04"= cols.nice[1],
                                                  "<1e-02"=  cols.nice[2],
                                                  "ns"= "grey")))
library(circlize)
col_fun = colorRamp2(c(0,0.8,1.5), c("white","red", "darkred"))

p_jaccard= Heatmap(as.matrix(plot_h),
                   col= col_fun,
        cluster_columns = F,
        cluster_rows = F,
        top_annotation = ha,
        rect_gp = gpar(col= "darkgrey"),
        name= "sacled\nJaccard\nIndex")

p_jaccard

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/jacc_patient_net_cluster_sig.pdf",
    height = 3.5,
    width=6)
p_jaccard
dev.off()


# add disease signatures --------------------------------------------------

comp_feat= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")


hfref= comp_feat%>% arrange(desc(estimate))%>% pull(PheCode)
hfpef= comp_feat%>% arrange((estimate))%>% pull(PheCode)
dis.sig= list("hfpef"= hfpef[1:50],
              "hfref"= hfref[1:50])
sap
diseasecluster

df= sapply(diseasecluster, function(x){
  sapply(dis.sig, function(y){
    length(intersect(x,y))/length(union(x,y))
  })
})

rownames(df)= c("HFpEF", "HFrEF")
colnames(df)= paste0("DC.", colnames(df))
p.jacc.dis.sig= df %>% as.data.frame()%>%
  rownames_to_column("comparison")%>%
  pivot_longer(-comparison)%>%
  mutate(comparison= factor(comparison, levels=c( "HFrEF","HFpEF") ))%>%
  ggplot(., aes(x= name, y = comparison, fill = value))+
  geom_tile(color= "black")+
  scale_fill_gradient(low= "white", high= "darkred")+
  theme_minimal()+
  #geom_text(aes(label= sig))+
  theme(axis.text= element_text(colour = "black"),
        axis.text.x = element_text(angle= 45, hjust= 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  coord_equal()+
  labs(fill ="Jaccard \n Index",x="",  y= "")


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/jacc_patient_sig.pdf",
    height = 3.5,
    width=6)
unify_axis(p.jacc.dis.sig)
dev.off()

lapply(dis.sig, function(x){
  GSE_analysis(x, diseasecluster)

})


GSE_analysis = function(geneList,Annotation_DB){
  library(dplyr)
  library(tidyr)
  library(tibble)

  geneList = geneList[geneList %in% unique(unlist(Annotation_DB))]

  ResultsDF = matrix(0,nrow = length(Annotation_DB),ncol = 5)
  rownames(ResultsDF) = names(Annotation_DB)
  colnames(ResultsDF) = c("GenesInPathway","GenesInList","GeneNames","p_value","corr_p_value")

  DB_genecontent = length(unique(unlist(Annotation_DB)))

  GenesDB = DB_genecontent
  SelectedGenes = length(geneList)

  for(gset in rownames(ResultsDF)){
    GP = length(((Annotation_DB[[gset]])))
    GL = length(intersect(Annotation_DB[[gset]],geneList))

    ResultsDF[gset,"GenesInList"] = GL
    ResultsDF[gset,"GenesInPathway"] = GP
    ResultsDF[gset,"GeneNames"] = paste(intersect(Annotation_DB[[gset]],geneList),collapse = ",")
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)
  }

  ResultsDF[,"corr_p_value"] = p.adjust(ResultsDF[,"p_value"],method = "BH")
  ResultsDF = data.frame(ResultsDF,stringsAsFactors = F)
  ResultsDF = ResultsDF[order(ResultsDF[,"p_value"]),]

  ResultsDF = ResultsDF %>%
    rownames_to_column("gset") %>%
    mutate_at(c("GenesInPathway","GenesInList",
                "p_value","corr_p_value"),
              as.numeric) %>%
    dplyr::arrange(corr_p_value,GenesInList)

  return(ResultsDF)

}

