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
library(ggpubr)
library(rstatix)
library(circlize)
library(lubridate)


source("analysis/utils/utils.R")

dd.net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet_clustered.rds")

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")
pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")
#link.data= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data2.rds")
pids= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/pids_endpoints.rds")

pid.df= read.csv( "T:/fsa04/MED2-HF-Comorbidities/data/RWH_September2022/raw/levinson_comorbidities_pids_2022-10-06.csv", sep = ";") %>%
  as_tibble


#add age calculation
pid.df= pid.df  %>%
  filter(pid %in% pids.list$hf_all)%>%
  mutate(birthday= lubridate::as_date(birthday),
         entry_date = lubridate::as_date(diag_date_min ),
         ICDint = lubridate::interval(birthday, entry_date),
         age.at.icd = ICDint/dyears(1)) %>%
  select(pid, sex, age.at.icd)

pids.interesting= list("HFpEF"= pids.list$hfpef,
                       "HFrEF"= pids.list$hfref,
                       "HFmrEF"= pids.list$hfmref,
                       "Female"=pid.df %>% filter(sex=="f")%>% pull(pid)%>% unique(),
                       "Male"= pid.df%>% filter(sex=="m")%>% pull(pid)%>% unique(),
                       "40_60"= pid.df%>% filter(age.at.icd<60)%>% pull(pid)%>% unique(),
                       "60_80"= pid.df%>% filter(age.at.icd>=60 & age.at.icd<=80)%>% pull(pid)%>% unique(),
                       "80+"= pid.df%>% filter(age.at.icd>80)%>% pull(pid)%>% unique()
                       )

map(pids.interesting, length)
pids.interesting= c(pids.interesting, pids[names(pids)!= "htx"])
names(pids.interesting)= c("HFpEF", "HFrEF", "HFmrEF", "Female", "Male", "40-60y","60-80y",
                           "+80y", "Intubated",
                           "Device.Implant", "PCI")

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

t(jaccs.res$hf_all)%>% as.data.frame()%>% rownames_to_column("pid")%>%
  pivot_longer(names_to= "DC", values_to= "jac", -pid)
boxplot(t(jaccs.res$hf_all))

boxplot(jaccs_sc)

#tidy df

#jacc.df = t(jaccs_sc) %>%
jacc.df = jaccs.res$hf_all %>%
 as_tibble() %>%
  mutate(cluster= factor(names(diseasecluster),levels = names(diseasecluster)))%>%
  pivot_longer(-cluster)

jacc.df.s = t(jaccs_sc) %>%
  as_tibble() %>%
  mutate(cluster= factor(names(diseasecluster),levels = names(diseasecluster)))%>%
  pivot_longer(-cluster)%>%
  dplyr::rename(value_scaled= value)

jacc.df= jacc.df %>%
  left_join(jacc.df.s)

for (i in names(pids.interesting)){
  print(i)
   jacc.df= jacc.df %>%
     mutate({{i}} := ifelse(name %in% pids.interesting[[i]], i, "n"
                         )
            )
}

## plot boxplots
jacc.df %>%
  mutate(age = ifelse(name %in% pids.interesting$`40-60y`, "40-60y", ""),
         age = ifelse(name %in% pids.interesting$`60-80y`, "60-80y", age),
         age = ifelse(name %in% pids.interesting$`+80y`, "80+y",age))%>%
  ggplot(aes(x= cluster, y= value_scaled , fill = age ))+
  geom_violin()+
  geom_boxplot()

jacc.df %>%
  mutate(HF = ifelse(name %in% pids.interesting$HFpEF, "HFpEF", ""),
         HF = ifelse(name %in% pids.interesting$HFrEF, "HFrEF", HF),
         HF = ifelse(name %in% pids.interesting$HFmrEF, "HFmrEF",HF))%>%
  filter(HF != "")%>%
  ggplot(aes(x= cluster, y= value_scaled , fill = HF ))+
  geom_violin()+
  geom_boxplot()


list.jaccs= lapply(names(pids.interesting), function(x){
  print(x)
  x= c(x, "cluster")
  df= jacc.df %>% group_by(across(all_of(x)))%>%
    summarise(mean= mean(value),
              median= median(value),
              mean.s= mean(value_scaled),
              median.s= median(value_scaled))

  colnames(df)[1]= "group"
  df

})

# #calc p-values one vs all
# list.jaccs.pvals= lapply(names(pids.interesting), function(x){
#   print(x)
#
#   df= jacc.df %>%
#     group_by(cluster)%>%
#     pairwise_wilcox_test(as.formula(paste(" value ~ ", x)), p.adjust.method = "BH")
#
#
# })

plot.df = do.call(rbind, list.jaccs)%>%
  filter(!group %in% c("n","htx"))%>%
  mutate(cluster = paste0("DC.", cluster))

p1 = plot.df%>%
  #mutate(group= factor(group, levels= rev(c("HFpEF", "HFrEF", "HFmrEF", "Female", "Male", "Intubated", "Device.Implant", "PCI"))))%>%
  ggplot(., aes(x= cluster,y= group , fill =mean.s))+
  geom_tile(col= "darkgrey")+
  scale_fill_gradient(low= "white", high= "darkred")+
  theme_minimal()+
  labs(fill= "scaled\nJaccard\nIndex")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle= 90, hjust= 1))+
  coord_equal()


p1

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/jacc_patient_net_cluster.pdf",
    height = 5,
    width=5)
print(unify_axis(p1))
dev.off()



plot.hmap= plot.df %>%
  filter(group!= "n") %>%
  select(-median, -median.s) %>%
  pivot_wider(-mean,  names_from = cluster, values_from = mean.s)


# plot.hmap= plot.df %>%
#   filter(group!= "n") %>%
#   select(-median, -median.s) %>%
#   pivot_wider(-mean,  names_from = cluster, values_from = median.s)

plot_h= column_to_rownames(plot.hmap, "group")

top.c= igraph::as_data_frame(dd.net, "vertices")%>%
  group_by(group_louv)%>%mutate(prev= 10^size)%>%
  summarise(x= median(prev))

top.c= jacc.df%>% group_by(cluster)%>%
  summarise(m.c= mean(value))

top.c

column_ha = HeatmapAnnotation("prevalence" = anno_barplot(top.c$x))

h1= Heatmap(plot_h[c("HFpEF", "HFrEF", "HFmrEF", "Male", "Female", "40-60y","60-80y", "+80y"),],
            cluster_columns = F,
            cluster_rows= F,
            rect_gp = gpar(col= "black"),
            split = c(rep("HF", 3), rep("Sex", 2), rep("Age", 3)),
            #rect_gp = gpar(col= "darkgrey"),
            name= "scaled\nJaccard\nIndex",
            border= T,
            gap = unit(0.5, "cm"),
            top_annotation = column_ha    )

h1

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/jacc_patient_net_cluster2.pdf",
    height = 3.5,
    width=6)
h1
dev.off()

range(plot_h)
col_fun = colorRamp2(c(-0.2,0, 1.31), c("darkblue","white", "darkred"), space= "RGB", transparency = 0)
max(plot_h)
p_jaccard= Heatmap(as.matrix(plot_h),
                   col= col_fun,
                   cluster_columns = F,
                   cluster_rows = F,
                   #top_annotation = ha,
                   row_names_side = "left",

                   rect_gp = gpar(col= "darkgrey"),
                   name= "sacled\nJaccard\nIndex")
p_jaccard
Heatmap(plot_h, cluster_columns = T,rect_gp = gpar(col= "black"))

# run wilcox to test difference for each cluster between hfpef-hfref
x= 1
test.col1 = "80+"
test.col2 = "60_80"
unique(jacc.df$cluster)

run_wilcox= function(jacc.df, test.col1, test.col2, alpha.level= 0.01){

  p.vals= map(unique(jacc.df$cluster), function(x){
    #print(x)
    unique(jacc.df[[test.col1]])
    x1= jacc.df%>% filter(cluster== x, (!!sym(test.col1))==test.col1)%>% pull(value)
    x2= jacc.df%>% filter(cluster== x, (!!sym(test.col2))==test.col2)%>% pull(value)
    c("p.val"= wilcox.test(x1, x2)$p.value,wilcox.test(x1, x2, conf.int= T)$estimate)

  })
  p.vals[[1]]
  enframe(p.vals) %>%
    unnest_longer(value)%>%
    pivot_wider(names_from = value_id, values_from = value)%>%
    mutate(logp= -log10(p.val),
           contrast=paste0(test.col1, " vs. ", test.col2))


  # p.vals= enframe(p.vals)%>%
  #   mutate(value = unlist(value),
  #          p.value = p.adjust(value, "BH"))%>%
  #   mutate(logp= -log10(p.value), sig = (p.value<alpha.level),
  #          contrast=paste0(test.col1, " vs. ", test.col2))
  # # p.vals %>%
  # #   mutate(sig= ifelse(value<0.0001, "<1e-04",
  # #                      ifelse(value<0.01,"<1e-02",
  # #                             "ns")))
}

df1= run_wilcox(jacc.df,
            test.col1 = "Female",
           test.col2 = "Male")

df2= run_wilcox(jacc.df,
               "+80y",
               "60-80y")

df7= run_wilcox(jacc.df,
                "+80y",
                "40-60y")

df3= run_wilcox(jacc.df,
                "60-80y",
               "40-60y")

df4= run_wilcox(jacc.df,
               "HFpEF",
               "HFmrEF")

df5= run_wilcox(jacc.df,
                test.col1= "HFpEF",
                test.col2= "HFrEF")

df6= run_wilcox(jacc.df,
                "HFmrEF",
                "HFrEF")





p.vals.compare= rbind(df1,df2,df3,df4,df5,df6)%>%
  mutate(p.adj= adjust_pvalue(p.val, method= "BH"),
         logp= -log10(p.adj),
         sig = (p.adj<0.01),
         s.effect= logp*sign(`difference in location` ),
         name= paste0("DC.", name),
         name= as.factor(name),
         )

p.p.vals= p.vals.compare%>%
  ggplot(aes(x= name, y= s.effect, fill= name))+
  geom_col()+
  scale_fill_manual(values= c(cols.nice, "darkred", "darkblue", "darkgreen", "black" ))+
  geom_hline(yintercept = c(-log10(0.01)), color= "darkgrey")+
  geom_hline(yintercept = c(log10(0.01)), color= "darkgrey")+
  #geom_hline(yintercept = c(-log10(0.00001)))+
  facet_grid(rows = vars(contrast))+
  theme_bw()+
  labs(x= "",
       y= "log10(adj.P) * sign of estimate",
       fill = "")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "left")

p.p.vals

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/jacc_patient_net_cluster_pvals.pdf",
    height = 12,
    width=4)
unify_axis(p.p.vals)
dev.off()

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

t(as.matrix(scale(t(plot_h), center = F)))

ha = HeatmapAnnotation(HFpEF_v_HFrEF =p.val.df$sig,
                       col= list(HFpEF_v_HFrEF= c("<1e-04"= cols.nice[1],
                                                  "<1e-02"=  cols.nice[2],
                                                  "ns"= "grey")))

col_fun = colorRamp2(c(0,1.5), c("white","yellow", "darkred"))

p_jaccard= Heatmap(as.matrix(plot_h),
                  col= col_fun,
        cluster_columns = F,
        cluster_rows = F,
        #top_annotation = ha,
        row_names_side = "left",

        rect_gp = gpar(col= "darkgrey"),
        name= "sacled\nJaccard\nIndex")

p_jaccard

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/jacc_patient_net_cluster_sig.pdf",
    height = 3.5,
    width=6)
p_jaccard
dev.off()


# ORA? --------------------------------------------------------------------
clust= 1
col= "HFpEF"


dfx= jacc.df %>% filter(cluster==clust)
vec1= dfx$value
names(vec1)= dfx$name

vec1= sort(vec1, decreasing = T)
list1= list("HFpEF"= dfx%>% filter(get(col)!= "n")%>% pull(name),
            "HFrEF"= dfx%>% filter(HFrEF!= "n")%>% pull(name)
)
fgseaSimple(pathways =list1 ,
            stats = vec1, nperm = 100, scoreType = "pos")

plotEnrichment(list1$HFrEF, stats = vec1)
# add disease signatures --------------------------------------------------

names(diseasecluster)= paste0("DC.", names(diseasecluster))
comp_feat= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFpEF_classifier_features.rds")


hfref= comp_feat%>% arrange(desc(estimate))%>% pull(PheCode)
hfpef= comp_feat%>% arrange((estimate))%>% pull(PheCode)

dis.sig= list("HFpEF"= hfpef[1:71],
              "HFrEF"= hfref[1:29])


##jaccard
df= sapply(dis.sig, function(x){
  sapply(diseasecluster, function(y){
    length(intersect(x,y))/length(x)
  })
})


col_fun = colorRamp2(c(0,0.32), c("white", "darkred"))

p.profile_percentage= Heatmap(t(df), col = col_fun, row_names_side = "left",
        cluster_rows = F,
        cluster_columns= F,
        name= "% of\ncomorbidity\nprofile",
        border= T,
        rect_gp = gpar(col = "black", lwd = 1)
        )

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/patient_cluster_profile.pdf",
    height = 1.5,
    width=4)
p.profile_percentage
dev.off()


##jaccard
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




# nomclust ----------------------------------------------------------------

s= nomclust(data= mca.df, measure= "goodall1", eval = T)

