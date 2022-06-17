## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-12-06
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## perform MCA analyis and associate dimensions with clinical co-variates
## ---------------------------

library(FactoMineR)
library(tidyverse)
library(GDAtools)
library(tidyverse)
library(qdapTools)
library(magrittr)
library(factoextra)
library(ggrepel)



# Read data and prepare ---------------------------------------------------

# read tables.
directory= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/"

icd = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe.rds")
#pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids.rds")
map(pids.list,length)
phecodes= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/top300_disease.rds")

pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/pidslist_oct2021.rds")

hf = read.csv("T:fsa04/MED2-HF-Comorbidities/data/RWH_March2020/levinson_comorbidities_in_hf_patients_2020-03-25.csv",
              sep = ";",
              na.strings=c("","NA")) %>% as_tibble

pheno.data= readRDS(file= "T:/fsa04/MED2-HF-Comorbidities/data/processed_data/full_clinic_info.rds")
source("~/GitHub/RWH_analysis/scripts/utils.R")


# Functions ---------------------------------------------------------------

## function to prepare df as input for MCA()
#' @param df,  df in tidy format with 2 columns (pid, PheCode)

MCAinput = function(df){
  #transform to count matrix with disease as columns, patients as rows and 0 and 1 defining the disease profile
  df= split(as.character(df$PheCode), df$pid) %>%
    lapply(function(x) strsplit(x, "-")) %>%
    mtabulate()%>%
    matrix2df("pid") %>%
    as_data_frame(check.names = FALSE) %>%
    column_to_rownames("pid")

  df2= data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
  rownames(df2) = rownames(df)
  return(df2)
}

## funciton to add gender to the mca.df, the output from MCAinput()
#' @param mca.df , ouput from MCAinput()

MCAinput_add_gender = function(mca.df){
  mca.df.gender = mca.df %>%
    rownames_to_column("pid") %>%
    left_join(hf %>%
                mutate(pid= as.character(pid)) %>%
                select(pid,sex), by= "pid") %>%
    filter(sex %in% c("m", "f")) %>%
    mutate(sex= factor(sex, levels=c("m", "f"))) %>%
    column_to_rownames("pid")
  return(mca.df.gender)
}

## function to create a ggplot friendly data frame from the output of MCA()
#' @param mca.res, takes the output from MCA function

plot_MCA_df = function(mca.res, ndim=3){
  var= as.data.frame(mca.res$var$coord[,1:ndim]) %>%
    rownames_to_column("ID") %>%
    mutate(VAR= "var")
  ind = as.data.frame(mca.res$ind$coord[,1:ndim])%>%
    rownames_to_column("ID") %>%
    mutate(VAR="ind")
  qual.sup= as.data.frame(mca.res$quali.sup$coord[,1:ndim])%>%
    rownames_to_column("ID") %>%
    mutate(VAR="quali.sup")

  df= do.call(rbind, list(var,ind,qual.sup))
  colnames(df) = str_replace_all(colnames(df), " ", "")

  return(as_tibble(df))
}

## wrapper for multiple plots from pacakge factoextra

plot_MCA_results = function(res.mca){
  x= fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 45))

  x2= fviz_mca_biplot(res.mca,
                      repel = F, # Avoid text overlapping (slow if many point)
                      ggtheme = theme_minimal(),
                      geom = "point",
                      label= c("var", "quali.sup"))

  x3= fviz_mca_ind(res.mca, col.ind = "cos2",
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = F, # Avoid text overlapping (slow if many points)
                   ggtheme = theme_minimal(),
                   geom= "point",
                   col.quali.sup = "black",
                   label= c("var", "quali.sup"))

  x4= fviz_contrib(res.mca, choice = "var", axes = 1, top = 15)
  return(list(x,x2,x3,x4))
}




# Data subsetting (pre MCA) -----------------------------------------------------------------------
icd_red = icd %>% drop_na %>% distinct(pid, PheCode) %>% filter(pid %in% c(pids.list$hfpef, pids.list$hfref),
                                                                PheCode %in% phecodes)

# icd_red is full dataset to work with. This will now be subsetted to answer specific questions of interest.
length(unique(icd_red$pid))
length(unique(icd_red$PheCode))

# MCA_1, all patients, all diagnosis -----------------------------------------------------------------------

#### use all patients, use all diagnostics and analyze impact of HF and Gender as qualitativ variables
## 1) transform full data set into input format for MCA function

mca.df=MCAinput(df = icd_red)

colnames(mca.df) = gsub("^X", "", colnames(mca.df))

mca.df[1:10, 1:10]

dim(mca.df)

ndims= length(unique(icd_red$PheCode))

mca.res = MCA(mca.df, ncp = ndims)

mca.res.df= plot_MCA_df(mca.res, ndim  =ndims )



# add clinical covariates -------------------------------------------------
#add.pheno.data

ind.df= mca.res.df %>%
  filter(VAR== "ind")
#%>%
#  left_join(pheno.data %>% dplyr::rename(ID= patient_id)%>% mutate(ID =as.character(ID)), by= "ID")

#add.patient.cohort:
ind.df= ind.df%>%   mutate(hf.type = ifelse(ID %in% pids.list$hfref,
                                       "hfref",
                                       ifelse(ID %in% pids.list$hfpef,
                                              "hfpef",
                                              ifelse(ID %in% pids.list$hfmref,
                                                     "hfmref",
                                                     "unlabeled")))
)

# add clinical endpoints:
endpoint= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/pids_endpoints.rds")
ind.df= ind.df%>%
  mutate(htx= ifelse(ID %in% endpoint$htx, "yes", "no"),
         intu= ifelse(ID %in% endpoint$intubation, "yes", "no"),
         defi= ifelse(ID %in% endpoint$defi, "yes", "no"),
         pci= ifelse(ID %in% endpoint$pci, "yes", "no"))

#add summary tableS=
sum.t= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/summary.table.pef.ref.rds")

sum.tfilt= sum.t %>%
  filter(pid %in% ind.df$ID) %>%
  distinct(pid,
           patient_cohort,
           sex,
           age.at.mean,
           nyha.max,
           median.BMI,
           mean.sys,
           mean.dias,
           median.LDL,
           median.HDL,
           median.Trigs,
           median.Chol,
           median.gfrcg,
           median.hba1c,
           median.bnp,
           PheCode_count,
           icd10gm_count,
           charlson_score,
           elixhauser_wscore,
           ef.min,
           edd.max,
           ea.min,
           ea.max,
           ee.max)


ind.df= ind.df %>% left_join(sum.tfilt%>% dplyr::rename(ID= pid)%>% mutate(ID= as.character(ID)), by= "ID")

## run wilcox to see if it associates dim1 with hf


dim_= paste0("Dim", seq(1:ndims-1))

df= ind.df
## loop over categorical variables , perform wilcox
cat.vars= c("hf.type", "sex", "htx", "intu", "defi", "pci")

tested.vars= map(cat.vars, function(x){
  print(x)
  levels= unique(df[[x]])
  print(levels)
  hf_t = df %>% filter((!!as.symbol(x))== levels[1])
  hf_f = df %>% filter((!!as.symbol(x))== levels[2])

  pvals= map(dim_, function(x){
    test. = wilcox.test(hf_t %>% pull(x),
                        hf_f%>% pull(x),
                        alternative = "two.sided")
    test.$p.value
  }) %>%unlist()

  #combine with variance explained:
  eig.val <- get_eigenvalue(mca.res)

  df.variance= eig.val[1:length(pvals),]

  associated_hf =cbind(df.variance, pvals) %>% as_tibble%>% mutate(sig= ifelse(pvals<0.05, "sig", "ns"))%>% group_by(sig)%>%
    summarise(s= sum(variance.percent)) %>% mutate(var= x)

})

p.cat= do.call(rbind, tested.vars)%>%filter(sig=="sig")%>%
  ggplot(., aes(x= reorder(var, -s), y= s))+
  geom_col()+
  labs(y= "percentage of variance explained (%)",
       x= "tested covariate")+
  theme(legend.position =  "none")#+
  #theme_minimal_hgrid()

## loop over continous variables, perform cor.test
cont.vars= c("ea.max", "median.bnp",
             "mean.sys",
             "mean.dias",
             "edd.max",
             "age.at.mean",
             "charlson_score",
             "median.BMI",
             "median.LDL",
             "median.Chol",
             "median.HDL",
           "median.Trigs",
           "median.gfrcg",
           "median.hba1c")

tested.vars.cont= map(cont.vars, function(x){
  print(x)

  pvals= map(dim_, function(y){
    test. = ind.df %>% select(!!as.symbol(x),
                              !!as.symbol(y))
    cor.test(test.[[x]], test.[[y]])$p.value
  }) %>%unlist()

  #combine with variance explained:
  eig.val <- get_eigenvalue(mca.res)

  df.variance= eig.val[1:length(pvals),]

  associated_hf =cbind(df.variance, pvals) %>%
    as_tibble%>%
    mutate(sig= ifelse(pvals<0.05, "sig", "ns"))%>%
    group_by(sig)%>%
    summarise(s= sum(variance.percent)) %>%
    mutate(var= x)

})

tested.vars.cont
#
# p.cont= do.call(rbind, tested.vars.cont)%>%
#   filter(sig=="sig")%>%
#   ggplot(., aes(x= reorder(var, -s), y= s))+
#   geom_col()+
#   labs(y= "percentage of variance explained (%)",
#        x= "tested covariate")+
#   theme(legend.position =  "none")+
#   theme_minimal_hgrid()
cowplot::plot_grid(p.cat, p.cont, align = "h")


#plot both tests together:

p.explained_V= do.call(rbind, c(tested.vars.cont, tested.vars))%>%
  mutate(var.type= ifelse(var %in% cat.vars, "cat","cont"))%>%
  filter(sig=="sig")%>%
  ggplot(., aes(x= reorder(var, s), y= s, fill = var.type))+
  geom_col()+
  labs(y= "estimated percentage of variance explained (%)",
       x= "tested covariate")+
  theme(legend.position =  "none")+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(color= "black"))

p.explained_V


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/explained_V.pdf",
    height = 3, width= 5)
p.explained_V

dev.off()


#DiM red---------------------------------------------------------

library(ggExtra)
source("~/GitHub/RWH_analysis/scripts/HF_classifier/classifier_functions.R")
set.seed(2)


#### UMAP
df= mca.df%>% mutate_all(as.numeric)
umap_res= umap(df)

umap.plot = as_tibble(cbind(as_tibble(umap_res$layout), pid=  rownames(df))) %>%
  mutate(hf= ifelse(pid %in% pids.list$hfpef,"hfpef", "hfref"),
         hf= factor(hf, levels= c("hfpef", "hfref") )) %>%
  left_join(sum.tfilt%>% mutate(pid= as.character(pid)), by= "pid")%>%
  mutate(HTX= ifelse(pid %in% endpoint$htx, "yes", "no"),
         Intubation= ifelse(pid %in% endpoint$intubation, "yes", "no"),
         ICD_implant= ifelse(pid %in% endpoint$defi, "yes", "no"),
         PCI= ifelse(pid %in% endpoint$pci, "yes", "no"))

p.umap.hf =ggplot(umap.plot, aes(x= V1, y= V2, color = hf))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_manual(values= col.set)+
  labs(color="HF cohort",
       x= "",
       y = "")+
  theme_classic()

p.umap.defi =ggplot(umap.plot, aes(x= V1, y= V2, color = ICD_implant))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_manual(values= col.set)+
  labs( x= "",
       y = "")+
  theme_classic()

p.umap.htx =ggplot(umap.plot, aes(x= V1, y= V2, color = HTX))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_manual(values= col.set)+
  labs( x= "",
        y = "")+
  theme_classic()

p.umap.sex =ggplot(umap.plot%>% mutate(sex=factor(sex, levels=c("m", "f"))), aes(x= V1, y= V2, color = sex))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_manual(values= col.set)+
  labs( x= "",
        y = "")+
  theme_classic()

p.umap.pci=ggplot(umap.plot, aes(x= V1, y= V2, color = PCI))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_manual(values= col.set)+
  labs( x= "",
        y = "")+
  theme_classic()

##

p.umap.age =ggplot(umap.plot, aes(x= V1, y= V2, color = age.at.mean))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_gradient(low= "orange", high= "black")+
  labs( x= "",
        y = "")+
  theme_classic()


p.umap.hba1c =ggplot(umap.plot, aes(x= V1, y= V2, color = median.hba1c))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_gradient(low= "orange", high= "black")+
  labs( x= "",
        y = "")+
   theme_classic()



p.umap.bmi =ggplot(umap.plot%>% filter(median.BMI<50), aes(x= V1, y= V2, color = median.BMI))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_gradient(low= "orange", high= "black")+
  labs( x= "",
        y = "")+
  theme_classic()


p.umap.bnp =ggplot(umap.plot, aes(x= V1, y= V2, color = log10(median.bnp)))+
  geom_point(alpha= 0.7, size= 1.)+
  scale_color_gradient(low= "orange", high= "black")+
  labs( x= "",
        y = "")+
  theme_classic()


p.umap.hba1c
main.umaps= cowplot::plot_grid( p.umap.defi, p.umap.htx, p.umap.hf,p.umap.sex, p.umap.age,p.umap.bmi, ncol = 2)

p.umap.age

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/manuscript/main/umaps_cat.pdf")
cowplot::plot_grid(p.umap.hf,p.umap.sex, p.umap.htx, p.umap.defi, p.umap.pci)
dev.off()

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/manuscript/main/umaps_cont.pdf")
    #heigh= 3, width= 6)
cowplot::plot_grid(p.umap.age, p.umap.hba1c, p.umap.bmi,p.umap.bnp)
dev.off()

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/umaps_main.pdf",
    height= 6, width= 6)
main.umaps
dev.off()

# umap with single disease feats ------------------------------------------
main_disease= c("585.3", "401.1", "272.13", "250.2", "296.22", "496",
              "280.1","327.3","411.4")

pls= map(main_disease, function(feature){
  pids.x= data %>% filter(PheCode == feature) %>% pull(pid)
  pheno.x= data%>%  filter(PheCode == feature) %>% pull(Phenotype)%>% unique()
  umap.plot %>% mutate(dis1= ifelse(pid %in% pids.x, "y", "n"))%>%
    ggplot(., aes(x= V1, y= V2, color = dis1))+
    geom_point(alpha= 0.7, size= .5)+
    scale_color_manual(values= col.set)+
    labs( x= "",
          y = "")+
    theme_classic()+
    ggtitle(paste0(feature, pheno.x))
})
pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/umaps_supp.pdf",
    height= 6, width= 6)
pls
dev.off()

print(pls[[1]])


## tSNE
## tSNE
df= mca.df[!duplicated(mca.df), ]
dim(df)

tsne_res= Rtsne(df, check_duplicates = F, perplexity = 50)

patient.df = tibble(pid= rownames(df)) %>%
  mutate(hf= ifelse(pid %in% pids.list$hfpef,"hfpef", "hfref"),
         hf= factor(hf, levels= c("hfpef", "hfref") )) %>%
  left_join(sum.tfilt%>% mutate(pid= as.character(pid)), by= "pid")%>%
  mutate(htx= ifelse(pid %in% endpoint$htx, "yes", "no"),
         intu= ifelse(pid %in% endpoint$intubation, "yes", "no"),
         defi= ifelse(pid %in% endpoint$defi, "yes", "no"),
         pci= ifelse(pid %in% endpoint$pci, "yes", "no"))


tsne.plot = as_tibble(cbind(as_tibble(tsne_res$Y), pid=  rownames(df))) %>%
  left_join(patient.df, by= "pid")

tsne= ggplot(tsne.plot, aes(x= V1, y= V2, color = hf))+
  geom_point()+
  scale_color_manual(values= c("red", "black","grey"))+
  ggtitle("tsne")

p.tsne= ggplot(tsne.plot, aes(x= V1, y= V2, color = htx))+
  geom_point()+
  scale_color_manual(values= c("red", "black","grey"))+
  ggtitle("tsne")

p.tsne


