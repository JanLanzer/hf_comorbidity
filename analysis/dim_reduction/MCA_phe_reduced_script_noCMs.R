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
library(factoextra)
library(ggrepel)
library(ggExtra)
library(umap)


# Read data and prepare ---------------------------------------------------

# read tables.
source("~/GitHub/RWH_analysis/scripts/utils.R")

directory= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/"

icd = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")

#load features
phecodes= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/topfreq_disease.rds")

pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")
map(pids.list,length)

pheno.data= readRDS(file= "T:/fsa04/MED2-HF-Comorbidities/data/processed_data/full_clinic_info.rds")
pheno.data2= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/patient_metadata_2022.rds")


# Functions ---------------------------------------------------------------

## function to prepare df as input for MCA()
#' @param df,  df in tidy format with 2 columns (pid, PheCode)

MCAinput = function(df){
  #transform to count matrix with disease as columns, patients as rows and 0 and 1 defining the disease profile
  df= split(as.character(df$PheCode), df$pid) %>%
    lapply(function(x) strsplit(x, "-")) %>%
    mtabulate()%>%
    matrix2df("pid") %>%
    dplyr::as_data_frame(check.names = FALSE) %>%
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
                      repel = F, # Avoids text overlapping (slow if many point)
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
icd_red = icd %>%
  drop_na %>%
  distinct(pid, PheCode) %>%
  filter(pid %in% c(pids.list$hfpef, pids.list$hfref, pids.list$hfmref),
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


# add clinical endpoints:
endpoint= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/pids_endpoints.rds")

ind.df= ind.df%>%
  mutate(htx= ifelse(ID %in% endpoint$htx, "yes", "no"),
         intu= ifelse(ID %in% endpoint$intubation, "yes", "no"),
         defi= ifelse(ID %in% endpoint$defi, "yes", "no"),
         pci= ifelse(ID %in% endpoint$pci, "yes", "no"))

#add summary tableS=
sum.t= readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/patient_metadata_2022.rds")

sum.tfilt= sum.t %>%
  filter(pid %in% ind.df$ID) %>%
  distinct(pid,
           patient_cohort,
           sex,
           age.at.icd,
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
           #icd10gm_count,
           charlson_score,
           elixhauser_score,
           ef.min,
           edd.max,
           ea.min,
           ea.max,
           ee.max)


ind.df= ind.df %>% left_join(sum.tfilt%>% dplyr::rename(ID= pid)%>% mutate(ID= as.character(ID)), by= "ID")

#add.patient.cohort:
ind.df= ind.df%>%   mutate(hf.type = ifelse(ID %in% pids.list$hfref,
                                            "hfref",
                                            ifelse(ID %in% pids.list$hfpef,
                                                   "hfpef",
                                                   ifelse(ID %in% pids.list$hfmref,
                                                          "hfmref",
                                                          "unlabeled"))))
## run wilcox to see if it associates dim1 with hf

dim_= paste0("Dim", seq(1:ndims-1))

df= ind.df
## loop over categorical variables with two levels and perform wilcox
cat.vars= c( "sex",  "intu", "defi", "pci")

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

  associated_hf =cbind(df.variance, pvals) %>%
    as_tibble%>%
    mutate(sig= ifelse(pvals<0.05, "sig", "ns"))%>%
    group_by(sig)%>%
    summarise(s= sum(variance.percent)) %>%
    mutate(var= x)

})

cat.df= do.call(rbind, tested.vars)%>%filter(sig=="sig")

#redo for hf.type (4 levels)
###
hf.var= data.frame(matrix(nrow = length(unique(df$hf.type)),
                          ncol = length(unique(df$hf.type))
                          ),
                   row.names = unique(df$hf.type)
                  )
colnames(hf.var)= unique(df$hf.type)


#old way (just differnet data structure, same analysis)
hftype= lapply(unique(df$hf.type), function(x){
  sapply(unique(df$hf.type), function(y){
    if(x==y){return(0)}

    print( paste(x,"vs.",  y))
    hf_t = df %>% filter(hf.type== x)
    hf_f = df %>% filter(hf.type== y)

    pvals= map(dim_, function(x){
      test. = wilcox.test(hf_t %>% pull(x),
                          hf_f%>% pull(x),
                          alternative = "two.sided")
      test.$p.value
    }) %>%unlist()

    #combine with variance explained:
    eig.val <- get_eigenvalue(mca.res)

    df.variance= eig.val[1:length(pvals),]

    associated_hf =cbind(df.variance, pvals) %>%
      as_tibble%>%
      mutate(sig= ifelse(pvals<0.05, "sig", "ns"))%>%
      group_by(sig)%>%
      summarise(s= sum(variance.percent))

    var.sum= associated_hf%>% filter(sig== "sig")%>% pull(s)
    print(var.sum)

    #add to matrix
    hf.var[x,y]= var.sum

    #names(var.sum)= paste(x,"vs.",  y)
    return(c(var.sum, paste(x,"vs.",  y)))
  })
})



cat.df= rbind(cat.df,
      c("sig",hftype[[1]]$hfpef[1], "HFrEF v. HFpEF" ),
      c("sig",hftype[[1]]$hfmref[1], "HFmrEF v. HFrEF" ),
      c("sig",hftype[[2]]$hfmref[1], "HFmrEF v. HFpEF" )
      )%>%
  mutate(s = as.numeric(s))


p.cat= cat.df%>%filter(sig=="sig")%>%
  ggplot(., aes(x= reorder(var, s), y= s))+
  geom_col()+
  labs(y= "percentage of variance explained (%)",
       x= "tested covariate")+
  theme(legend.position =  "none")#+
  #theme_minimal_hgrid()

## loop over continous variables, perform cor.test
cont.vars= c(
             "median.bnp",
             "mean.sys",
             "mean.dias",
              "age.at.icd",
             "elixhauser_score",
             "median.BMI",
             "median.LDL",
             "median.Chol",
             "median.HDL",
           "median.Trigs")


tested.vars.cont= map(cont.vars, function(x){
  print(x)

  pvals= map(dim_, function(y){
    test. = ind.df %>% select(!!as.symbol(x),
                              !!as.symbol(y))%>%
      drop_na()
    #cor.test(test.[[x]], test.[[y]])$p.value
    fit= lm(formula= paste0(y, " ~ ", x), data= test.)
    summary(fit)$coefficients[2,4]
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

df.explained_V= rbind(do.call(rbind, tested.vars.cont) , cat.df)%>%
  mutate(var.type= ifelse(var %in% cat.vars, "categorical","continuous"))%>%
  mutate(var.type = ifelse(grepl("HF", var), "categorical",var.type))%>%
  filter(sig=="sig")


#clean var names

df.explained_V= df.explained_V%>% mutate(var= ifelse(var == "sex", "Sex", var),
                         var= ifelse(var == "age.at.icd", "Age", var),
                         var= ifelse(var == "mean.sys", "Systolic RR", var),
                         var= ifelse(var == "mean.dias", "Diastolic RR", var),
                         var= ifelse(var == "median.Chol", "Cholesterol", var),
                         var= ifelse(var == "elixhauser_score", "Elixhauser", var),
                         var= ifelse(var == "pci", "PCI", var),
                         var= ifelse(var == "defi", "Device implant", var),
                         var= ifelse(var == "median.Trigs", "Triglycerides", var),
                         var= ifelse(var == "intu", "Intubation", var),
                         var= ifelse(var == "median.hba1c", "HbA1c", var),
                         var= ifelse(var == "median.BMI", "BMI", var),
                         var= ifelse(var == "median.HDL", "HDL", var),
                         var= ifelse(var ==  "median.bnp","NT-ProBNP", var),
                         var= ifelse(var == "median.LDL", "LDL", var))


p.explained_V=
  df.explained_V%>%
  ggplot(., aes(x= reorder(var, s), y= s, fill = var.type))+
  geom_hline(yintercept = c(20,30,40,50, 60,70), color= "darkgrey")+
  geom_col()+
  scale_fill_manual(values=cols.nice)+
  labs(y= "Estimated explained variance (%)",
       x= "Tested covariate",
       fill= "Variable type")+
  theme(legend.position =  "none")+
  #scale_y_continuous(breaks = seq(0, 80, 10))+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(color= "black"),
        legend.position = "right",
        legend.text = element_text(size= 8),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
        )+
  ylim(c(0,70))

p.explained_V


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/explained_V.pdf",
    height = 4, width= 6)
unify_axis(p.explained_V)

dev.off()


#DiM red---------------------------------------------------------



source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_classifier_ML.R")

set.seed(10)


#### UMAP

df= mca.df%>% mutate_all(as.numeric)
umap_res= umap(df,preserve.seed= T)

umap.plot = as_tibble(cbind(as_tibble(umap_res$layout), pid=  rownames(df))) %>%
  mutate(hf = ifelse(pid   %in% pids.list$hfref,
                                "HFrEF",
                                ifelse(pid   %in% pids.list$hfpef,
                                       "HFpEF",
                                       ifelse(pid   %in% pids.list$hfmref,
                                              "HFmrEF",
                                              "unlabeled"))))%>%
  mutate(hf= factor(hf, levels= c("HFpEF","HFmrEF", "HFrEF", "unlabeled"))) %>%
  left_join(sum.tfilt%>% mutate(pid= as.character(pid)), by= "pid")%>%
  mutate(HTX= ifelse(pid %in% endpoint$htx, "yes", "no"),
         Intubation= ifelse(pid %in% endpoint$intubation, "yes", "no"),
         Device_implant= ifelse(pid %in% endpoint$defi, "yes", "no"),
         PCI= ifelse(pid %in% endpoint$pci, "yes", "no"))

p.umap.hf =ggplot(umap.plot, aes(x= V1, y= V2, color = hf))+
  geom_point(alpha= 0.7, size= .5)+
  #scale_color_manual(values=c( col.set[2], "#7FC6A4", col.set[3]))+
  scale_color_manual(values= cols.nice[-3])+
  labs(color="HF cohort",
       x= "",
       y = "")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))

p.umap.hf

p.umap.defi =ggplot(umap.plot, aes(x= V1, y= V2, color = Device_implant))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_manual(values= cols.nice[-2])+
  labs( x= "",
       y = "",
       col = "Device Implant")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))

p.umap.intu =ggplot(umap.plot, aes(x= V1, y= V2, color = Intubation ))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_manual(values= cols.nice[-2])+
  labs( x= "",
        y = "",
        col= "Intubation")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))


p.umap.sex =ggplot(umap.plot%>% mutate(sex=factor(sex, levels=c("m", "f"))), aes(x= V1, y= V2, color = sex))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_manual(values= cols.nice[-2])+
  labs( x= "",
        y = "",
        col= "Sex")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))

p.umap.pci=ggplot(umap.plot, aes(x= V1, y= V2, color = PCI))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_manual(values= cols.nice[-2])+
  labs( x= "",
        y = "",
        col= "PCI")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))

##

p.umap.age =ggplot(umap.plot, aes(x= V1, y= V2, color = age.at.icd))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_gradient(low= cols.nice[3], high= cols.nice[1])+
  labs( x= "",
        y = "",
        col = "Age")+
  theme_classic()

p.umap.charlson =ggplot(umap.plot, aes(x= V1, y= V2, color = elixhauser_wscore  ))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_gradient(low= cols.nice[3], high= cols.nice[1])+
  labs( x= "",
        y = "",
        col = "Elixhauser Index")+
  theme_classic()


p.umap.hba1c =ggplot(umap.plot, aes(x= V1, y= V2, color = median.hba1c))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_gradient(low= cols.nice[3], high= cols.nice[1])+
  labs( x= "",
        y = "")+
   theme_classic()



p.umap.bmi =ggplot(umap.plot%>% filter(median.BMI<50), aes(x= V1, y= V2, color = median.BMI))+
  geom_point(alpha= 0.7, size= .5)+
  scale_color_gradient(low= cols.nice[3], high= cols.nice[1])+
  labs( x= "",
        y = "")+
  theme_classic()


p.umap.bnp =ggplot(umap.plot, aes(x= V1, y= V2, color = log10(median.bnp)))+
  geom_point(alpha= 0.7, size= 1.)+
  scale_color_gradient(low= cols.nice[3], high= cols.nice[1])+
  labs( x= "",
        y = "")+
  theme_classic()+
  guides(fill  = guide_legend(override.aes = list(size = 0.1)))

p.umap.bnp
p.umap.hba1c
main.umaps= cowplot::plot_grid( p.umap.defi+theme(axis.text = element_blank(),
                                                  axis.ticks = element_blank(),
                                                  panel.border = element_rect(colour = "black", fill=NA, size=1)),
                                p.umap.charlson+theme(axis.text = element_blank(),
                                                       axis.ticks = element_blank(),
                                                       panel.border = element_rect(colour = "black", fill=NA, size=1)),
                                p.umap.pci+theme(axis.text = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 panel.border = element_rect(colour = "black", fill=NA, size=1)),
                                p.umap.hf+theme(axis.text = element_blank(),
                                                axis.ticks = element_blank(),
                                                panel.border = element_rect(colour = "black", fill=NA, size=1)),
                                p.umap.sex+theme(axis.text = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 panel.border = element_rect(colour = "black", fill=NA, size=1)),
                                p.umap.age+theme(axis.text = element_blank(),
                                                  axis.ticks = element_blank(),
                                                  panel.border = element_rect(colour = "black", fill=NA, size=1)),
                                p.umap.intu+theme(axis.text = element_blank(),
                                                  axis.ticks = element_blank(),
                                                  panel.border = element_rect(colour = "black", fill=NA, size=1)),
                                p.umap.bmi+theme(axis.text = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 panel.border = element_rect(colour = "black", fill=NA, size=1)),
                                ncol = 3, align = "hv")


main.umaps
pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/manuscript/main/umaps_cat.pdf")
#cowplot::plot_grid(p.umap.hf,p.umap.sex, p.umap.htx, p.umap.defi, p.umap.pci)
main.umaps
dev.off()

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/manuscript/main/umaps_cont.pdf")
    #heigh= 3, width= 6)
cowplot::plot_grid(p.umap.age, p.umap.hba1c, p.umap.bmi,p.umap.bnp)
dev.off()

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/umaps_main.pdf",
    height= 9, width= 14)
main.umaps
dev.off()

# umap with single disease feats ------------------------------------------
main_disease= c("585.3", "401.1", "272.13", "250.2", "296.22", "496",
              "280.1","327.3","411.4", "425.1","394.2"    )
main_disease= c("394.2","395.1"   )
pls= map(main_disease, function(feature){
  pids.x= icd_red %>% filter(PheCode == feature) %>% pull(pid)
  pheno.x= icd%>%  filter(PheCode == feature) %>% pull(Phenotype)%>% unique()
  umap.plot %>% mutate(dis1= ifelse(pid %in% pids.x, "y", "n"))%>%
    ggplot(., aes(x= V1, y= V2, color = dis1))+
    geom_point(alpha= 0.7, size= .5)+
    scale_color_manual(values= col.set)+
    labs( x= "",
          y = "")+
    theme_classic()+
    ggtitle(paste0(feature,pheno.x))
})

plot_grid(plotlist = pls)
pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/umaps_supp.pdf",
    height= 6, width= 6)
pls
dev.off()

print(pls[[1]])


## tSNE

df= mca.df[!duplicated(mca.df), ]
dim(df)

tsne_res= Rtsne(df, check_duplicates = F, perplexity = 50)

patient.df = tibble(pid= rownames(df)) %>%
  mutate(hf = ifelse(pid   %in% pids.list$hfref,
                     "HFrEF",
                     ifelse(pid   %in% pids.list$hfpef,
                            "HFpEF",
                            ifelse(pid   %in% pids.list$hfmref,
                                   "HFmrEF",
                                   "unlabeled"))))%>%
  mutate(hf= factor(hf, levels= c("HFpEF","HFmrEF", "HFrEF", "unlabeled")))%>%
  left_join(sum.tfilt%>% mutate(pid= as.character(pid)), by= "pid")%>%
  mutate(htx= ifelse(pid %in% endpoint$htx, "yes", "no"),
         intu= ifelse(pid %in% endpoint$intubation, "yes", "no"),
         defi= ifelse(pid %in% endpoint$defi, "yes", "no"),
         pci= ifelse(pid %in% endpoint$pci, "yes", "no"),
         tavi= ifelse(pid %in% pids.oi, "yes", "no"))


tsne.plot = as_tibble(cbind(as_tibble(tsne_res$Y), pid=  rownames(df))) %>%
  left_join(patient.df, by= "pid")

tsne= ggplot(tsne.plot, aes(x= V1, y= V2, color = sex))+
  geom_point()+
  scale_color_manual(values= c("red", "black","grey"))+
  ggtitle("tsne")

p.tsne= ggplot(tsne.plot, aes(x= V1, y= V2, color = pci))+
  geom_point()+
  scale_color_manual(values= c("red", "black","grey"))+
  ggtitle("tsne")

p.tsne

p.tsne= ggplot(tsne.plot, aes(x= V1, y= V2, color = tavi))+
  geom_point()+
  scale_color_manual(values= c("red", "black","grey"))+
  ggtitle("tsne")


# select patients ---------------------------------------------------------
p.umap.age

pids.oi = umap.plot %>% filter(V1 < -3) %>% pull(pid)
pids.oi = icd_red %>% filter(pid %in% htx.cohort$defi)%>% pull(pid)

pids.oi

x= disease_frequencies(pids = pids.oi, icd = icd)
xy= disease_frequencies(pids = unlist(pids.list[2:3])[!unlist(pids.list[2:3]) %in% pids.oi], icd = icd)
icd %>% filter(PheCode =="425.1", !pid %in% pids.oi) %>% distinct(icd4, PheCode)


# ####### full cohort umap ------------------------------------------------




icd_red = icd %>%
  drop_na %>%
  distinct(pid, PheCode) %>%
  filter(pid %in% unlist(pids.list),
         PheCode %in% phecodes)

# icd_red is full dataset to work with. This will now be subsetted to answer specific questions of interest.
length(unique(icd_red$pid))
length(unique(icd_red$PheCode))

#### use all patients, use all diagnostics and analyze impact of HF and Gender as qualitativ variables
## 1) transform full data set into input format for MCA function

mca.df=MCAinput(df = icd_red)

colnames(mca.df) = gsub("^X", "", colnames(mca.df))

df= mca.df%>% mutate_all(as.numeric)
umap_res= umap(df,preserve.seed= T)

umap.plot = as_tibble(cbind(as_tibble(umap_res$layout), pid=  rownames(df))) %>%
  mutate(hf = ifelse(pid   %in% pids.list$hfref,
                     "HFrEF",
                     ifelse(pid   %in% pids.list$hfpef,
                            "HFpEF",
                            ifelse(pid   %in% pids.list$hfmref,
                                   "HFmrEF",
                                   "unlabeled"))))%>%
  mutate(hf= factor(hf, levels= c("HFpEF","HFmrEF", "HFrEF", "unlabeled")))# %>%

p.umap.hf =ggplot(umap.plot, aes(x= V1, y= V2, color = hf))+
  geom_point(alpha= 0.7, size= .5)+
  #scale_color_manual(values=c( col.set[2], "#7FC6A4", col.set[3]))+
  scale_color_manual(values= cols.nice[-3])+
  labs(color="HF cohort",
       x= "",
       y = "")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=4)))

