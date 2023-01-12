## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-11-18
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
##  DEA
library(tidyverse)
library(biomaRt)
library(edgeR)
library(limma)
library(ggrepel)

gex.obj= readRDS("output/gex.obj.bulk.rds")

#camkk2.dko= gex.obj$meta %>% filter(CaMK2!= "CON") %>% pull(Sample_ID)

gex.obj$meta
# funcs ---------------------------------------------------------------------------------------


get_toptable= function(g.type = "J",
                       trtmnt = c("HFD.LNAME", "CON"),
                       contrast.name,
                       remove.DKO = T){
  meta= gex.obj$meta
  #exp= gex.obj$exp

  if(remove.DKO){
    meta = meta %>% filter(CaMK2== "CON")

    }

  meta= meta %>% filter(NNT %in% g.type,
                                Diet %in% trtmnt)

  #NNT= factor(meta$NNT, levels= c("WT", "S100a1", "Strit", "Strit_S100a1"))
  genotype= factor(meta$NNT, levels= g.type)

  #treatment= factor(meta$treatment, levels=c("Sham", "TAC"))
  treatment= factor(meta$Diet, levels= trtmnt)

  if(length(g.type)==1){
    design = model.matrix(~0+treatment)
  }else if(length(trtmnt)==1){
    design = model.matrix(~0+genotype)
  }
  print(design)
  #line= factor(meta$lane_id)
  print(head(gex.obj$exp[,meta$Sample_ID]))

  fit = lmFit(gex.obj$exp[,meta$Sample_ID], design)
  #colnames(gex.obj$exp)

  # #Define contrasts
  # cont.matrix = makeContrasts(c.name= contrast.name,
  #                             levels=design)

  #makeContrast doesnt accept string, so weird work around here:
  cmd <- paste("cont.matrix <- makeContrasts(", contrast.name, ", levels = design)", sep =
                 '"')
  eval(parse(text = cmd))

  print(cont.matrix)
  #Empirical Bayes
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)

  #obtain differentially expressed genes
  DE_results2 = as.data.frame(topTable(fit2,adjust.method = "fdr",number = Inf)) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(desc(abs(t))) %>%
    as_tibble()%>%
    mutate(gene= str_replace_all(gene, " ", ""))


}

# perform DEA on different contrasts: ---------------------------------------------------------

J_HFD= get_toptable(g.type = c("J"),
                  trtmnt = c("HFD.LNAME", "CON"),
                  contrast.name = "treatmentHFD.LNAME-treatmentCON",
                  remove.DKO = T

)

N_HFD= get_toptable(g.type = c("N"),
                    trtmnt = c("HFD.LNAME", "CON"),
                    contrast.name = "treatmentHFD.LNAME-treatmentCON",
                    remove.DKO = F

)

HFD_JN= get_toptable(g.type = c("N", "J"),
                    trtmnt = c("HFD.LNAME"),
                    contrast.name = "genotypeJ-genotypeN"
)

CON_JN= get_toptable(g.type = c("N", "J"),
                     trtmnt = c("CON"),
                     contrast.name = "genotypeJ-genotypeN"
)
library(rstatix)
HFD= list("J"= J_HFD,
          "N"= N_HFD,
          "CON"= CON_JN,
          "HFD"= HFD_JN)

saveRDS(HFD, "output/dea.results.list.rds")
HFD = readRDS("output/dea.results.list.rds")


# plot basic overview -------------------------------------------------------------------------


# plot for N strain ---------------------------------------------------------------------------

#PCA
samps= gex.obj$meta %>% filter(NNT == "N",
                               CaMK2 == "CON")%>%
  pull(Sample_ID)

PCA <- prcomp(t(gex.obj$exp[,samps]) ,center = TRUE, scale. = F)

plot.pca = PCA$x %>%
  as.data.frame %>%
  rownames_to_column("Sample_ID") %>%
  as_tibble()%>%
  left_join(gex.obj$meta)

p.pca1 = ggplot(plot.pca,aes(x= PC1, y= PC2,color = Diet))+
  geom_point(size=2)+
  theme_minimal()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))+
  theme_bw()+
  theme(axis.text = element_text(size = 11, color ="black"),
        rect = element_rect(fill = NA, size = 1))
p.pca1

p.pca2 = ggplot(plot.pca,aes(x= PC3, y= PC4,color = Diet))+
  geom_point(size=2)+
  theme_minimal()+
  labs(x= paste0("PC3 (",as.character(round(PCA$sdev[3]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC4 (",as.character(round(PCA$sdev[4]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))+
  theme_bw()
p.pca2
p.pca1

pdf("output/figures/supp/bulk_dasetal_report.pdf",
    width= 5,
    height= 5)
p.pval.dist
p.volcano
p.pca
dev.off()

cowplot::plot_grid(p.pval.dist, p.volcano, p.pca)

p1=HFD$N %>%
  mutate(sig= ifelse(P.Value<.05 & logFC>1.5, "up",
                     ifelse(P.Value<.05 & logFC< -1.5, "dn", "")),
         label = ifelse( sig== "", "", gene)
  )%>%
  ggplot(aes(x=logFC, y= -log10(P.Value), color =sig, label= label))+
  geom_point()+
  scale_color_manual(values=c("darkgrey", "red", "blue"))+
  geom_text_repel(max.overlaps = 20, show.legend = F)+
  theme_bw()+
  theme(axis.text = element_text(size = 11, color ="black"),
        rect = element_rect(fill = NA, size = 1),
        legend.position = "none")

null_model <- pnorm(rnorm(length(HFD$N$gene)))
qq.plot= cbind("null"= sort(null_model), "p.val"= sort(HFD$N$P.Value))%>%
  as_tibble()%>%
  ggplot(aes(x= null, y= p.val))+
  geom_point()+
  xlim(c(1,0))+
  ylim(c(1,0))+
  geom_abline(slope = 1, intercept= 0)+
  theme_bw()+
  theme(axis.text = element_text(size = 11, color ="black"),
        rect = element_rect(fill = NA, size = 1),
        legend.position = "none")


# enrich it -----------------------------------------------------------------------------------

library(enrichR)


## prepare DBs

dbs <- listEnrichrDbs()

GO_dbs= dbs$libraryName[grepl("GO_", x = dbs$libraryName)]
GO_dbs= GO_dbs[grepl("2021", x = GO_dbs)]

## function to apply enrichr. for GO terms to a list of gene sets
enrich.wrap= function(sets){

  enrich.r= lapply(sets, function(x){
    y= enrichr(genes = x,
               databases =GO_dbs)

  })


  bio.process= map(names(enrich.r), function(x){
    enrich.r[[x]]$GO_Biological_Process_2021%>% mutate(cluster= x)
  })%>% do.call(rbind, .)%>% as_tibble()

  mol.func= map(names(enrich.r), function(x){
    enrich.r[[x]]$GO_Molecular_Function_2021%>% mutate(cluster= x)
  })%>% do.call(rbind, .)%>% as_tibble()

  return(list("mf"= mol.func,
              "bp"= bio.process))
}

r= enrich.wrap(sets=list("up"= HFD$N %>% filter(P.Value<0.05 & logFC >0)%>% pull(gene),
                       "dn"= HFD$N %>% filter(P.Value<0.05 & logFC <0)%>% pull(gene)
                       )
            )

clean_names= function(df, col= "Term"){
  df[[col]]=  gsub(pattern = "\\(.*",replacement = "",x = df[[col]])
  df
}
r$bp= clean_names(r$bp)
r$mf= clean_names(r$mf)
p.bp= r$bp %>%
  filter(Adjusted.P.value<0.02)%>%
  mutate(OR_signed = ifelse(cluster== "up", Odds.Ratio, Odds.Ratio*-1))%>%
  ggplot(., aes(x= OR_signed,y =  reorder(Term,Odds.Ratio), fill = -log10(Adjusted.P.value)))+
  geom_col()+
  theme_bw()+
  geom_col(col = "black")+
  scale_fill_gradient(low= "white", high = "darkred")+
  theme(axis.text = element_text(size = 11, color ="black"),
        rect = element_rect(fill = NA, size = 1),
        legend.position = "right")+
  labs(x= "Odds ratio * sign of regulation",
       y= "",fill = "-log10(q-val)")

p.mf= r$mf %>% filter(Adjusted.P.value<0.2)%>%
  mutate(OR_signed = ifelse(cluster== "up", Odds.Ratio, Odds.Ratio*-1))%>%
  ggplot(., aes(x= OR_signed,y =  reorder(Term,Odds.Ratio), fill = -log10(Adjusted.P.value)))+
  geom_col(col = "black")+
  scale_fill_gradient(low= "white", high = "darkred")+
  theme_bw()+
  theme(axis.text = element_text(size = 11, color ="black"),
        rect = element_rect(fill = NA, size = 1),
        legend.position = "right")+
  labs(x= "Odds ratio * sign of regulation",
       y= "", fill = "-log10(q-val)")

px=cowplot::plot_grid(p.bp, p.mf, ncol = 1 ,labels =  c("D", "E"),
                      rel_widths = c(1,0.8),
                      align= "v")
pdf("output/bulk_overview_GO.pdf",
    width= 8,
    height=8)
px
dev.off()


py= cowplot::plot_grid(p.pca1, qq.plot, p1,
                       rel_widths = c(1,0.8),
                       align = "v",
                       labels= "AUTO")

pdf("output/bulk_overview.pdf",
    width= 8,
    height=8)
py
dev.off()
# plot ----------------------------------------------------------------------------------------

plot.v= function(top.table){
  top.table %>%
    mutate(sig= ifelse(adj.P.Val<.05, "y", "n"))%>%
    ggplot(aes(x=logFC, y= -log10(P.Value), color =sig))+
    geom_point()
}


plot.h= function(top.table){
  top.table %>%
    ggplot(aes(x=P.Value))+
    geom_histogram()
}



pls= map(HFD, plot.h)
cowplot::plot_grid(plotlist = pls)

pls= map(HFD, plot.v)
cowplot::plot_grid(plotlist = pls)


plot.gene.exp= function(gene, excl = T){
  x= enframe(as.matrix(gex.obj$exp)[gene, ],
          name = "Sample_ID",
          value = "exp")%>%
    full_join(gex.obj$meta, by= "Sample_ID")

  if(excl){x =x %>% filter(CaMK2!= "DKO")}
  x%>%
    ggplot(., aes(x= Group, y= exp, color= NNT,shape= Diet, size= CaMK2 ))+
    geom_point()+
    ggtitle(gene)

}
plot.gene.exp(" Postn")
plot.gene.exp(" Acta2")
plot.gene.exp(" Fap")
plot.gene.exp(" Col1a1")
plot.gene.exp(" Col1a2")
plot.gene.exp(" Angptl4")
plot.gene.exp(" Myh7")


# ORA -----------------------------------------------------------------------------------------
##
## gene predictions from comorbidities
gene.p=read.csv("output/predicted_HF_genes.csv")%>% as_tibble()%>%
  arrange(rank.hfpef.prio)

## fibroblast sigs scRNAseq
gene_sigs= readRDS("../scell_hfpef/output/fib_integration/marker_list/DEG_per_study_in_fibs_SET_downsampled.rds")

## full tissue Das. et al 2019
hfpef.h= readRDS("../scell_hfpef/data/_Hfpef_bulk_study/HF_pEF_processed.rds")

## reheat data can be downloaded https://zenodo.org/record/3797044#.XsQPMy2B2u5
reheat= readRDS("/home/jan/R-projects/HF_meta-analysis/data/shiny/directed_ranking.rds")
reheat.p= readRDS("/home/jan/R-projects/HF_meta-analysis/data/shiny/fisher_rank.rds")

reheat.df= enframe(reheat, name= "gene", value = "t")

reheat.df= enframe(reheat.p, name= "gene", value = "t")%>% mutate(t= -log10(t))

df1= comp_feats(HFD$N, pathways = gene_sigs$total, human = F)
df2= comp_feats(HFD$J, pathways = gene_sigs$total, human = F)


comp_feats= function(toptable, gene.p,pathways, human= F, absolute= F){
  toptable= toptable %>% mutate(gene= str_replace_all(gene, " ", ""))
  if(human){
    toptable= toptable %>% mutate(gene= str_to_upper(gene))
  }

  if(absolute){
    vec= abs(toptable$t)
    }else{
      vec= (toptable$t)
  }

  names(vec)= toptable$gene
  vec= sort(vec, decreasing = T)
  fgsea::fgsea(stat= vec, pathways= pathways)
}

GSE_analysis = function(geneList,Annotation_DB, n= NULL){
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
    if(is.null(n)){
      ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)

    }else{
      ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=n, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)

    }
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


# loop over cut-offs:

res= map(seq(25,100, by = 25), function(x){
  hfpef= gene.p %>% arrange(rank.hfpef.prio) %>% slice_head(.,n= x ) %>% pull(gene)
  hfref= gene.p %>% arrange(rank.hfref.prio) %>% slice_head(.,n= x)%>% pull(gene)

  gene.pred= list("hfpef"= hfpef,
                  "hfref"= hfref)

  #enrich genes in HFD N
  df1= comp_feats(HFD$N, pathways =gene.pred, human  = T , absolute = F)
  df2= comp_feats(HFD$J, pathways =gene.pred, human  = T , absolute = F)

  df3= comp_feats(hfpef.h$DEA,
                  pathways = gene.pred, human = F, absolute = F)

  df4= comp_feats(reheat.df,
                  pathways = gene.pred, human = F, absolute = F)

  df = rbind(df1%>% mutate(model = "mouse_HFD_N"),
             df2%>% mutate(model = "mouse_HFD_J"),
             df3 %>% mutate(model = "human_HFpEF"),
             df4 %>% mutate(model = "human_HFrEF")
  )%>% mutate(cutoff= x)

  return(df)
})

full.df= do.call(rbind, res) %>%
  #mutate(padj= p.adjust(pval))%>%
  mutate(sig = ifelse(padj<0.05, "*", ""),
         sig = ifelse(padj<0.01, "**", sig),
         sig = ifelse(padj<0.001, "***", sig),
         pathway= ifelse(pathway== "hfpef", "HFpEF", "HFrEF"))

p.enrich= full.df %>% filter(grepl("mouse" , model))%>%
  ggplot(., aes(x= cutoff, y= pathway, fill = -log10(pval), label  = sig))+
  facet_grid(rows=vars(model))+
  geom_tile(color= "darkgrey")+
  geom_text(aes(text= sig))+
  scale_fill_gradient2(low= "blue", mid= "white", high = "red")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text = element_text(color= "black", size = 11))+
  labs(y= "predicted_genes")

p.enrich2= full.df %>% filter(grepl("human" , model))%>%
  mutate(cutoff= as.factor(cutoff))%>%
  ggplot(., aes(x= cutoff, y= pathway, fill = -log10(padj), label  = sig))+
  facet_grid(rows=vars(model))+
  geom_tile(color= "darkgrey")+
  geom_text(aes(text= sig))+
  scale_fill_gradient2(low= "blue", mid= "white", high = "red")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text = element_text(color= "black", size = 11))+
  labs(y= "predicted gene set")


p.enrich2


p.enrich3= full.df %>%
  filter(grepl("mouse" , model) & grepl("N" , model))%>%
  mutate(cutoff= as.factor(cutoff))%>%
  ggplot(., aes(x= cutoff, y= pathway,
                fill = -log10(padj),

                label  = sig))+
  facet_grid(rows=vars(model))+
  geom_tile(color= "black")+
  geom_text(aes(text= sig))+
  scale_fill_gradient2(low= "blue", mid= "white", high = "red")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text = element_text(color= "black", size = 11))+
  labs(y= "predicted gene set")
p.enrich3

pdf(file = "output/fig.enrich_mice.pdf",
    width= 4.2,
    height=3)
p.enrich3+coord_equal()
dev.off()

pdf(file = "output/fig.enrich.human.pdf",
    width= 4,
    height= 3)
p.enrich2
dev.off()

pdf(file = "output/fig.enrich.human.pdf",
    width= 4,
    height= 3)
p.enrich2
dev.off()

####### volcanos
hfpef= gene.p %>% arrange(rank.hfpef.prio) %>% slice_head(.,n= 100 ) %>% pull(gene)
hfref= gene.p %>% arrange(rank.hfref.prio) %>% slice_head(.,n= 100)%>% pull(gene)
library(ggrepel)


volc.plot.df= HFD$N %>%
#volc.plot.df= N_HFD%>%
  mutate(sig= ifelse(adj.P.Val<.1, "y", "n"))%>%
  mutate(gene = str_replace_all(gene, " ", ""))%>%
  mutate(pred= ifelse(gene %in% str_to_title(hfpef), "HFpEF",
                        ifelse(gene %in% str_to_title(hfref), "HFrEF", "none")),
         pred= factor(pred, levels = c("HFrEF","HFpEF","none" )),
         #lab= ifelse(pred != "bg", gene, "")
         )%>%
  arrange(desc(pred))

hits= volc.plot.df%>%
  filter(pred %in% c("HFpEF"))%>%
  slice_max(order_by = abs(logFC), n=25)%>%
  pull(gene)

p.volc = volc.plot.df %>%
  mutate(lab= ifelse(gene %in% hits, gene, ""))%>%
  ggplot(aes(x=logFC, y= -log10(P.Value), color =pred, label = pred))+
  geom_point()+
  scale_color_manual(values=c( "#E54F6D", "#451F55","darkgrey"))+
  theme_bw()+
  geom_label_repel(aes(label=lab ),alpha= 0.8,
                   max.overlaps = 1040,segment.alpha= 0.5,
                   show.legend = FALSE)+
  theme(axis.text = element_text(size= 11, color ="black"),
        rect = element_rect(size= 1, fill = NA))+
  labs(col= "Gene\npredictions")+
  guides(colour = guide_legend(override.aes = list(size=4)))

p.volc
pdf(file = "output/fig.mice.volc.predicted",
    width = 6,
    height= 5)
p.volc
dev.off()

vec= (HFD$N$t)

names(vec)= HFD$N$gene
names(vec)= str_to_upper(names(vec))

vec= sort(vec, decreasing = T)
hfpef= gene.p %>% arrange(rank.hfpef.prio) %>% slice_head(.,n= x ) %>% pull(gene)
hfref= gene.p %>% arrange(rank.hfref.prio) %>% slice_head(.,n= x)%>% pull(gene)

gene.pred= list("hfpef"= hfpef,
                "hfref"= hfref)

fgsea::fgsea(stat= vec, pathways= gene.pred)

p.pef= fgsea::plotEnrichment(stats= vec, pathway= gene.pred$hfpef)
p.ref = fgsea::plotEnrichment(stats= vec, pathway= gene.pred$hfref)

px= cowplot::plot_grid(p.pef, p.ref)

pdf(file = "output/fig.mice.enrichplot.predicted",
    width = 4.5,
    height= 2.5)
px
dev.off()



## human
volc.plot.df2= hfpef.h$DEA %>%
  #volc.plot.df= N_HFD%>%
  mutate(sig= ifelse(adj.P.Val<.1, "y", "n"))%>%
  mutate(gene = str_replace_all(gene, " ", ""))%>%
  mutate(pred= ifelse(gene %in% (hfpef), "HFpEF",
                      ifelse(gene %in% (hfref), "HFrEF", "bg")),
         pred= factor(pred, levels = c("HFrEF","HFpEF","bg" )),
         #lab= ifelse(pred != "bg", gene, "")
  )%>%
  arrange(desc(pred))

hits2= volc.plot.df2%>%
  filter(pred %in% c("HFpEF"))%>%
  slice_max(order_by = abs(logFC), n=25)%>%
  pull(gene)

p.volc2 = volc.plot.df2 %>%
  mutate(lab= ifelse(gene %in% hits2, gene, ""))%>%
  ggplot(aes(x=logFC, y= -log10(P.Value), color =pred, label = pred))+
  geom_point()+
  scale_color_manual(values=c( "green", "blue","darkgrey"))+
  theme_bw()+
  geom_label_repel(aes(label=lab ), max.overlaps = 1040,show.legend = FALSE)+
  theme(axis.text = element_text(size= 11, color ="black"))

cowplot::plot_grid(p.volc2, p.volc)

pdf("output/fig.volcano_mouse.pdf",
    width= 6,
    height= 6)
p.volc
p.volc2
dev.off()



l.edge= full.df%>%
  filter(cutoff== 51,
         pathway== "hfpef")%>%
  #       model %in% c( "mouse_HFD_N", "human_HFpEF"))%>%
  pull(leadingEdge)
names(l.edge)= c("N" , "J", "pef", "ref")
l.edge
ggvenn::ggvenn(l.edge[[1]], l.edge[[2]])
intersect(l.edge[[1]], l.edge[[2]])
VennDiagram::venn.diagram(l.edge)

HFD$N%>% arrange(desc(t))%>% filter(gene %in% str_to_title(l.edge$N))
HFD$J%>% arrange(desc(t))%>% filter(gene %in% str_to_title(l.edge$J))

hfpef= gene.p %>% arrange(rank.hfpef.prio) %>% slice_head(.,n= 50 ) %>% pull(gene)
hfref= gene.p %>% arrange(rank.hfref.prio) %>% slice_head(.,n= 50)%>% pull(gene)

gene.pred= list("hfpef"= hfpef,
                "hfref"= hfref)

comp_feats(toptable = HFD$N, pathways = list("hfpef"= hfpef,
                                  "hfref"= hfref), human = T)




df1= comp_feats(HFD$N, pathways =gene.pred, human  = T , absolute = F)
# df2= comp_feats(HFD$J, pathways =gene.pred, human  = T , absolute = T)
#
# df3= comp_feats(hfpef.h$DEA%>% mutate(gene= str_to_title(gene)),
#                 pathways = gene_sigs$total, human = F)
df4= comp_feats(hfpef.h$DEA,
                pathways = gene.pred, human = F, absolute = F)


#df2= comp_feats(HFD$J, pathways = gene_sigs$total, human = F)

df = rbind(df1%>% mutate(model = "murine"),
      df4 %>% mutate(model = "human")
)

df= df%>%mutate(sig = ifelse(padj<0.05, "*", ""),
            sig = ifelse(padj<0.01, "**", sig)
)
df%>%
  #ggplot(., aes( x= strain, y= pathway, fill = -log10(padj), label = sig))+
  ggplot(., aes( x= model, y= pathway, fill = NES, label = round(padj, 4)))+
  geom_tile()+
  geom_text(aes(text= sig))+
  scale_fill_gradient(low= "white", high = "red")+
  labs(y= "predicted genes",
       x= "Bulk RNAseq")+
  theme_minimal()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text = element_text(color= "black", size = 11))


hfpef.h$DEA %>%
  left_join(HFD$J%>% mutate(gene= str_to_upper(gene)),by= "gene")%>%
  mutate(gene.predicted= ifelse(gene %in% gene.pred$hfpef, "hfpef",
                                ifelse(gene %in% gene.pred$hfref, "hfref", "n")))%>%
  filter(gene.predicted != "n")%>%
  ggplot(., aes(x= t.x, y= t.y, col= gene.predicted))+
  geom_point()+
  geom_text_re

#ora with fib sigs
gene_sigs.h= lapply(gene_sigs$total, str_to_upper)
GSE_analysis(geneList = gene.pred$hfpef, Annotation_DB = gene_sigs.h, n= 15000)
GSE_analysis(geneList = gene.pred$hfref, Annotation_DB = gene_sigs.h, n= 15000)
GSE_analysis(geneList = gene.pred$hfpef, Annotation_DB = list("r"= names(reheat[1:500])))

