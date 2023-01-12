## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-11-17
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## RNAseq processing

library(tidyverse)
library(biomaRt)
library(rtracklayer)
library(edgeR)
library(limma)
library(ggrepel)


# Annotations ---------------------------------------------------------------------------------
## gtf from https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus_c57bl6nj/
# gtf <- rtracklayer::import('data/M016_RNA.Seq_Counts/Mus_musculus_c57bl6nj.gtf')
# gtf <- rtracklayer::import('~/Downloads/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.108.gtf.gz')
# gtf_df=as.data.frame(gtf)

gtf2 <- read.table('data/M016_RNA.Seq_Counts/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.108.gtf',
                   header = FALSE, sep = '\t')
mgi= as_tibble(gtf2)
mgi2= mgi%>% separate(V9, sep = ";", into= c(paste0("col", 1:6)))%>%
  select(V2, V3, col1, col3, col4)%>%
  filter(V3 == "gene") %>%
  filter(grepl(col3, pattern = "gene_name"))%>% ## here we filter out for protein coding genes, revisit if transcripts/peseudogenes are of interest
  mutate(col1= str_replace_all(col1, "gene_id ", ""),
         col3= str_replace_all(col3, "gene_name ", ""))

#example for a highly expressed pseudo gene  =
mgi %>% filter(grepl("MGP_C57BL6NJ_G0040439", V9))

#mgi.protein= as_tibble(mig2)%>% filter(V3 == "gene")





# get count table -----------------------------------------------------------------------------

#we will use unstranded data

file.vec= list.files("data/M016_RNA.Seq_Counts/")
samples= c()

samples= map(file.vec, function(i){
  sample.name= paste0("ID_",substr(i, 1, 4))

  sample= read_table(paste0("data/M016_RNA.Seq_Counts/",i), skip = 3)
  colnames(sample)= c("gene", sample.name, "forward.s", "reverse.s")

  return(sample)

})

df= Reduce(function(...) merge(..., by='gene', all.x=TRUE), samples)

df= cbind("gene"= df$gene, df[,grepl("ID", colnames(df))])

df= df %>% left_join(mgi2%>% select(col3, col1) %>% dplyr::rename(gene = col1,
                                                              geneSymbol= col3))

df2= df[,-1] %>%filter(!is.na(geneSymbol))%>% group_by(geneSymbol)%>% summarise_all(mean)

df2= column_to_rownames(.data = df2, "geneSymbol")



# target.file ---------------------------------------------------------------------------------

target= read.csv("data/colData.csv") %>% mutate(Sample_ID= paste0("ID_", Sample_ID))


# filter & normalize --------------------------------------------------------------------------

#order df
df= df2[, match( target$Sample_ID, colnames(df2))]

colnames(df)== target$Sample_ID

group= target$Group
dge <- DGEList(counts=df, group=group)

#detect and remove low expressed gene
keep <- filterByExpr(dge,group = group, min.prop =0.5,min.count = 10)
#?filterByExpr
dge <- dge[keep,,keep.lib.sizes=FALSE]
#dge$counts= dge$counts + 5

dge <- calcNormFactors(dge)
# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)

### save R-objects
voom_count= v$E

gex.obj= list(dge= dge,
              exp= v$E,
              meta= target)

saveRDS(gex.obj, "output/gex.obj.bulk.rds")
gex.obj=readRDS("output/gex.obj.bulk.rds")
#PCA

PCA <- prcomp(t(gex.obj$exp) ,center = TRUE, scale. = F)

plot.pca = PCA$x %>%
  as.data.frame %>%
  rownames_to_column("Sample_ID") %>%
  as_tibble()%>%
  left_join(gex.obj$meta)

p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2,color = Diet, shape= NNT, size= CaMK2))+
  geom_point()+
  theme_minimal()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))+
  geom_text_repel(aes(label= Sample_ID),show.legend = FALSE)
p.pca

pdf("output/figures/supp/bulk_dasetal_report.pdf",
    width= 5,
    height= 5)
p.pval.dist
p.volcano
p.pca
dev.off()

cowplot::plot_grid(p.pval.dist, p.volcano, p.pca)

