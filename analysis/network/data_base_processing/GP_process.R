## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-23
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
##  get GO ontology based network representation
## ---------------------------

library(ontologyIndex)
library(ontologySimilarity)
library(GOSemSim)
library(tidyverse)


GO = get_ontology("http://purl.obolibrary.org/obo/go.obo")
GO_slim= get_ontology("http://purl.obolibrary.org/obo/go/go-basic.obo")

hp = get_human_heart_proteome()
hg= get_GTEX_heart_genes()

ggVennDiagram::ggVennDiagram(list(hp, hg))
genes= union(hp, hg)

GO.mf <- godata( ont="MF")
GO.bp <- godata( ont="BP")

genes <- c("CDC45", "MCM10", "CDC20", "NMU", "MMP1")

geneSim(genes, semData=hsGO, measure="Wang", combine="BMA", verbose=FALSE)


hetnet= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/HETNET_edges.rds")

GO.mf <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="MF", computeIC=TRUE)
genes= unique(c(hetnet[[2]]$nodeA, hetnet[[2]]$nodeB))
# Lin
pairwise_map = mgeneSim(genes, semData=GO.mf, measure="Resnik", combine="BMA", verbose=TRUE)



#######
information_content <- descendants_IC(GO)
?get_sim_grid

sim_maT = get_term_sim_mat(GO,  information_content = information_content)
sim_mat <- get_sim_grid(ontology=GO,
                        term_sets=genes,
                        information_content = information_content,
                        term_sim_method = "lin")
