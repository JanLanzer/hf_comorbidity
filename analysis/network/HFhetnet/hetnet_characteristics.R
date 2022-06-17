## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-11-24
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## use Pubscore to analyse pubmed bias of genes
## ---------------------------

library(igraph)
library(tidyverse)
library(PubScore)
library(annotate)
library(org.Hs.eg.db)


#get edge lists of hetnet layers:
edge.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/multilayer_edge_list.rds")


source("analysis/utils/utils_network.R")

# add pubscore ------------------------------------------------------------


data("gene2pubmed_db")
genedf= gene2pubmed_db %>%
  drop_na%>%
  mutate(GeneID= as.character(GeneID))

y= getSYMBOL(x = as.character(unique(genedf$GeneID)),data = "org.Hs.eg") %>%
  enframe(name = "GeneID", value = "gene")

gene_literature_count= genedf %>%
  left_join(y) %>%
  as_tibble() %>%
  dplyr::distinct(PubMed_ID, gene) %>% # count only publication once per gene
  group_by(gene) %>%
  count

#node degree_correlation with pubmed2gene score:
degrees= lapply(edge.list$gene, function(x){
  igraph::degree(
    x,
    v = V(x),
    mode = "all",
    loops = TRUE,
    normalized = FALSE
  )

})

cor.res= lapply(degrees, function(x){
  df= enframe(x, name = "gene", value = "degree") %>% left_join(gene_literature_count)
  cor.= cor.test(df$degree, df$n)
})

corr= lapply(cor.res, function(x){
  data.frame("corr"= x$estimate,
             "p"= x$p.value)
})%>% do.call(rbind,.) %>% rownames_to_column("layer")


p.corr=corr %>%
  ggplot(., aes(x= layer, y= corr))+
  geom_col(width= 0, colour= "black", lwd= 2)+
  geom_point(size= 5)+
  coord_flip()+
  labs(y= "Pearson r")+
  theme_cowplot()+
  geom_hline(yintercept = 0, size= 1)

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/manuscript/main/pub_bias_gene_layer.pdf",
    height =4, width= 4)
p.corr
dev.off()



# compare network sizes ---------------------------------------------------

gene_net.sizes= sapply(edge.list$gene, function(x){
  c(length(V(x)),  length(E(x)))
})

disease_net.sizes= sapply(edge.list$disease, function(x){
  c(length(V(x)),  length(E(x)))
})

df.size= cbind(disease_net.sizes, gene_net.sizes)

rownames(df.size)= c("# nodes","# edges")

df.size= df.size %>% as.data.frame() %>% rownames_to_column("feature") %>% pivot_longer(-feature)

ggplot(df.size, aes(x= name, y= value ))+
  facet_grid(rows=vars(feature), scales = "free_y")+
  geom_col()



# compare toplogy ---------------------------------------------------------

gene_net.feats= lapply(names(edge.list$gene), function(x){
  getProps_global(edge.list$gene[[x]], directed = F)%>% mutate(net= x)
})
df.g= do.call(rbind, gene_net.feats)%>% mutate(value = as.numeric(value)) %>%
  mutate(feature= paste(measure, method))%>% dplyr::rename(name= net)

dis.net.feats= lapply(names(edge.list$disease), function(x){
  getProps_global(edge.list$disease[[x]], directed = F)%>% mutate(net= x)
})
df.d= do.call(rbind, dis.net.feats)%>% mutate(value = as.numeric(value))%>%
  mutate(feature= paste(measure, method))%>% dplyr::rename(name= net)



?assortativity_degree
#
#
# ggplot(df, aes(x= net, y= value))+
#   facet_grid(rows= vars(method), scales = "free_y")+
#   geom_col()
#
#
# edge_density(edge.list$gene$goMF, loops = F
# )
#
# df.size
# df$name

# combine and plot --------------------------------------------------------

corr.df= corr %>%
  mutate(feature= "pubmed_correlation") %>%
  dplyr::rename(name= layer, value= corr)

df= rbind(df.g[,c("feature", "name", "value")],
           df.d[,c("feature", "name", "value")],
          df.size[,c("feature", "name", "value")],
          corr.df[,c("feature", "name", "value")])

unique(df$name)
unique(df$feature)
df= df %>%
  mutate("layer"= ifelse(name %in% c("heidelberg", "hpo"), "disease", "gene"),
         name =ifelse(name== "heidelberg", "HFnet",
                      ifelse(name=="hpo", "HPOnet", name))
  )

# plot net sizes
p1= ggplot(df%>% filter(feature %in% c("# nodes" )), aes(x= name, y= value))+
  facet_grid(cols= vars(layer), scales = "free")+
  geom_col(width= 0, colour= "black", lwd= 1)+
  geom_point(size= 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  labs(x= "", y= "# nodes")

p2= ggplot(df%>% filter(feature %in% c("# edges" )), aes(x= name, y= value))+
  facet_grid(cols= vars(layer), scales = "free")+
  geom_col(width= 0, colour= "black", lwd= 1)+
  geom_point(size= 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  labs(x= "", y= "# edges")
p3= ggplot(df%>% filter(feature %in% c("edge_density percent" )), aes(x= name, y= value))+
  facet_grid(cols= vars(layer), scales = "free")+
  geom_col(width= 0, colour= "black", lwd= 1)+
  geom_point(size= 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  labs(x= "", y= "edge_density")

p4= ggplot(df%>% filter(feature %in% c("centrality degree" )), aes(x= name, y= value))+
  facet_grid(cols= vars(layer), scales = "free")+
  geom_col(width= 0, colour= "black", lwd= 1)+
  geom_point(size= 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  labs(x= "", y= "centrality")

p5= ggplot(df%>% filter(feature %in% c("transitivity global" )), aes(x= name, y= value))+
  facet_grid(cols= vars(layer), scales = "free")+
  geom_col(width= 0, colour= "black", lwd= 1)+
  geom_point(size= 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  labs(x= "", y= "transitivity")


p6= ggplot(df%>% filter(feature %in% c("pubmed_correlation" )), aes(x= name, y= value))+
  facet_grid(cols= vars(layer), scales = "free")+
  geom_col(width= 0, colour= "black", lwd= 1)+
  geom_point(size= 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  labs(x= "", y= "pubmed_corr")

p6= cowplot::plot_grid(NULL, p6, ncol = 2)

p7= ggplot(df%>% filter(feature %in% c("degree_assortativity global" )), aes(x= name, y= value))+
  facet_grid(cols= vars(layer), scales = "free")+
  geom_col(width= 0, colour= "black", lwd= 1)+
  geom_point(size= 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  labs(x= "", y= "degree_assort")
cowplot::plot_grid(p1,p2,p3,p4,p5 ,p7,p6, ncol = 2)

df2= df%>% dplyr::filter(feature %in% c("# edges","# nodes","edge_density percent", "centrality degree","transitivity global", "pubmed_correlation"))

ggplot(df2, aes(x= name, y= value))+
  facet_grid(cols= vars(feature), scales = "free")+
  geom_col(width= 0, colour= "black", lwd= 1)+
  geom_point(size= 2)+
  #theme_minimal_hgrid()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  labs(x= "")

ggplot(df%>% filter(layer %in% c("disease" )), aes(x= name, y= value))+
  facet_grid(rows= vars(feature), scales = "free_y")+
  geom_col(width= 0, colour= "black", lwd= 1)
  geom_point(size= 2)+
  #theme_bw()+
  theme(axis.text.x = element_text(angle= 40, hjust= 1))+
  labs(x= "")

