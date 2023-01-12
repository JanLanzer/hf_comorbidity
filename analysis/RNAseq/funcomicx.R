## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-11-21
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## funcomicx
library(tidyverse)
library(ComplexHeatmap)
library(decoupleR)
library(progeny)
library(dorothea)

data("dorothea_mm")

dea= readRDS( "output/dea.results.list.rds")

get_matrix= function(c.list){

  y = map(c.list, function(x){
    x$t
  })

  M= do.call(cbind, y)
  rownames(M)= c.list[[1]]$gene
  return(M)
}

M.T= get_matrix(dea)

regulons= dorothea_mm %>%
  filter(confidence %in% c("A", "B", "C", "D")) %>%
  distinct(tf, target, mor) %>%
  mutate(likelihood=1)

dec.res = decouple(M.T,
                   network = regulons,

                   .source ="tf",
                   .target = "target",
                   statistics = "wmean",

                   consensus_score = F

)

unique(dec.res$statistic)

dec.res %>% filter(condition == "N",
                   statistic== "corr_wmean",
                   grepl("At", source))
dec.res%>% filter(condition== "TAC.dko.vs.wt",
                  statistic== "corr_wmean")%>%
  arrange(desc(score))

tf.oi= dec.res%>%
  group_by(condition)%>%
  filter(any(p_value<0.05),
         statistic== "norm_wmean")%>%
  slice_max(order_by = score, n=10)%>%
  pull(source)

dec.res%>%
  group_by(condition)%>%
  filter(source %in% tf.oi,
         statistic== "norm_wmean")%>%
  ggplot(., aes(x= condition, y= source, fill = score))+
  geom_tile()+
  scale_fill_gradient2(low= "blue",mid = "white",  high= "red")


df= dec.res %>%
  group_by(condition) %>%
  filter(source %in% tf.oi,
         statistic== "norm_wmean") %>%
  dplyr::select(condition, source, score) %>%
  pivot_wider(names_from = "condition",
              values_from = "score")

p1= df %>% column_to_rownames("source")%>%
  Heatmap()

pdf("output/figures/tf_act.pdf",
    width= 3.5,
    height= 11
)
p1
dev.off()


x= progeny(M.T, organism = "Mouse")

p2= Heatmap(t(x))

pdf("output/figures/pathway_act.pdf",
    width= 3.5,
    height= 5
)
p2
dev.off()


# progeny -------------------------------------------------------------------------------------

x= progeny(expr = M.T, organism = "Mouse", top = 200)
p2= Heatmap(t(x))

PathwayActivity_zscore <- progeny(M.T,
                                  scale=TRUE, organism="Mouse", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscore) <- "NES"

Heatmap(PathwayActivity_zscore[,1:2], name = "PROGNEy z-score")

gex.obj$meta$sampleID== colnames(t(x))
