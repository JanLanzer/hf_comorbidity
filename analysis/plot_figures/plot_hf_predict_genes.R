## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-11-14
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## plot hf predictions
## ---------------------------
library(tidyverse)
library(ggrepel)
source("analysis/utils/utils_hetnet.R")

sets= load_validation_genes(disgenet_value= 0.29)
names(sets)= c("DisGeNET", "PheWAS", "Kegg_DCM", "ReHeaT", "r2", "Cardiomyopathy", "Top_single_variants", "Top_common_variants")
sets= sets[-5]
val.set = unlist(sets)

source("analysis/utils/utils.R")
g.ranks= readRDS( file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HF_gene_ranks.rds")

top_n = 25

g.pef.top= g.ranks %>% arrange(desc(hfpef.prio))%>% dplyr::slice(1:top_n) %>% pull(gene)
g.ref.top= g.ranks %>% arrange(desc(hfref.prio))%>% dplyr::slice(1:top_n) %>% pull(gene)

p.full.ranks= g.ranks %>% ggplot(., aes(x= rank.hfpef, y= rank.hfref))+
  geom_point()


p.pef.prio= g.ranks%>% mutate(label = ifelse(gene %in% g.pef.top, gene, ""),
                  label2= ifelse(RW.value.hfpef >1e-04, gene, ""),
                  col= ifelse(rank.diff >500, "yes", "no"))%>%
  ggplot(., aes(x= RW.value.hfpef, y= abs(rank.diff), col = hfpef.prio))+
  geom_point()+
  scale_color_gradient2(low= "black",mid= "#b30000",  high = "#ff0000", midpoint= 0.09)+
  scale_y_log10()+
  geom_text_repel(aes(label= label ), box.padding = 2,segment.alpha= 0.3,
                  max.overlaps = 130,
                  show.legend = FALSE)+
  theme_classic()+
  labs(x= "RW probability HFpEF",
       y= "gene ranking difference",
       col ="HFpEF\nprioritization")+
  theme(axis.text = element_text(size= 14))+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))
unify_axis(p.pef.prio)


p.ref.prio= g.ranks%>% mutate(label = ifelse(gene %in% g.ref.top, gene, ""),
                              label2= ifelse(RW.value.hfref >1e-04, gene, ""),
                              col= ifelse(rank.diff <-500, "yes", "no"))%>%
  ggplot(., aes(x= RW.value.hfref, y= abs(rank.diff), col = hfref.prio))+
  geom_point()+
  scale_color_gradient2(low= "black",mid= "#b30000",  high = "#ff0000", midpoint= 0.35)+
  scale_y_log10()+
  geom_text_repel(aes(label= label ), box.padding = 2,segment.alpha= 0.3,
                  max.overlaps = 1000,
                  show.legend = FALSE)+
  theme_classic()+
  labs(x= "RW probability HFrEF",
       y= "gene ranking difference",
       col ="HFpEF\nprioritization")+
  theme(axis.text = element_text(size= 14))+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))
unify_axis(p.ref.prio)

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/HFpEF.gene.rankings.score2.pdf",
    height= 5,
    width= 9)

unify_axis(p.ref.prio)
unify_axis(p.pef.prio)


dev.off()



## PLOT RW PROBs
p.ef = ggplot(g.pef, aes(x= rank,y= value))+
  geom_point()+
  theme_classic()+
  geom_vline(xintercept = 500, col= "darkgrey")+
  labs(x= "gene ranking hfpef",
       y= "RW probability")+
  ggtitle("HFpEF")

p.ef

p.ref = ggplot(g.ref, aes(x= rank,y= value))+
  geom_point()+
  theme_classic()+
  geom_vline(xintercept = 500, col= "darkgrey")+
  labs(x= "gene ranking hfref",
       y= "RW probability")+
  ggtitle("HFrEF")

p.ref

p.cutoff= plot_grid(unify_axis(p.ef),
                    unify_axis(p.ref), labels = "AUTO")

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/hetnet/gene.rankings.RW.prob.pdf",
    height= 5,
    width= 8)
p.cutoff
dev.off()

############ plot pef and ref gene candidates

p.pef.top250 =
  g.ranks %>%
  filter(rank.hfpef<500)%>%
  mutate(label = ifelse(gene %in% g.pef.top, gene, ""),
         col= ifelse(rank.diff >500, "yes", "no")) %>%
  ggplot(., aes(x= rank.hfpef,
                y= rank.hfref,
                col= hfpef.prio))+
  geom_point()+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  #geom_hline( yintercept = log10(100), type = 2, color= "grey")+
  labs(x= "gene ranking HFpEF",
       y= "gene ranking HFrEF",
       col ="RW probability * \n ranking difference")+
  #geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()+
  scale_color_gradient(low= "darkgrey", high= "blue")+
  theme(axis.text = element_text(size= 14))+
  ggtitle("HFpEF Ranking top 250")

p.pef.top250

p.pef.top100 = g.ranks %>%
  filter(rank.hfpef<200)%>%
  mutate(label = ifelse(gene %in% g.pef.top, gene, ""),
         col= ifelse(rank.diff >500, "yes", "no")) %>%
  ggplot(., aes(x= rank.hfpef,
                y= rank.hfref,
                col= col))+
  geom_point()+
  scale_color_manual(values= cols.nice[c(5,1)])+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  #geom_hline( yintercept = log10(100), type = 2, color= "grey")+
  labs(x= "gene ranking HFpEF",
       y= "gene ranking HFrEF",
       col ="ranking\ndifference \n >500")+
  geom_text_repel(aes(label= label ), box.padding = 2,segment.alpha= 0.3,max.overlaps = 50,show.legend = FALSE)+
  theme_classic()+
  theme(axis.text = element_text(size= 14))


p.pef.top100

p.pef.HFgenes = g.ranks %>%
  filter(rank.hfpef<500,
         rank.hfref<500)%>%
  mutate(label = ifelse(gene %in% val.set, gene, ""),
         label2 = ifelse(!gene %in% val.set, "other", "hf"),
         col= ifelse(rank.diff >500, "yes", "no")) %>%
  ggplot(. ,aes(x= rank.hfpef,
                y= rank.hfref,
                col= label2))+
  geom_point( size= 2)+
  scale_color_manual(values= cols.nice[c(1,4)])+
  geom_abline(slope =1, intercept = 0)+
  #geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  labs(x= "gene ranking HFpEF",
       y= "gene ranking HFrEF",
       col ="HF gene set")+
  geom_text_repel(aes(label= label ),segment.alpha= 0.3,  box.padding = 2, max.overlaps = 40,show.legend = FALSE)+
  theme_classic()+
  theme(panel.border  = element_rect(size= 1, fill  = NA))

unify_axis(p.pef.HFgenes)

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/HFpEF.gene.rankings.score.pdf",
    height= 5,
    width= 5)
unify_axis(p.pef.top250)+theme( panel.border = element_rect(colour = "black", fill=NA, size=1))

dev.off()

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/HFpEF.gene.rankings.dg2.pdf",
    height= 4.5,
    width= 9)

unify_axis(p.pef.top100)+theme( panel.border = element_rect(colour = "black", fill=NA, size=1))
unify_axis(p.pef.HFgenes)+theme( panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()


########## Plot same for ref

p.ref.top250 = g.ranks %>%
  filter(rank.hfref<250)%>%
  mutate(label = ifelse(gene %in% g.ref.top, gene, ""),
         col= ifelse(rank.diff < (-500), "yes", "no")) %>%
  ggplot(., aes(x= rank.hfref,
                y= rank.hfpef,
                col= col))+
  geom_point()+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  #geom_hline( yintercept = log10(100), type = 2, color= "grey")+
  labs(x= "ranking hfref",
       y= "ranking hfpef",
       col ="ranking difference \n >500")+
  geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()+
  theme(axis.text = element_text(size= 14))+
  ggtitle("HFrEF Ranking top 250")

p.ref.top100 = g.ranks %>%
  filter(rank.hfref<100)%>%
  mutate(label = ifelse(gene %in% g.ref.top, gene, ""),
         col= ifelse(rank.diff < (-500), "yes", "no")) %>%
  ggplot(., aes(x= rank.hfref,
                y= rank.hfpef,
                col= col))+
  geom_point()+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  #geom_hline( yintercept = log10(100), type = 2, color= "grey")+
  labs(x= "ranking hfref",
       y= "ranking hfpef",
       col ="ranking difference \n >500")+
  geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()+
  theme(axis.text = element_text(size= 14))+
  ggtitle("HFrEF Ranking top 100")

p.ref.top250
p.ref.top100

p.ref.HFgenes = g.ranks %>%
  filter(rank.hfref<250)%>%
  mutate(label = ifelse(gene %in% val.set, gene, ""),
         label2 = ifelse(!gene %in% val.set, "other", "hf"),
         col= ifelse(rank.diff >(-500), "yes", "no")) %>%
  ggplot(. ,aes(x= rank.hfref,
                y= rank.hfpef,
                col= label2))+
  geom_point( size= 2)+
  scale_color_manual(values= cols.nice[c(1,4)])+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  labs(x= "ranking hfref",
       y= "ranking hfpef",
       col ="HF gene set")+
  geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()

p.ref.HFgenes

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/figures/manuscript/main/HFrEF.gene.rankings.dg2.pdf",
    height= 8,
    width= 12)

p.ref.top100
p.ref.top250
p.ref.HFgenes

dev.off()


p.ref.HFgenes2 = g.ranks %>%
  filter(rank.hfref<250)%>%
  mutate(label = ifelse(abs(rank.diff)<50, name, ""),
         label2 = ifelse(!name %in% val.set, "other", "hf"),
         col= ifelse(rank.diff >(-500), "yes", "no")) %>%
  ggplot(. ,aes(x= rank.hfref,
                y= rank.hfpef,
                col= label2))+
  geom_point( size= 2)+
  scale_color_manual(values= cols.nice[c(1,4)])+
  geom_abline(slope =1, intercept = 0)+
  geom_abline(slope =1, intercept = 500, col= "darkgrey")+
  labs(x= "ranking hfref",
       y= "ranking hfpef",
       col ="HF gene set")+
  geom_text_repel(aes(label= label ), box.padding = 2, max.overlaps = 50,show.legend = FALSE)+
  theme_classic()

p.ref.HFgenes2

# compare network modules -------------------------------------------------


monet= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/gene.modules.MONET.M1.rds")

names(monet)= paste0("c.", seq(1:length(monet)))

monet_length= lapply(monet, length)

names.5= names(monet_length[unlist(map(monet_length,function(x)(x>5)))])

enframe(monet, name= "cluster", value= "name")%>%
  unnest(name)

after_clustering =
  g.ranks %>%
  left_join(enframe(monet, name= "cluster", value= "name")%>%
              unnest(name))%>% print(n=100)

after_clustering%>%
  arrange(desc(hfpef.prio))

t= after_clustering%>% group_by(cluster)%>% mutate(mean_diff= median(hfpef.prio))%>%
  filter(cluster %in% names.5)

t%>% arrange(desc(mean_diff)) %>% distinct(cluster, mean_diff) %>% print(n=100)


after_clustering %>% filter(cluster==  "c.373")%>% print(n=100)


t%>% distinct(cluster, mean_diff) %>% arrange(desc(mean_diff))%>% print(n=500
)
after_clustering%>% filter(grepl("ANGP", name))
g.ranks %>%
  arrange(desc(hfpef.prio)) %>% slice(1:50)




