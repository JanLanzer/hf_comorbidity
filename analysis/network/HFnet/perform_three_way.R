## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-10-11
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##  Performing a three way test between hfref and hfpef
##
## ---------------------------
library(igraph)
library(tidyverse)
library(reshape2)
library(cowplot)
library(RANKS)
library(RandomWalkRestartMH)

source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_network.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_pairwise_disease.R")
source("~/GitHub/hf_comorbidity_genes/analysis/utils/utils_hetnet.R")

data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/ICD10_labeled_phe2022.rds")

# PheCode dictionary:
Phe_dic= data %>%
  distinct(PheCode, Phenotype,category) %>%
  drop_na

#patients of interest

# define both cohorts -----------------------------------------------------

hf.pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")
phecodes= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/topfreq_disease.rds")


# Calculate disease networks ----------------------------------------------

dd_links= create_edge_table_for_cohorts(data = data,
                                        pid.list = hf.pids.list,
                                        topn_disease = 0,
                                        phecodes = phecodes )

links.hfpef= create.links(pids= hf.pids.list$hfpef,
                     phecodes = phecodes,
                     data= data)
links.hfref= create.links(pids= hf.pids.list$hfref,
                          phecodes = phecodes,
                          data= data)
links.hfref%>% filter(!is.finite(odds.ratio))




saveRDS(dd_links,"T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/links_threeway.rds" )
dd_links= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/links_threeway.rds" )
dd_hfnet= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_table_hf_cohort.rds" )

# Functions  -----------------------------------


# perform Breslow Day test:

# takes to link tables with count data filling cont. table.
# merges both and performs breslow day test for difference of both, returns merged data frame
# with added p.val and adjust p.val
links1= dd_links$links$hfpef
links2= dd_links$links$hfref

add.breslow.test= function(links1, links2 , name1= "links1", name2= "links2", ...){

  require(DescTools)

  ref= get_edges_with_edgeID(links= links1)
  pef= get_edges_with_edgeID(links=links2)

  comb.= ref %>% inner_join(pef, by= c("edge_IDs"))

  links.comb= comb. %>%
    rowwise()%>%
    mutate(p.breslow= DescTools::BreslowDayTest(array(
      c(a.x+1,
        c.x+1,
        b.x+1,
        d.x+1,
        a.y+1,
        c.y+1,
        b.y+1,
        d.y+1),
      dim=c(2, 2, 2)
    )
    )$p.value) %>%
    ungroup()%>%
    mutate(p.bres.adj = p.adjust(p.breslow, method= "BH")) %>%
    mutate(hf_link= ifelse(abs(log(odds.ratio.x))>abs(log(odds.ratio.y)), name1, name2))%>%
    dplyr::select(-disease1.y, -disease2.y, -dis1_phenotype.y, -dis2_phenotype.y)

  return(links.comb)
}

#set infinite OR to max OR
reset.OR= function(links){
  max.OR= max(links$odds.ratio[is.finite(links$odds.ratio)], na.rm = T)
  links = links %>%
    mutate(odds.ratio = ifelse(is.infinite(odds.ratio), max.OR, odds.ratio))

}


#wrapper, takes same
get_breslow_disease_links= function(links1, links2,name1= "links1", name2= "links2",
                                    alpha.breslo= 0.05,
                                    alpha.fisher= 0.05,
                                    ...){

  #before calculating breslow p, we will set infinite OR to max OR per cohort

  links1= reset.OR(links1)
  links2= reset.OR(links2)

  sig.links = add.breslow.test(links1 , links2, name1, name2)

  sig.links= sig.links  %>%
    mutate(sig= ifelse(p.breslow<alpha.breslo, "s", "ns"),
           sig.x= ifelse(fisher.p.adj.x<alpha.fisher, "s", "ns"),
           sig.y= ifelse(fisher.p.adj.y<alpha.fisher, "s","ns"),
           sig.comb= ifelse(fisher.p.adj.x<alpha.fisher & fisher.p.adj.y<alpha.fisher,
                        paste0("sig", name1, name2),
                        ifelse(fisher.p.adj.x<alpha.fisher,
                               paste0("sig", name1),
                               ifelse(fisher.p.adj.y<alpha.fisher,
                                      paste0("sig", name2),
                                      "ns")
                               )
                        )
           )


  colnames(sig.links) = str_replace_all(colnames(sig.links), "\\.x", name1)
  colnames(sig.links) = str_replace_all(colnames(sig.links), "\\.y", name2)
  return(sig.links)

  # # filter for links that are significiantly different between both cohorts
  # # but also have 2-way significance
  # links.x= sig.links %>% filter(sig.x== "sig.x" & sig== "s") %>% pull(edge_IDs)
  # links.y= sig.links %>% filter(sig.x== "sig.y" & sig== "s") %>% pull(edge_IDs)
  #
  # ## use those edge IDs to filter each link table:
  # netx= get_edges_with_edgeID(links= links1)
  # nety= get_edges_with_edgeID(links= links2)
  #
  # netx= netx %>% filter(edge_IDs %in% links.x)
  # nety = nety %>% filter(edge_IDs %in% links.y)
  #
  # net.t1= table_to_links(netx, p.val= 1, weight_dd = 0.05)
  # # net.t2= table_to_links(nety, p.val= 1,  weight_dd = 0.05)
  #
  # return(list("merged"= sig.links,
  #             "separate"= list("x" = netx,
  #                              "y"= nety
  #                              )
  #             )
  #        )
}


#wrapper for making plots
visualize_test_results= function(sig.links){

  # plot all links and show which disease links are shared based on fisher p-adj =

  p1= ggplot(sig.links, aes(x= log10(odds.ratio.hfpef), y= log10(odds.ratio.hfref), col= sig.comb))+
    geom_point(alpha= 0.8)+
    scale_color_manual(values = cols.nice)

  #plot only significant links from p.breslow


  df= sig.links %>% filter(log10(odds.ratio.hfref)>4)
  p2= sig.links %>%
    filter(p.bres.adj< 0.05)%>%
    ggplot(. , aes(x= log10(odds.ratio.hfpef), y= log10(odds.ratio.hfref), col= -log10(p.bres.adj)))+
     geom_point(alpha= 0.8)+
    scale_color_gradient(low= "red", high = "black")
    # scale_x_continuous(limits = c(-2, 10))+
    # scale_y_continuous(limits = c(-2, 10))

  test.df= sig.links %>% filter(p.breslow< 0.05)
  # of all links, how many are significiant in either one?
  tab1= prop.table(table(test.df$fisher.p.adj.hfref<0.05  | test.df$fisher.p.adj.hfpef<0.05  ))

  # how many are significinat in both?
  tab2= prop.table(table(test.df$fisher.p.adj.hfref<0.05  & test.df$fisher.p.adj.hfpef<0.05 ))

  return(list("plot.all"= p1,
              "plot.sig"= p2,
              "tab.1"= tab1,
              "tab.2"= tab2)
  )

}



# Run breslow dayes test for independent links from hfref and hfpef net ------------------------------------------------

hfpef.hfref= get_breslow_disease_links(dd_links$links$hfref,
                          dd_links$links$hfpef,
                          name1= ".hfref",
                          name2= ".hfpef"

                          )

hfpef.hfref= get_breslow_disease_links(links.hfpef,
                                       links.hfref,
                                       name1= ".hfref",
                                       name2= ".hfpef",
                                       alpha.fisher = 0.0001
                                       )

saveRDS(hfpef.hfref, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/three_way_obj.rds")
hfpef.hfref= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/three_way_obj.rds")


#correlated:
cor.test(hfpef.hfref$odds.ratio.hfpef[is.finite(hfpef.hfref$odds.ratio.hfpef)],
         hfpef.hfref$odds.ratio.hfref[is.finite(hfpef.hfref$odds.ratio.hfref)])

prop.table(table(hfpef.hfref$sig.comb))

p.= visualize_test_results(sig.links = hfpef.hfref)
p.$plot.all


p.1= hfpef.hfref %>%
  filter(odds.ratio.hfpef< 10000000)%>% # remove one outlier for plotting
  ggplot(., aes(x= log10(odds.ratio.hfpef), y= log10(odds.ratio.hfref), col= sig.comb))+
  geom_point(alpha= 0.8)+
  scale_color_manual(values = cols.nice)+
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(col = "Significant\nin cohort:")



pcomb.= plot_grid(p.$plot.all+
            theme_minimal(),
          p.$plot.sig+
            theme_minimal(),
          rel_widths = c(1,0.8),
          labels = "AUTO")



df.breslow2= hfpef.hfref %>% select(edge_IDs, dis1_phenotype.hfref, dis2_phenotype.hfref,
                                   odds.ratio.hfpef,odds.ratio.hfref,
                                   fisher.p.adj.hfpef, fisher.p.adj.hfref,
                                   p.breslow, p.bres.adj)

max(hfpef.hfref$odds.ratio.hfpef)
log10(max(hfpef.hfref$odds.ratio.hfref))

df.breslow2 %>% mutate(loghfpef= log10(odds.ratio.hfpef),
                       loghfref= log10(odds.ratio.hfref))%>%
  filter(loghfpef> 3.9
         )

# hist(df.breslow2$p.breslow, bins= 100)
# hist(df.breslow2$p.bres.adj, bins= 100)
df.breslow2%>% filter(p.bres.adj<0.01)%>% mutate(loghfpef= log10(odds.ratio.hfpef),
                                                 loghfref= log10(odds.ratio.hfref))
p.heat= df.breslow2 %>%
  arrange(p.bres.adj) %>%
  mutate(log.p= -log10(p.breslow+0.0000001))%>%
  filter(p.bres.adj<0.01) %>%
  mutate(label = paste0(dis1_phenotype.hfref, "+",  dis2_phenotype.hfref))%>%
  dplyr::rename(hfpef= odds.ratio.hfpef,
                hfref= odds.ratio.hfref)%>%
  pivot_longer(cols = c(hfpef, hfref), names_to = "cohort", values_to = "odds.ratio")%>%
  ggplot(., aes(x= cohort, y=reorder(label,-p.bres.adj), fill = log10(odds.ratio)))+
  geom_tile(color ="black")+
  scale_fill_gradient2(low= "blue", mid= "white" ,high = "red")+
  theme_minimal()+
  theme(        legend.position = "left")+
  labs(y= "")

p.heat
p.hist= df.breslow2 %>%
  arrange(p.bres.adj) %>%
  mutate(log.p= -log10(p.bres.adj+0.0000001))%>%
  filter(p.bres.adj<0.01) %>%
  ggplot(., aes(x= reorder(edge_IDs,-p.bres.adj),  y= log.p))+
  geom_col(color= "black")+
  geom_hline(yintercept = -log10(0.05))+
  coord_flip()+
  theme_minimal()+
  theme(axis.text.y = element_blank(),
       axis.title.y = element_blank())+
  labs(y= "-log10 p-value")

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/breslow.hmap.pdf",
    width= 10,
    height= 10)
plot_grid(unify_axis(p.heat),
          unify_axis(p.hist),
          rel_widths = c(1, 0.1))
dev.off()

pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/pairwise_disease/breslow.OR.comp.pdf",
    width= 8,
    height= 4)
unify_axis(p.1)+coord_equal()
dev.off()



#save all
df.breslow2 %>%
  write.csv2("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/manuscript/diseasenet/link.table.breslow.hfpef.hfref.csv")
