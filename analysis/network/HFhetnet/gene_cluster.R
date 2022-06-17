## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-05-27
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## multilayer cluster
## ---------------------------


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




