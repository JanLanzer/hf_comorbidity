## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-11-08
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## compare PPIs
## ---------------------------

ppi= process_human_ppi_multiverse()
ppi2= process_ppi_v2()
ppi3= process_ppi_lit()

net.list= lapply(list(ppi, ppi2, ppi3), function(x){
  graph_from_data_frame(x, directed= F)
})

v.list= map(net.list, function(y){
  unlist(V(y)$name)
})

vs = intersect(v.list[[3]],intersect(v.list[[2]], v.list[[1]]))


net.list= lapply(net.list, function(x){
  induced_subgraph(x, vids = V(x)$name %in% vs)
})


df= sapply(net.list, function(x){
  sapply(net.list, function(y){
    deltacon_igraph(x,y)
  })
})

rbind(ppi, ppi3)
