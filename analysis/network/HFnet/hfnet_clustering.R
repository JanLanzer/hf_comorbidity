## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2021-12-09
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2021
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##  compare different cluster algorithms in igraph for the disease network
##
## ---------------------------

library(WriteXLS)
library(tidyverse)
library(igraph)
library(pheatmap)
library(ComplexHeatmap)
library(cowplot)

net= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")

set.seed(11)
######## run different cluster algos:

info = igraph::infomap.community(net,
                  e.weights = E(net)$weight,
                  v.weights = NULL,
                  nb.trials = 500,
                  modularity = T)

louv= cluster_louvain(net,
                      weights = E(net)$weight,resolution = 1.0

                        )

leid= cluster_leiden(net,objective_function = "modularity",
                                  weight= E(net)$weight,
                                  resolution_parameter = 1.0)
fast= cluster_fast_greedy(net,
                          weights = E(net)$weight)

walktrap = cluster_walktrap(net,
                            weights = E(net)$weight)

# label= cluster_label_prop(graph = net,
 #                          #initial= V(net)$group_cat,
  #                         weights =NULL)# E(net)$weight)

eigen= cluster_leading_eigen(net,
                             weight= E(net)$weight)

spinglass= cluster_spinglass(net,
                             weight= E(net)$weight)


# Plot and compare --------------------------------------------------------------------

modules= list(#"infomap"= info,
              "louvain"= louv,
              "fast_greedy"= fast,
              "walktrap" = walktrap,
              #"label_prop"= label,
              "eigen"= eigen,
              "spinglass"= spinglass,
              "leiden"= leid)


A= matrix(NA,
       nrow=length(modules),
       ncol= length(modules))

colnames(A) = rownames(A)= names= names(modules)

# adjusted_rand
A.adj= A
for (i in names ){
  for(j in names){
    A.adj[i,j] = compare(modules[[i]], modules[[j]], method = "adjusted.rand")
  }
}

A.nmi= A
for (i in names ){
  for(j in names){
    A.nmi[i,j] = compare(modules[[i]], modules[[j]], method = "nmi")
  }
}



p.compare = pheatmap(A.nmi)#, main = "normalized mutual information between cluster algorithms for HF-net")
Heatmap((A.adj))

#use order based on hierarchical clustering:
p.compare$tree_col$order

ordered= c("spinglass", "louvain","leiden", "walktrap", "fast_greedy", "eigen", "label_prop")#, "infomap")
ordered=  p.compare$tree_col$labels[p.compare$tree_col$order]

heatmap.df.nmi= as.data.frame(A.nmi) %>% rownames_to_column("algorithm1") %>%
  pivot_longer(cols = !algorithm1,
               names_to = "algorithm2",
               values_to= "nmi") %>%
  mutate(algorithm1 = factor(algorithm1, levels = ordered)) %>%
  mutate(algorithm2 = factor(algorithm2, levels = ordered))

p.nmi = ggplot(heatmap.df.nmi, aes(x= algorithm1, y= algorithm2, fill = nmi))+
  geom_tile()+
  scale_fill_gradient2(low = "white", high = "darkred")+
  theme_minimal()+
  labs(fill = "NMI")

# calc row max

dfx= rbind(enframe(apply(A.nmi,1, mean))%>% mutate(feat= "NMI"),
  enframe(apply(A.adj,1, mean))%>% mutate(feat= "adj.RI")
  )

p.mean.comps= ggplot(dfx, aes(x= name, y= value, fill = feat))+
  geom_col()+
  facet_grid(rows= vars(feat))+
  theme(axis.text.x = element_text(angle= 45, hjust= 1))+
  labs(fill = "")

p.mean.comps=  unify_axis(p.mean.comps+
    theme(axis.title = element_blank())
 )

barplot(apply(A.nmi,1, mean))

barplot(apply(A.adj,1, mean))

heatmap.df.adj= as.data.frame(A.adj) %>% rownames_to_column("algorithm1") %>%
  pivot_longer(cols = !algorithm1,
               names_to = "algorithm2",
               values_to= "adjusted.rand") %>%
  mutate(algorithm1 = factor(algorithm1, levels = ordered)) %>%
  mutate(algorithm2 = factor(algorithm2, levels = ordered))

p.adj = ggplot(heatmap.df.adj, aes(x= algorithm1, y= algorithm2, fill = adjusted.rand))+
  geom_tile()+
  scale_fill_gradient2(low = "white", high = "darkblue")+
  theme_minimal()+
  labs(fill = "adjusted\nRand Index")


# compare for cluster sizes
df = do.call(rbind,
        lapply(names(modules), function(x){
            df= enframe(sizes(modules[[x]]))%>% mutate(algo = x)
        })
        )
df%>%
  left_join(df %>% count(algo), by= "algo") %>%
  mutate(algo = factor(algo, levels= ordered))


p.clustersize= ggplot(df, aes(x=algo, y= value))+
  geom_point()+
  labs(y= "nodes per cluster",
       x= "community detection algorithm")+
  theme_minimal()

p.size = ggplot(df, aes(x=algo))+
  geom_bar()+
  labs(y= "number of clusters",
       x= "community detection algorithm")+
  theme_minimal()

# compare for modularity:

modularity= sapply(modules, function(x){modularity(net , membership(x), weights= E(net)$weight)})
modularity_df= as_tibble(modularity) %>%
  mutate(algo =names(modularity)) %>%
  mutate(modul= "modularity")%>%
  mutate(algo = factor(algo, levels= ordered))

p.mod= ggplot(modularity_df, aes(x= algo,y= modularity))+
  geom_col()+
  theme_minimal()


p1= plot_grid(p.size+
                labs(x="")+
                theme(axis.text.x = element_text(angle = 70, hjust =1)),
              p.clustersize+
                labs(x="")+
                theme(axis.text.x = element_text(angle = 70, hjust =1)),
              p.mod+
                labs(x="")+
                theme(axis.text.x = element_text(angle = 70, hjust =1))
              , nrow = 1,
              rel_heights = c(0.5, 1,0.5) ,labels= "AUTO"
              )

p.cluster= plot_grid( unify_axis(p.nmi+
                        labs(x="", y= "")+
                        theme(axis.text.x = element_text(angle = 70, hjust =1))),
                      unify_axis(p.adj+
                        labs(x="", y="")+
                        theme(axis.text.x = element_text(angle = 70, hjust =1))),
                        nrow= 1,
                      labels= c("D", "E"))

p.comb= plot_grid(p1, (p.cluster), nrow= 2, rel_heights =c(0.8, 1) )
p.comb

pdf(file= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/cluster.algo.compare.pdf",
    width = 8, height= 6)
p.comb
dev.off()

##check resolution for leiden

# resolution of leiden ----------------------------------------------------

test.seq= seq(0.1,2, 0.1)

df= sapply(test.seq, function(x){

  leid= cluster_leiden(net,objective_function = "modularity",
                       weight= E(net)$weight,
                       resolution_parameter = x)
  c1= modularity(net , membership(leid), weights= E(net)$weight)
  c2= sapply(modules, function(x){
   compare(x, leid, method = "nmi")

  })%>%
    median()
  #leid$nb_clusters
  return(c(c1,c2,leid$nb_clusters))

})

colnames(df)= test.seq
rownames(df)= c("modularity", "NMI", "nclust")

p.res= df %>% as.data.frame() %>% rownames_to_column("feat")%>% as_tibble()%>%
  pivot_longer(-feat)%>%
  ggplot(., aes(x= name, y= value, col = feat))+
  geom_point(size= 3)+
  facet_grid(rows= vars(feat), scales= "free")+
  labs(x= "leiden_resolution",
       col= "feature")+
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle= 45, hjust= 1))
p.res

pdf(file= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/leide.res.pdf",
    width = 6, height= 5)
unify_axis(p.res)

dev.off()




# HF net cluster ----------------------------------------------------------
## we will use leiden clustering at resolution 1.1

set.seed(2)

leidenpart= cluster_leiden(net,objective_function = "modularity",
                     weight= E(net)$weight,
                     resolution_parameter = 1.1)

leidenpart$nb_clusters

table((leidenpart$membership))

## add cluster to node features
V(net)$group_louv= leidenpart$membership

nodes = as_tibble(igraph::as_data_frame(net, what = "vertices"))

nodes_count= nodes %>%
  complete(group_louv, group_cat)%>%
  mutate(group_louv = factor(group_louv))%>%
  group_by(group_louv, group_cat) %>%
  count()


mat = nodes_count %>%
  ungroup()%>%
  group_by(group_louv)%>%
  mutate( prop = (n/sum(n))*100)

cat.count= mat%>%
  ungroup()%>%
  group_by(group_cat)%>%
  mutate(x= sum(n))%>%
  distinct(group_cat, x)

library(circlize)
col_fun = colorRamp2(c(0,50, 100), c("white", "red", "darkred"))
col_fun(seq(0, 1))

cluster.size= mat%>% group_by(group_louv) %>% summarise(sum= sum(n))

ha = HeatmapAnnotation("cluster\nsize" = anno_barplot(cluster.size$sum,
                                                      add_numbers = TRUE,
                                                      height = unit(1.2, "cm"),
                                                      border= F))

ra = rowAnnotation("    category\nsize" = anno_barplot(cat.count$x,
                                                   add_numbers = TRUE,border = F,
                                                   #show_annotation_name= F,
                                                   height = unit(2, "cm")))

h.map= mat %>%
  select(-n)%>%
  mutate(group_louv= paste("DC.",group_louv))%>%
  pivot_wider(., names_from = group_louv, values_from= prop)%>%
  column_to_rownames( "group_cat")%>%
  Heatmap(., cluster_rows = F,
          cluster_columns = F,
          col= col_fun,
          top_annotation = ha,
          right_annotation = ra,
          show_heatmap_legend = T,
          name= "%",
          row_names_side = "left",
          rect_gp = gpar(col = "darkgrey", lwd = 1))

saveRDS(net, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet_clustered.rds")


pdf("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/main/DC.cluster.hmap.pdf",
    width= 5,
    height=4)
h.map

dev.off()


#net.mod= modularize(net, method= "leiden",resolution_parameter = 1.3)

links= igraph::as_data_frame(net, "edges")
nodes= igraph::as_data_frame(net, "vertices")

#Hfnet= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/hfnet.rds")
write_delim(links, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFnet_links.tsv", delim = "\t")
write_delim(nodes, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFnet_nodes.tsv",delim = "\t")

##
library(WriteXLS)
n.list= split(nodes, nodes$group_louv)
saveRDS(n.list, "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFnet_nodes.rds")

WriteXLS(n.list, SheetNames = names(n.list), ExcelFileName = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/output/HFnet_nodes.xlsx")


