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

set.seed(10)
######## run different cluster algos:

info = igraph::infomap.community(net,
                  e.weights = E(net)$weight,
                  v.weights = NULL,
                  nb.trials = 500,
                  modularity = T)

louv= cluster_louvain(net,
                      weights = E(net)$weight,resolution = 1.1

                        )

leid= cluster_leiden(net,objective_function = "modularity",
                                  weight= E(net)$weight,
                                  resolution_parameter = 0.9)
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


 louv$membership[louv$names=="hfpef"]
 louv$membership[louv$names=="hfref"]
 leid$membership[leid$names=="hfref"]
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
  theme_minimal()

# calc row max

apply(A.nmi,1, mean)

apply(A.adj,1, mean)


heatmap.df.adj= as.data.frame(A.adj) %>% rownames_to_column("algorithm1") %>%
  pivot_longer(cols = !algorithm1,
               names_to = "algorithm2",
               values_to= "adjusted.rand") %>%
  mutate(algorithm1 = factor(algorithm1, levels = ordered)) %>%
  mutate(algorithm2 = factor(algorithm2, levels = ordered))

p.adj = ggplot(heatmap.df.adj, aes(x= algorithm1, y= algorithm2, fill = adjusted.rand))+
  geom_tile()+
  scale_fill_gradient2(low = "white", high = "darkblue")+
  theme_minimal()


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
  labs(y= "number of nodes per cluster",
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

p.cluster= plot_grid( p.nmi+
                        labs(x="", y= "")+
                        theme(axis.text.x = element_text(angle = 70, hjust =1)),
                      p.adj+
                        labs(x="", y="")+
                        theme(axis.text.x = element_text(angle = 70, hjust =1))
                      , nrow= 1,
                      labels= c("D", "E"))

p.comb= plot_grid(p1, p.cluster, nrow= 2, rel_heights =c(0.8, 1) )
p.comb

pdf(file= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/figures/supp/comorbidity_net/cluster.algo.compare.pdf",
    width = 10, height= 8)
p.comb
dev.off()


