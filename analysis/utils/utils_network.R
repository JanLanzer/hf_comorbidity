## these are useful functions when analysing comorbidity networks:


## selects the largest connected component from an igraph object
#' @param net, igraph graph object to select the LCC from

selectLCC = function(net){
  # is not connected
  message(paste("your network was connected:", is_connected(net)))

  comp = igraph::components(net)
    # identify nodes that are not part of the LCC and delete from the network
  which(comp$csize == max(comp$csize))
  #select node  that don belong to the biggest cluster
  deleteme = V(net)[V(net)$name %in% names(comp$membership[comp$membership != which(comp$csize == max(comp$csize))])] #those nodes are not memberes of the big cluster #1
  #delete those nodes from the network
  net= delete_vertices(net, deleteme)

  message(paste("your network now is connected:", is_connected(net)))
  plot(hist(comp$csize, breaks= 100))
  return(net)

}

## get global network properties from igraph network
#' @param net, igraph graph object
#' @param directed, logical to analyse directed or undirected networks

getProps_global = function(net, directed = T){
  centr_score = data.frame("degree"= centr_degree(net, mode="all", normalized=T)$centralization,
                           "eigen"= centr_eigen(net, directed=directed, normalized=T)$centralization,
                           "betweenness"=centr_betw(net, directed=directed, normalized=T)$centralization,
                           "closeness"= centr_clo(net, mode="all", normalized=T)$centralization
                           )


  props = centr_score %>%
    gather("method","value") %>%
    dplyr::mutate(measure = "centrality") %>%
    dplyr::select(measure, method, value)

  props = rbind(props,
                c("diameter","mean_distance", mean_distance(net, directed=directed, unconnected = F)),
                c("diameter", "longest_geodesic",diameter(net, directed=directed)),
                c("transitivity", "global", transitivity(net, type="global")),
                c("edge_density", "percent", edge_density(net, loops=F)),
                c("degree_assortativity", "global",  assortativity_degree(net, directed = F))

  )%>%
    as_tibble %>%
    #dplyr::mutate(value = round(as.numeric(value),2))
  return(props )

}


# Node analysis -----------------------------------------------------------

# function to retrieve node properties including the rank within the network
#' @param net, igraph graph object
#' @param directed, logical to analyse directed or undirected networks
#' @param nodes, the name of the node to analyze. ATM only singular nodes possible

getProps_local = function(net,directed= T,  node){
  nodes= V(net)[V(net)$name == node]

  props = data.frame( measure = c("degree",
                                  "closeness",
                                  "betweenness",
                                  "transitivity",
                                  "strength"),
                      value= c( degree(net, node, mode= "all"),
                                closeness(net, node),
                                betweenness(net, node),
                                transitivity(net, type ="local", vids= c(node)),
                                strength(net, vids= node)),
                      rank = c( rank(sort(desc(degree(net, mode= "all"))))[node],
                                rank(sort(desc(closeness(net))))[node],
                                rank(sort(desc(betweenness(net))))[node],
                                which(sort(transitivity(net, type ="local")) == transitivity(net, type ="local", vids= c(node)))[1],
                                rank(sort(desc(strength(net))))[node]
                                )

                        )%>%
  as_tibble()
  if(directed ==T){
    hs <- hub_score(net, weights=NA)$vector
    as <- authority_score(net, weights=NA)$vector
    props= props %>% rbind(c("hub", hs[node], which(names(hs)== node)),
                           c("authority",as[node],which(names(as)== node)),
                           c("degree_out",degree(net, node, mode = "out"),rank(sort(desc(degree(net, mode ="out"))))[node]),
                           c("degree_in",degree(net, node, mode = "in"),rank(sort(desc(degree(net, mode ="in"))))[node])
                           )
  }
  props = props %>%
    mutate(value= round(as.numeric(value), 2)) %>%
    arrange(measure)
  return(props)
}


# function to recieve incoming and outcoming nodes (as tables and plots)
#' @param net, igraph objec
#' @param node, node of interest


get_edges = function(net, node){

  if(!exists("Phe_dic")){

    Phe_dic = readRDS("T:/fsa04/MED2-HF-Comorbidities/data/RWH_March2020/output_Jan/ICD10_labeled_phe.rds") %>%
      filter(pid %in% patients) %>%
      distinct(PheCode, Phenotype) %>%
      drop_na

  }


  dataframe = as_tibble(as_data_frame(net, what = c("edges")))

  outgoing = dataframe %>%
    filter(from == node) %>%
    arrange(desc(partial_cor)) %>%
    left_join(Phe_dic %>% rename(to= PheCode), by= "to")
  outgoing= outgoing %>%
    mutate(logp= ifelse(is.finite(logp),logp, max(outgoing$logp[is.finite(outgoing$logp)])))

  p1= ggplot(outgoing, aes(x= reorder(Phenotype, partial_cor),
                           y= partial_cor,
                           size = logp))+
    geom_point()+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    ggtitle(paste("outgoing edges from", node))


  incoming = dataframe %>%
    filter(to == node) %>%
    arrange(desc(partial_cor)) %>%
    left_join(Phe_dic %>% rename(from= PheCode), by= "from")

  incoming= incoming %>% mutate(logp= ifelse(is.finite(logp),
                                             logp,
                                             max(incoming$logp[is.finite(incoming$logp)])))


  p2= ggplot(incoming, aes(x= reorder(Phenotype, partial_cor),
                           y= partial_cor,
                           size = logp))+
    geom_point()+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    ggtitle(paste("incoming edges to",node))

  return(list("in_df"= incoming,"plot_in"= p2,"out_df" = outgoing, "plot_out"= p1))


}

# function to subset network to specific nodes
#' @param net, igraph graph object
#' @param order, numerical indicating the order of nodes to subset to
#' @param nodes, the name of the nodes to be included

get_subnet = function(net, order, nodes){

  selnodes <- V(net)[name %in% nodes]

  selegoV <- ego(net, order=order, nodes = selnodes, mode = "all", mindist = 0)

  subnet <- induced_subgraph(net,unlist(selegoV))
  return(subnet)
}

#function to compare two networks based on jaccard, worrks with nodes or edges
#' @param g1, graph 1
#' @param g2, graph 2
#' @param type, select vertex or edge to compare with jaccard
#'
jaccard_index<-function(g1,g2, type = "vertex") {
  if(type == "vertex"){
    g1_vert= igraph::as_data_frame(g1, what= "vertices")
    g2_vert = igraph::as_data_frame(g2, what= "vertices")

    u= union(g1_vert$name, g2_vert$name)
    i= intersect(g1_vert$name, g2_vert$name)

    jacc= length(i)/length(u)
    return(jacc)

  }else if(type== "edge"){
    g1_edg= igraph::as_data_frame(g1, what= "edges") %>% select(from,to) %>% as_tibble
    g2_edg = igraph::as_data_frame(g2, what= "edges")%>% select(from,to)%>% as_tibble

    # create an edgeID that is the same in both networks (sort the 2 nodes and than paste)
    g1_IDs= map2(g1_edg$from, g1_edg$to, function(x,y){
      paste(sort(c(x,y)),collapse="_")
    })
    g1_edg= g1_edg %>% mutate(edgeID = unlist(g1_IDs))

    g2_IDs= map2(g2_edg$from, g2_edg$to, function(x,y){
      paste(sort(c(x,y)),collapse="_")
    })
    g2_edg= g2_edg %>% mutate(edgeID = unlist(g2_IDs))

    # now calculate jaccard based on this edge id union and intersect
    u= union(g1_edg$edgeID, g2_edg$edgeID)
    i= intersect(g1_edg$edgeID, g2_edg$edgeID)

    jacc= length(i)/length(u)
    return(jacc)

  }else{
    message("Please select vertex or edge")
  }
}

# Network attributes -------------------------------------------------------------------------
network.df = function(seq, column,pids, links, data){

  netdfs= map(seq, function(x){

    #1. create networks
    if(column == "adj.p"){
      links_red = links[links[[column]]<x,] %>%
        select(disease1, disease2, estimate, partial_cor) %>%
        drop_na %>%
        filter(partial_cor>0) #important to keep out negative weights
    }else{
      links_red = links[links[[column]]>x,]
      # %>% ## define the sign !
      #   select(disease1, disease2, estimate, partial_cor) %>%
      #   drop_na %>%
      #   filter(partial_cor>0) #important to keep out negative weights
    }

    ## create network
    net= quick_base_net(links_red, pids, data,weight_col = column )

    ## if not connected, select LCC
    if(!is.connected(net)){
      message(paste("largest component selectet at value", as.character(x)))
      net= selectLCC(net)
    }
    # general network properties
    results_props= getProps_global(net, directed= F)

    #modularity
    louvain_partition <- igraph::cluster_louvain(net, weights = E(net)$weight)

    walk_partition= cluster_walktrap(net, weights = E(net)$weight, steps = 10,
                                     merges = TRUE, modularity = TRUE, membership = TRUE)


    x1 = c("modularity", "walk_partition", modularity(net, walk_partition$membership))
    x2=  c("modularity", "louvain_partition", modularity(net, louvain_partition$membership))
    x3 = c("n.nodes", "n.nodes", vcount(net))
    x4 = c("n.edges", "n.edges", ecount(net))
    rbind(results_props, x1,x2,x3,x4)

  })

  names(netdfs)= seq

  netdfs= melt(netdfs) %>%
    as_tibble() %>%
    mutate(value= as.numeric(value),
           L1= as.numeric(L1)) %>%
    pivot_wider(id_cols= c(method, value, L1), names_from= method, values_from = value) %>%
    rename( transitivity = global,
            edge_density=  percent)

  pls= map(colnames(netdfs)[-1], function(x1){
    ggplot(netdfs,aes_string(x= "L1", y= x1 ))+
      geom_point()+
      theme_minimal()+
      ggtitle(x1)+
      labs(x= column)

  }
  )

  return(list("df"= netdfs, "plots"= pls))
}

## select an node attribute and create a color attribute with a color assignment
#' @param net,
#' @param attr
#'
attr= "description"

create_color = function(net, attr, col.set = c("#fb8500",
                                               "#0099cc",
                                               "#f60000",
                                               "#fde800",
                                               "#00b600",
                                               "#00896e",
                                               "#0053c8",
                                               "#ff3399",
                                               "#b374e0",
                                               "#00cc99",
                                               "#fbc9e2",
                                               "#8a7437",
                                               "#2d3b61",
                                               "#612d3f",
                                               "#c9386b",
                                               "#357785",
                                               "#91b5bd",
                                               "#564e73")){

nodes= igraph::as_data_frame(net, what = "vertices")

  if(length(col.set)< length(unique(nodes[[attr]]))){
    return(message("col.set too small for number of categories, please expand col.set"))

  }else{
    df_color = data.frame("label" = unique(nodes[[attr]]),
                        "color" = col.set[1:length(unique(unique(nodes[[attr]])))]
    )

    color_dic = sapply(unique(nodes[[attr]]), function(x){
      color = df_color %>% filter(label == x) %>% pull(color)
    })
    color_vec= color_dic[nodes$group_cat]
    V(net)$color = color_vec

    return(net)
  }
}





# Network plotting --------------------------------------------------------

#function to export to cytoscape

R.to.Cyto = function(net, path){
  nodes = as_tibble(igraph::as_data_frame(net, what = "vertices"))
  edges = as_tibble(igraph::as_data_frame(net, what = "edges"))
  write.table(nodes, file=paste0(path,"_nodes.tsv"), quote=FALSE, sep='\t',row.names = F)
  write.table(edges, file=paste0(path,"_edges.tsv"), quote=FALSE, sep='\t',row.names = F)
  #write_delim(nodes, paste0(path,"_nodes.tsv"))
  #write_delim(edges, paste0(path, "_edges.tsv"), delim = "\t")
  #write.csv(nodes, paste0(path,"_nodes.csv"))
  #write.csv(edges, paste0(path, "_edges.csv"))
}

# function to get a quick vis net presentation:

quickvisnet= function(net, save_comp= T){
  require(visNetwork)
  data <- toVisNetworkData(net)

  data$edges= data$edges %>%
    mutate(width = weight*10) %>%
    #mutate(arrows = "to") %>%
    as_tibble#%>% filter(estimate >2)

  data$nodes = data$nodes %>%
    as_tibble  %>%
    dplyr::select(-label)%>%
    rename(value = size,
           label = description
    )

  network =visNetwork(data$nodes, data$edges, main = "",
                      height = "1000px", width = "100%") %>%
    visEdges(shadow = F,
             color = list(color = "darkgrey", highlight = "black")) %>%
    visPhysics(stabilization = T) %>%
    #visNodes()%>%
    visEdges(smooth = T) %>%
    addFontAwesome()%>%
    visLegend(useGroups = F,
              #addNodes = df_color,
              ncol = 1, stepY=40) %>%
    #visHierarchicalLayout(direction = "LR", levelSeparation = 100)%>%
    visOptions(highlightNearest =  list(enabled = TRUE, algorithm = "hierarchical",
                                        degree = list(from = 1, to = 1)), collapse = TRUE) %>%
    visInteraction(dragNodes = T,
                   dragView = T,
                   zoomView = T)

  if(save_comp==T){
    network =network %>%
      visPhysics(stabilization = F) %>%
      visEdges(smooth = FALSE) %>%
      visInteraction(dragNodes = F)
    return(network)
  }
  return(network)

}


# community analysis ------------------------------------------------------
#function to create a basic network (UW),with partial cor as weights, node size as log(freq)
#' @param links, link table
#' @param pids, patient ids used to create the network
#' @param data, full data table of all diagnostic entries in patients
#'
quick_base_net = function(links, pids, data, weight_col = "corr.phi"){
  #get nodes stats
  nodes=  unique(c(links %>% pull(disease1),
                   links %>% pull(disease2)))

  freqs = disease_frequencies(pids, data, "PheCode")%>%
    filter(PheCode %in% nodes) %>%
    left_join(Phe_dic) %>%
    distinct(PheCode, freq, Phenotype, category)

  # 2. create network
  net= graph_from_data_frame(d= links,
                             vertices = freqs$PheCode,
                             directed = F)
  ## attributes:
  V(net)$size= log(freqs$freq)
  V(net)$description = freqs$Phenotype
  V(net)$group_cat= freqs$category

  #edges, weight
  if (weight_col == "pcorr.corpor"){
    E(net)$weight = E(net)$pcorr.corpor
    }else if(weight_col == "cor2pcor"){
      E(net)$weight = E(net)$pcor.corpor
    }else if(weight_col == "rho"){
      E(net)$weight = E(net)$rho
    }else if(weight_col == "corr.phi"){
      E(net)$weight = E(net)$corr.phi
    }

  return(net)
}


#function to perform louvain_partion based on maximal modularity
#' @param , net =igraph net, (UW) to for community analysis

module= function(net){

  #lovain
  louvain_partition <- igraph::cluster_louvain(net, weights = E(net)$weight)

  print(modularity(net, louvain_partition$membership))

  V(net)$group_louv= louvain_partition$membership
  nodss = as_tibble(igraph::as_data_frame(net, what = "vertices"))

  # plot clusters
  plots =map(sort(unique(nodss$group_louv)), function(x){
    clust= nodss %>% filter(group_louv == x)
    clust %>%
      count(group_cat) %>%
      mutate(prop = n/sum(n)) %>%
      arrange(desc(n)) %>%
      ggplot(., aes(x= reorder(group_cat, -prop), y= prop, fill = group_cat))+
      geom_col()+
      labs(x= "",y = "")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle=60, hjust= 1),
            legend.position = "none",
            plot.title = element_text(size = 10, face = "bold"))+
      ggtitle(paste("cluster:",as.character(x), "; ndisease:", as.character(dim(clust)[1])))
  })
  plot =plot_grid(plotlist= plots, nrow= 3)

  if("428.1" %in% nodss$name){
  HF.g= nodss %>%  filter(name == "428.1") %>% pull(group_louv)
  hfsub = nodss %>% filter(group_louv == HF.g) %>% pull(name)
  net2 = induced_subgraph(net,hfsub , impl = c("auto"))

  print(paste0("hf cluster has number:", HF.g))
  }else{
    net2 = NULL
  }

  return(list(plot = plot, fullnet = net, HFmodule= net2))
}




# Over representation analysis --------------------------------------------

double_enrich = function(fullnet, pids.list, genesets){

  #extract nodes and edges from the network:
  edge.module_net = as_tibble(igraph::as_data_frame(fullnet, what = "edges"))
  nodes.module_net = as_tibble(igraph::as_data_frame(fullnet, what = "vertices"))

  edge_IDs= map2(edge.module_net$from, edge.module_net$to, function(x,y){
    paste(sort(c(x,y)),collapse="_")
  })

  edge.module_net= edge.module_net %>% mutate(edgeID = unlist(edge_IDs))

  clusters= sort(unique(nodes.module_net$group_louv))


  # for every cluster calculate 2 enrichments,  return both
  listed.gsea= lapply(clusters, function(clust){

    #clust= 4
    print(paste0("enriching in cluster #",clust))

    #get all nodes in cluster x
    module.x= nodes.module_net %>%
      filter(group_louv == clust) %>%
      pull(name)

    # filter edges of the network to edges of module (keep only edges WITHIN)
    g1_edg = edge.module_net %>%
      filter(from %in% module.x & to %in% module.x)

    # #into a sorted graph ID , this is unqiue for a disease pair
    # g1_IDs= map2(g1_edg$from, g1_edg$to, function(x,y){
    #   paste(sort(c(x,y)),collapse="_")
    # })
    #
    # g1_edg= g1_edg %>% mutate(edgeID = unlist(g1_IDs))
    #
    # create the ranked vector based on the edge ID

    g1_edg = g1_edg %>% arrange(desc(weight))
    ranked_vec= g1_edg%>% pull(weight)
    names(ranked_vec) = as.character(g1_edg$edgeID)

    ##### 1. GSEA, enrich patients in  each module
    ## "genesets" =patients
    ## "genes" = EdgeIDs of a patients comorbidity profile
    ## ranks = EdgeIDS of Edges in a module , ranked by weight (part_corr)
    gseares= fgsea(pathways=genesets,
                   stats= ranked_vec,
                   scoreType= "pos",
                   minSize =2,
                   maxSize=500)%>%
      as_tibble()

    ##### 2. GSEA, enrich patient labels in patients
    ## "genesets" =patient cohorts (male, female etc.)
    ## "genes" = patient IDs (belonging to that cohort)
    ## ranks = enriched patients in a network module, ranked by enrichment pval.

    gseares= gseares %>% drop_na %>% arrange(desc(NES)) %>% mutate(pval = -log10(pval))
    ranked_vec2 = gseares %>% pull(pval)

    names(ranked_vec2) = gseares$pathway


    if(length(ranked_vec2)>= 10){ # only perform a second enrichment if a pat. vector of at least 10 exists

      gsea_final = fgsea(pathways= pids.list, stats= ranked_vec2, scoreType= "pos")

      plot.list= map(unique(gsea_final$pathway), function(x){
        plotEnrichment(pathway= as.character(unlist(pids.list[x])), stats = ranked_vec2)+
          labs(title=x)
      })

      p.combined = plot_grid(plotlist = plot.list, ncol = 1)

    }else{
      gsea_final = NULL
      p.combined= NULL
    }

    return(list("2nd"= gsea_final, "1st"= gseares, "enrichplots"= p.combined))

  })

  names(listed.gsea)  = clusters

  return(listed.gsea)
}


plot_enrich_res= function(enrich_res){
  first_enrichement = map(names(enrich_res), function(i){
    enrich_res[[i]]$`1st`
  })

  names(first_enrichement) = names(enrich_res)

  x= map(first_enrichement, function(x){dim(x)[1] })

  ## get full numbers of vectors:
  real_nums= sapply(names(first_enrichement), function(x){
    df= first_enrichement[[x]]
    x1 = df %>% filter(pathway %in% pids.list$hfref)
    x2= df %>% filter(pathway %in% pids.list$hfpef)
    x3= df %>% filter(pathway %in% pids.list$hfmref)
    #c(x1,x2)
    return(c(dim(x1)[1], dim(x2)[1], dim(x3)[1]))
  })

  real_nums= real_nums %>%
    as_tibble() %>%
    mutate(group = c("hfref", "hfpef", "hfmref")) %>%
    pivot_longer(cols= -group)

  real_nums2= real_nums %>%
    group_by(name) %>%
    summarise("pat_vec"= sum(value)) %>%
    left_join(real_nums) %>%
    mutate(name= factor(name, levels= names(enrich_res)))

  plot.vectorlength = ggplot(real_nums2%>% filter(pat_vec >10),
                             aes(x= name,
                                 y= value,
                                 fill = group))+
    geom_col()+
    #scale_y_log10()+
    coord_flip()+
    scale_fill_manual(values =  cols.nice)

  second_enrichments= map(names(enrich_res), function(i){
    enrich_res[[i]]$`2nd`
  })

  names(second_enrichments)= names(enrich_res)

  plot.data =bind_rows(second_enrichments,.id = "clust")

  plot.data$stars <- cut(plot.data$pval, breaks=c(-Inf, 0.05, 0.01, 0.1, Inf), label=c("***", "**", "*", ""))


  p.hmap = plot.data %>%
    mutate(clust= factor(clust, levels= names(enrich_res))) %>%
    ggplot(., aes(x= pathway, y= reorder(clust,clust), fill= NES))+
    geom_tile()+
    scale_fill_gradient(low = "white", high = "red" )+#,limits = c(0,max(plot.data$NES)))+
    geom_text(aes(label=stars), color="black", size=5)+
    theme_minimal()+
    labs(x= "")

  plot_grid(p.hmap,
            plot.vectorlength,
            rel_widths = c(0.8,1),
            nrow = 1) #+ labs(y= "fraction of hfpef+hfref patients in first enrichment"))

}
