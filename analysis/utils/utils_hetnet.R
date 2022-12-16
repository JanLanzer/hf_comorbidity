## UTIL  funcs for creating het net
# load additional data
#patients =readRDS(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/time_range_patientIDs.rds")
HF_allcause= c("I11.0", "I13.0", "I13.2", "I25.5", "I42.0", "I42.5", "I42.8", "I42.9", "I50.0", "I50.1", "I50.9")



# randowm walk function ---------------------------------------------------

# multilayer:

run_multi_RW= function(dd,
                       dg,
                       gg,
                       p.val= 0.1,
                       weight_dd= 0.1,
                       weight_gg= 0.1,
                       weight_dg= 0.1,
                       seed_diseases,
                       prefiltered=F,
                       r= 0.7,
                       weighted= F,
                       ...){

  require(RandomWalkRestartMH)

  if(!prefiltered){
    dd =  dd %>%
      filter(fisher.p.adj<p.val) %>%
      rename(nodeA= disease1,
             nodeB= disease2,
             weight= corr.phi) %>%
      select(nodeA, nodeB, weight) %>%
      rowwise()%>%
      mutate(weight_sc= weight,
             type= "d-d") %>%
      filter(weight>weight_dd)

  }


  # filter bipartite disease-gene for for nodes in the disease and gene layers:
  #print(table(unlist(sets) %in% c(gg$nodeA, gg$nodeB, dg$nodeA)))

  dd = dd %>% filter(weight_sc > weight_dd)
  gg = gg %>% filter(weight_sc > weight_gg)
  dg= dg %>% filter(nodeA %in% c(gg$nodeA, gg$nodeB),
                    nodeB %in% c(dd$nodeA, dd$nodeB))%>% filter(weight_sc > weight_dg)

  #print(table(unlist(sets) %in% c(gg$nodeA, gg$nodeB, dg$nodeA)))
  #create two networks:
  ppi_net = graph_from_data_frame(gg, directed = F)
  dd_net = graph_from_data_frame(dd, directed = F)

  if(!weighted){
    ppi_net= delete_edge_attr(ppi_net, "weight")
    dd_net= delete_edge_attr(dd_net, "weight")

  }
  ppi_net_multiplex= create.multiplex(list("PPI"= ppi_net))
  dd_net_multiplex <- create.multiplex(list("DD"=dd_net))
  #?create.multiplex
  # create het net:
  if(!weighted){
    PPI_Disease_Net <- create.multiplexHet(ppi_net_multiplex,
                                           dd_net_multiplex,
                                           as.data.frame(dg[,c(1,2)]), "disease")

  }else{
    PPI_Disease_Net <- create.multiplexHet(ppi_net_multiplex,
                                           dd_net_multiplex,
                                           as.data.frame(dg[,c(1,2,3)]), "disease")

  }
    # PPI_Disease_Net <- create.multiplexHet(ppi_net_multiplex,
    #                                        dd_net_multiplex,
    #                                        as.data.frame(dg[,c(1,2)]), "disease")
    #
  PPIHetTranMatrix <- compute.transition.matrix(PPI_Disease_Net)

  diseases = unique(c(dd$nodeA, dd$nodeB))

  RWRH_PPI_Disease_Results <-
    Random.Walk.Restart.MultiplexHet(x= PPIHetTranMatrix,
                                     MultiplexHet_Object = PPI_Disease_Net,
                                     #SecondNet_Seed_Nodes= unique(seed_diseases[seed_diseases %in% diseases]),
                                     #Multiplex_Seed_Nodes= c(),
                                     Multiplex1_Seeds= c(),
                                     Multiplex2_Seeds = unique(seed_diseases[seed_diseases %in% diseases]),
                                     r=r)



   gene_res= RWRH_PPI_Disease_Results$RWRMH_Multiplex1 %>%
     rename(value= Score,
            name= NodeNames)%>% as_tibble()

  return(gene_res)
}

# break function up in two parts to save computational time:
# 1 only create multiplex nets

Compute_multiplexes= function(dd,
                              ppi,
                              gg= NULL,
                              p.val= 0.1,
                              weight_dd= 0.1,
                              weight_gg= 0.1,
                              weight_dg= 0.1,
                              weight_ppi= 0.1,

                              prefiltered=F,

                              weighted= F,
                              ...){

  require(RandomWalkRestartMH)

  if(!prefiltered){
    dd =  dd %>%
      filter(fisher.p.adj<p.val) %>%
      rename(nodeA= disease1,
             nodeB= disease2,
             weight= corr.phi) %>%
      select(nodeA, nodeB, weight) %>%
      rowwise()%>%
      mutate(weight_sc= weight,
             type= "d-d") %>%
      filter(weight>weight_dd)

  }

  dd = dd %>% filter(weight_sc > weight_dd)

  ppi = ppi%>% filter(weight_sc > weight_ppi)

  #print(table(unlist(sets) %in% c(gg$nodeA, gg$nodeB, dg$nodeA)))
  #create two networks:

  ppi_net = graph_from_data_frame(ppi, directed = F)

  dd_net = graph_from_data_frame(dd, directed = F)

  if(!weighted){

    ppi_net= delete_edge_attr(ppi_net, "weight")

    dd_net= delete_edge_attr(dd_net, "weight")
  }


  # prepare the gene-gene net ( if available in the same manner:)

  if(!is.null(gg)){

    gg = gg %>% filter(weight_sc > weight_gg)

    func_net= graph_from_data_frame(gg, directed = F)

    if(!weighted){
      func_net= delete_edge_attr(func_net, "weight")
    }

    ppi_net_multiplex= create.multiplex(list("PPI"= ppi_net,
                                             "GG"= func_net))
  } else {
    ppi_net_multiplex= create.multiplex(list("PPI"= ppi_net))
  }

  dd_net_multiplex <- create.multiplex(list("DD"=dd_net))

  return(list("ppi"= ppi_net_multiplex,
              "dd"= dd_net_multiplex))
}



# 2 add d-g links and compute het net + run multi layer RW

run_multi_RW_2= function(multi_list,
                       dg,
                       weight_dg= 0.1,
                       seed_diseases,
                       r= 0.7,
                       weighted= T,
                       ...){

  require(RandomWalkRestartMH)

  if(is.igraph(multi_list$ppi[[2]])){

    genes= unique(c(V(multi_list$ppi[[1]])$name,
             V(multi_list$ppi[[2]])$name
    ))
  } else {
    genes= unique(c(V(multi_list$ppi[[1]])$name
    ))
  }

  diseases= V(multi_list$dd$DD)$name

  dg_fil= dg %>% filter(nodeA %in% genes,
                    nodeB %in% diseases)#%>%
    #filter(weight_sc > weight_dg)

  # create het net:
  if(!weighted){
    PPI_Disease_Net <- create.multiplexHet(multi_list$ppi,
                                           multi_list$dd,
                                           as.data.frame(dg_fil[,c(1,2)]), "disease")

  }else{
    PPI_Disease_Net <- create.multiplexHet(multi_list$ppi,
                                           multi_list$dd,
                                           as.data.frame(dg_fil[,c(1,2,3)]), "disease")

  }

  #message(PPI_Disease_Net$Multiplex1$Number_of_Nodes_Multiplex)
  #message(PPI_Disease_Net$Multiplex2$Number_of_Nodes_Multiplex)

  PPIHetTranMatrix <- compute.transition.matrix(PPI_Disease_Net)

  RWRH_PPI_Disease_Results <-
    Random.Walk.Restart.MultiplexHet(x= PPIHetTranMatrix,
                                     MultiplexHet_Object = PPI_Disease_Net,
                                     Multiplex1_Seeds= c(),
                                     Multiplex2_Seeds = unique(seed_diseases[seed_diseases %in% diseases]),
                                     r=r)



  gene_res= RWRH_PPI_Disease_Results$RWRMH_Multiplex1 %>%
    rename(value= Score,
           name= NodeNames)%>% as_tibble()

  return(gene_res)
}


#aggregated:
get_genes_for_net= function(.net= NULL,
                            diseaseoi,
                            gg= gg,
                            weighted= T,
                            weight_dg= 0.1,
                            weight_gg= 0.1,
                            weight_dd= 0.1,
                            gamma.= 0.7,
                            .links= NULL,
                            col.weight= "corr.phi",
                            ...){

  # connect that net to the rest
  if(!is.null(.net)){
    edges.dd= igraph::as_data_frame(.net , what= "edges") %>%
      as_tibble() %>%
      dplyr::rename(nodeA= from,
                    nodeB= to)%>%
      dplyr::select(nodeA, nodeB, weight)%>%
      mutate(weight_sc= weight,
             type= "d-d" ) %>%
      filter(weight>weight_dd)

  }else if(!is.null(.links)){
    edges.dd= .links %>%
      dplyr::rename(nodeA= disease1,
                    nodeB= disease2)%>%
      rename(weight= {{col.weight}}) %>%
      dplyr::select(nodeA, nodeB, weight)%>%
      mutate(weight_sc= weight,
             type= "d-d" ) %>%
      filter(weight>weight_dd)
  }
  #filter(nodeA %in% dis. & nodeB %in% dis.)
  # get disease nodes in this net:
  disease= unique(c(edges.dd$nodeA, edges.dd$nodeB))
  genes= unique(c(gg$nodeA, gg$nodeB))

  #filter bipartite relationships
  dg= dg %>% filter(nodeA %in% genes,
                    nodeB %in% diseaseoi) %>%
    filter(weight> weight_dg)

  gg = gg %>% filter(weight_sc> weight_gg)


  .hetnet= rbind(gg, dg %>% filter(nodeB %in% disease), edges.dd)
  .gX= graph_from_data_frame(.hetnet, directed= F)

  if(weighted){
    A_g1= get.adjacency(.gX, attr = "weight_sc")
  }else{
    A_g1= get.adjacency(.gX, attr = NULL)
  }

  x1= RWR(as.matrix(A_g1),
          gamma = gamma.,
          ind.positives = diseaseoi[diseaseoi %in% disease],
          norm = T)

  enframe(x1$p) %>% arrange(desc(value)) %>% filter(!name %in% disease)

}

# processing layers -------------------------------------------------------


process_morbinet= function(){
  #load morbidnet
  morbinet= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/morbinet_edges.csv")

  # translate icd3 to PheCode
  morbinet= morbinet %>% left_join(data %>% dplyr::distinct(PheCode, icd3) %>%
                                     rename(PheCode_A = PheCode,
                                            Disease_A= icd3), by= "Disease_A") %>%
    left_join(data %>% dplyr::distinct(PheCode, icd3) %>% rename(PheCode_B = PheCode,
                                                          Disease_B= icd3), by= "Disease_B") %>% as_tibble()
  ## disease disease edges
  dd_edges = morbinet %>%filter(OR>2) %>% rename(nodeA= PheCode_A,
                                                 nodeB= PheCode_B,
                                                 weight= OR)%>%
    select(nodeA, nodeB, weight) %>%
    drop_na

  # scale the weight /max
  dd_edges=  dd_edges %>%
    filter(nodeA != nodeB) %>%
    group_by(nodeA, nodeB) %>%
    summarise(weight = mean(weight))%>%
    mutate(weight_sc= weight/max(weight))

}


get_net_link= function(links.sub, pids, col.weight= "corr.phi",col.weight.cutoff= 0.1, fisherp= 0.05){
  .net= quick_base_net(links.sub %>% filter(fisher.p.adj<fisherp,
                                            #odds.ratio>2,
                                            !!as.symbol(col.weight)>col.weight.cutoff) %>%
                         rename(weight= {{col.weight}}) , weight_col = "weight",pids= pids,  data= data)

  E(.net)$weight

  selectLCC(.net)

  .net= create_color(.net, "group_cat")

  .net= module(.net)

}

create_edge_table_for_cohorts= function(pid.list, data, topn_disease,phecodes){

  .links = map(pid.list, function(pids){

    .links= create.links(pids= pids,
                         phecodes = phecodes,
                         data= data)


  })
  message("links created for all cohorts")
  names(.links)=names(pid.list)

  nets= map(names(.links), function(x){
    message("processing")
    print(head(.links[[x]]))
    print(pid.list[[x]])
    get_net_link(links.sub = .links[[x]],
                 pids = pid.list[[x]]
                 )
    }
    )
  names(nets)=names(pid.list)

  return(list(nets= nets,
              links= .links))

}


################## LAYER 2-3 gene- disease
#load data and dictionary:
#data = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/ICD10_labeled_phe.rds")
#pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/pidslist_v2021.rds")

process_gd = function(weight_cutoff= 0, clingen= "none"){
  icd3_dic= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/icd3_disgenet_dictionary.rds")
  # load the manual curated mappings:

  ## get the table from  :
  # "https://docs.google.com/spreadsheets/d/1cZ7WZz8MF5fchrVnUrx4r4Sps62Y3O9gK4k88gVZZsQ/edit#gid=985996458"
  icd3_dic_manual = read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/disgenet.mapping.for.manual - disgenet.mapping.for.manual.csv")%>%
    as_tibble

  manual_dic= rbind(icd3_dic_manual %>%
                      dplyr::select(icd3, X.1)%>%
                      dplyr::rename(diseaseId= X.1),
                    icd3_dic_manual %>%
                      dplyr::select(icd3, X.2)%>%
                      dplyr::rename(diseaseId= X.2))%>%
    dplyr::filter(diseaseId!= "")%>%
    dplyr::mutate(diseaseId= str_replace_all(diseaseId, pattern = " ", ""))

  icd3_dic= rbind(icd3_dic, manual_dic)

  disgenet= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/all_gene_disease_associations.tsv", sep ='\t') %>%
    as_tibble()

  # combine the data:
  # patientmap_df= data %>%
  #   filter(pid %in% pids.list$hf_all) %>%
  #   distinct(pid, icd3, PheCode) %>%
  #   drop_na %>%
  #   left_join(icd3_dic)

  if(clingen!= "none"){
    disgenet= disgenet %>% dplyr::filter(grepl(clingen, source))

  }

  disgenetmap_df= disgenet %>%
    dplyr::filter(score>weight_cutoff) %>%
    dplyr::select(geneSymbol, diseaseId, diseaseName, score)

  #pat_dis_gene = patientmap_df %>% left_join(disgenetmap_df, by= "diseaseId")

  # add PheCode to disgenet
  disgenet.map = disgenetmap_df %>%
    left_join(icd3_dic, by= "diseaseId") %>%
    dplyr::filter(!is.na(icd3)) %>%
    left_join(data %>% dplyr::distinct(icd3, PheCode))

  # create edge table:
  dg_edges = disgenet.map %>%
    dplyr::distinct(geneSymbol, PheCode, score)  %>%
    dplyr::rename(nodeA= geneSymbol, nodeB= PheCode,weight= score) %>%
    drop_na %>%
    dplyr::mutate(weight_sc= weight/max(weight)) %>%
    dplyr::mutate(type = "d-g")%>%
    dplyr::filter(weight_sc>weight_cutoff)
}

#get gwas catalogue based diseae gene mappings if a dg is provided it will be merged and returned
process_gwas= function(dg=NULL){
  gwas.phe= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/gene-disease/Gwas.cat.processed.rds")

  if(!is.null(dg)){

    dg.gwas= gwas.phe %>% rename(nodeA= REPORTED.GENE.S.,
                                 nodeB= PheCode)%>%
      select(nodeA, nodeB) %>%
      mutate(weight= median(dg$weight),
             weight_sc= weight,
             type = "d-gwas") %>%
      drop_na()

    dg.comb= rbind(dg, dg.gwas)

    dg.comb= dg.comb %>%
      distinct(nodeA, nodeB)%>% mutate(weight= median(dg$weight),
                                                          weight_sc= weight,
                                                          type = "d-gwas")

  }
}
################## LAYER 3 gene-gene (PPI)

process_omnipath = function(weight_cutoff= 0,
                       n_res= 2,
                       n_cur= 2){


  ppi= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/int_omni_12072021.rds")

  ## PPI


  ppi_red= ppi %>% filter(n_resources >= n_res,
                          curation_effort>= n_cur)

  gg_edges = ppi_red %>%
    dplyr::select(source_genesymbol, target_genesymbol, n_resources) %>%
    dplyr::rename(nodeA= source_genesymbol, nodeB= target_genesymbol, weight= n_resources) %>% drop_na %>%
    dplyr::mutate(weight_sc= weight/max(weight)) %>%
    dplyr::mutate(type = "g-g")%>%
    filter(weight_sc>weight_cutoff)

}

process_ppi_v2 = function(weight_cutoff= 0){

  ppi = read_table2("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/PPI_2016-11-23.gr.txt",
                   col_names = F)
  colnames(ppi)= c("nodeA", "nodeB")
  gg_edges = ppi %>%
    #mutate(weight_sc= weight/max(weight)) %>%
    mutate(type = "g-g",
           weight= 1,
           weight_sc= 1)
}

process_human_ppi_multiverse = function(weight_cutoff= 0){

  ppi = read_table2("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/HI-union.tsv",
                    col_names = F)
  colnames(ppi)= c("nodeA", "nodeB")


  gene_translate= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/biomart_genetranslation.rds" )%>%
    as_tibble()

  ppi= ppi %>%
    left_join(., gene_translate %>% dplyr::rename(nodeA= ensembl_gene_id ,
                                                             nodeA1= hgnc_symbol))%>%
    left_join(., gene_translate %>% dplyr::rename(nodeB= ensembl_gene_id ,
                                                  nodeB1= hgnc_symbol)) %>%
    drop_na %>%
    dplyr::select(nodeA1, nodeB1) %>%
    dplyr::rename(nodeA= nodeA1,
                  nodeB= nodeB1)
  gg_edges = ppi %>%
    #mutate(weight_sc= weight/max(weight)) %>%
    dplyr::mutate(type = "g-g",
           weight= 1,
           weight_sc= 1)


}

process_ppi_lit= function(){
  ppi= read_table2("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Lit-BM.tsv", col_names = F)
  colnames(ppi)= c("nodeA", "nodeB")

  gene_translate= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/biomart_genetranslation.rds" )%>%
    as_tibble()

  ppi= ppi %>% left_join(., gene_translate %>% dplyr::rename(nodeA= ensembl_gene_id ,
                                                             nodeA1= hgnc_symbol))%>%
    left_join(., gene_translate %>% dplyr::rename(nodeB= ensembl_gene_id ,
                                                  nodeB1= hgnc_symbol)) %>%
    drop_na %>%
    select(nodeA1, nodeB1) %>%
    dplyr::rename(nodeA= nodeA1,
                  nodeB= nodeB1)
  gg_edges = ppi %>%
    #mutate(weight_sc= weight/max(weight)) %>%
    mutate(type = "g-g",
           weight= 1,
           weight_sc= 1)
}
# load data from study: https://www.nature.com/articles/s41467-017-01747-2#Sec2

process_human_net= function(){
  #(http://www.inetbio.org/humannet)
  human_net= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/HumanNet-XN.tsv",
                      sep = "\t",
                      skip = 0)%>%as_tibble()

  gene_tran= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/biomart_genetranslation_human_net.rds")%>%
    select(-ensembl_gene_id)

  human_net %>% left_join(gene_tran%>% rename(EntrezGeneID1=entrezgene_id,
                                              nodeA= hgnc_symbol))%>%
    left_join(gene_tran%>% rename(EntrezGeneID2=entrezgene_id,
                                  nodeB= hgnc_symbol)) %>%
    rename(weight= LLS) %>% select(nodeA, nodeB, weight)%>%
    mutate(type = "g-g",
           weight_sc= weight/max(weight))

}

process_GO_layer= function(layer){
  if(layer=="BP"){
    x= read.delim("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/edge_list_raredisease_pred/GOBP.tsv") %>%
      as_tibble()

  }else if(layer== "MF"){
    x= read.delim("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/edge_list_raredisease_pred/GOMF.tsv")%>%
      as_tibble()

  }

  x%>% dplyr::rename(nodeA= name1,
              nodeB= name2)%>%
    mutate(weight= 1,
           weight_sc= 1,
           type= "g-g")
}

get_human_heart_proteome=function(){
  heart.prot= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/proteinGroups_Human Heart Map.csv", sep = ";")
  #print(head(x))
  unique(heart.prot$Leading.Gene.Name)


}

get_GTEX_heart_genes = function(){
  x = gzfile("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")
  y= read.table(x, skip = 2, sep = '\t')
  colnames(y)= y[1,]
  y= y[-1,]

  heart.gtex= as_tibble(y) %>% mutate(Heart.A= as.numeric(`Heart - Atrial Appendage`),
                                      Heart.V= as.numeric(`Heart - Left Ventricle`))%>%
    dplyr::select(Description, Heart.A, Heart.V)


  heart.genes = heart.gtex %>% filter(Heart.A>1 | Heart.V >1) %>% pull(Description)

}


process_pathways= function(){
  data("Pathway_Network")
  pg_edges= igraph::as_data_frame(Pathway_Network, what= "edges")
  pg_edges= pg_edges %>% mutate(weight= 0.9,
                                weight_sc= weight,
                                type= "p-g") %>%
    dplyr::rename(nodeA= from,
                  nodeB = to)

}

#remove doubled edges, loops and select largest connected component
simplify_link_table= function(links, type = "d-d"){
  links.red= graph_from_data_frame(links, directed = F) %>% selectLCC(.) %>%
    igraph::simplify(remove.multiple = T, remove.loops = T) %>%
    igraph::as_data_frame(.,  "edges") %>%
    dplyr::rename(nodeA = from,
                  nodeB= to
    )%>%
    dplyr::mutate(weight_sc= 1,
                  type= type)

}

# validation mapping functions --------------------------------------------


randomize_links_by_layer= function(links, rewireprob= 0.9){

  .gX= graph_from_data_frame(links, directed= F)

  g.rand= rewire(.gX, with = keeping_degseq(niter= 100000))
  #E(g.rand)$weight_sc= runif(length(E(g.rand)), 0,1)
  E(g.rand)$weight= sample(E(.gX)$weight)

  g.rand2= rewire(.gX, with = each_edge(prob= rewireprob))
  E(g.rand2)$weight= sample(E(.gX)$weight)

  return(list(
    "random.degree"= igraph::as_data_frame(g.rand, "edges") %>%
      dplyr::rename(nodeA= from,
             nodeB= to) %>%
      arrange(desc(weight_sc)) %>%
      as_tibble(),
    "random.complete"= igraph::as_data_frame(g.rand2, "edges")%>%
      dplyr::rename(nodeA= from,
             nodeB= to)%>%
      arrange(desc(weight_sc))%>%
      as_tibble(),
    "real"= links
  )
  )
}

table_to_links= function(table,
                         p.val= 0.05,
                         weight_dd= 0.1, absolute_filter= F,
                         weight_col = "corr.phi"){
  table =  table %>%
    dplyr::filter(fisher.p.adj<p.val) %>%
    dplyr::rename(nodeA= disease1,
           nodeB= disease2,
           weight= {{weight_col}}) %>%
    dplyr::select(nodeA, nodeB, weight) %>%
    dplyr::mutate(weight_sc= weight,
           type= "d-d")
  if(absolute_filter){
    table= table  %>%
      dplyr::filter(abs(weight)>weight_dd)%>%
      dplyr::arrange(desc(weight))

  }else{
    table= table %>%
    dplyr::filter(weight>weight_dd)%>%
    dplyr::arrange(desc(weight))
  }

  return(table)
}

get.auc.for.combi= function(layer_name, layer_table, seed_diseases= dis.random,
                            val.set, ...){
  #assign all the real values
  xdd =randomize_links_by_layer(dd)
  xdg= randomize_links_by_layer(dg)
  xgg= randomize_links_by_layer(gg)

  dd= xdd$real
  gg= xgg$real
  dg= xdg$real

  #re-assign the layer that was randomized
  assign(x = layer_name, value = layer_table)

  g.rand = run_multi_RW(dd = dd,
                        dg = dg,
                        gg = gg,
                        seed_diseases = dis.hf.lr,
                        prefiltered = T, weight_dd = 0, weight_gg = 0, weight_dg = 0,
                        r= 0.9  ) %>%
    validate_results2(gene_results = .,set = val.set)

  g.rand$roc$auc

}

xxx= function(dd, dg, gg){
  xdd = randomize_links_by_layer(dd)
  xdg= randomize_links_by_layer(dg)
  xgg= randomize_links_by_layer(gg)

  linkS= list("dd"= xdd,
              "gg"= xgg,
              "dg"= xdg)

  res= sapply(names(linkS), function(x){

    sapply(names(linkS[[x]]), function(y){
      print(x)
      print(y)
      get.auc.for.combi(x, linkS[[x]][[y]])

    })
  })
}

random_geneset_prediction= function(gene_results, genes){
  sample_set= sample(x = genes, size = length(unlist(sets)), replace = F)
  val.r=validate_results2(sample_set, gene_results = gene_results)
  val.r$roc$auc
}



get_random_disease_set = function(dd,
                                  n.sample = 20
){
  disease= unique(c(dd$nodeA, dd$nodeB))

  sample(disease, n.sample, replace= F)

}



# validation via gene set--------------------------------------------------------------


get_Disgenet_overlaps= function(dg){

  df= dg %>% distinct(nodeB, nodeA)  %>% drop_na()

  # data transformation into n x p Matrix
  df= split(as.character(df$nodeB), df$nodeA) %>%
    lapply(function(x) strsplit(x, "-")) %>%
    mtabulate()%>%
    qdapTools::matrix2df("genes") %>%
    column_to_rownames("genes")

  jacc.dist.phecodes=
    dist(t(df), method= "binary") %>%
    as.matrix()
}

load_validation_genes = function(disgenet_value= 0.4, top_reheat= 500){
  #genes_disgenet =
  disgenet= read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/all_gene_disease_associations.tsv", sep ='\t') %>%
    as_tibble()


  disgenet_hflinks= disgenet %>%
    filter(score>disgenet_value) %>%
    dplyr::select(geneSymbol, diseaseId, diseaseName, score) %>%
    filter(grepl("C0018801", diseaseId)) %>%
    pull(geneSymbol) %>%
    unique()
  disgenet_hflinks

  #genes phe=
  phecode = read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/phewas-catalog.csv")
  phewas_hflinks= as_tibble(phecode)%>% filter(grepl("428", phewas.code))%>%
    filter(p.value<0.05, odds.ratio>1)

  kegg_DCM = c("LMNA", "CDCD1", "CDDC", "CMD1A", "CMT2B1", "EMD2", "FPL", "FPLD", "FPLD2", "HGPS",
               "IDC", "LDP1", "LFP", "LGMD1B", "LMN1", "LMNC", "LMNL1", "MADA", "PRO1")
  #source= https://www.nature.com/articles/s41467-019-13690-5/tables/1
  naturecomm_GWAS= c("FTO", "ATXN2", "BAG3", "SYNPO2L", "AGAP5", "ABO", "SURF1", "LPA", "CDKN1A", "KLHL3", "PITX2", "FAM241A", "CELSR2")

  # source = https://cvd.hugeamp.org/phenotype.html?phenotype=HF
  top_single_variant= c("PITX2", "CDKN2B", "LPA", "BAG3", "SURF6", "PSRC1", "CDKN1A", "SYNPO2L")
  top_common_variant= read_csv(file= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/genesets/gene_table.csv")


  #cardiomyopathies:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5933969/table/ehf212267-tbl-0001/
  set_cardiomyopathy= c("LMNA", "MYH7", "TTN", "TNNT2",
                        "MYBPC3", "TNNI3", "TPM1", "MYL3",
                        "DES",
                        "DSC2", "DSG2", "DSP", "JUP", "PKP2")

  hfpef= c("LOX", "SLC25A24", "MTERF1", "SCLT1", "USP1", "ARHGAP18", "SH3BGRL2", "MTURN", "NEIL3" )
  hfpef= c("MYH7", "MYBPC3", "TCAP", "PTGDS", "BSG", "ALDOA", "EEF2", "IDH2", "CRIP2", "EEF1A2", "PLN", "CYCS", "MYOZ2", "MDH1", "DCN")
  reheat= as_tibble(read.csv("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/METArank_March2020.csv"))
  reheat_up= reheat %>% arrange(desc(mean_t))


  return(list(set_dis= sort(disgenet_hflinks),
              set_phe= sort(unique(phewas_hflinks$gene_name)),
              #set_natureGWAS= naturecomm_GWAS,
              set_keggDCM= kegg_DCM,
              set_reheat = reheat$gene[1:top_reheat],
              set_reheat_up = reheat_up$gene[1:top_reheat],
              set_CM= set_cardiomyopathy,
              set_top_sing= top_single_variant,
              set_top_comm_var= top_common_variant$gene[1:100])
  )
}


#sets2= list("set_dis"= c(""))

validate_results= function(sets, gene_results, col, node= "hf"){
  require(pROC)
  library(PRROC)

  x= map(sets, function(set){

    roc.df= gene_results %>%
      arrange(desc(value)) %>%
      mutate(validation_genes = ifelse(name %in% set, 1, 0))

    roc.res= roc(response= roc.df$validation_genes, predictor = roc.df$value)

  })


  plot.df = data.frame(sensitivity = c(x$set_dis$sensitivities,
                                       x$set_phe$sensitivities),
                       specificity = c(x$set_dis$specificities,
                                       x$set_phe$specificities),
                       validation.set= c(rep("set_dis", length(x$set_dis$sensitivities)),
                                         rep("set_phe", length(x$set_phe$sensitivities))
                       ),
                       auc = c(rep(x$set_dis$auc, length(x$set_dis$sensitivities)),
                               rep(x$set_phe$auc, length(x$set_phe$sensitivities))
                       )
  )


  label_auc= paste("AUROC_dis = ",round(x$set_dis$auc,2),"\n",  "AUROC_phe = ",  round(x$set_phe$auc,2))
  roc.plot= plot.df %>% ggplot(., aes(x= specificity, y= sensitivity, color = validation.set))+
    geom_path()+
    geom_abline(slope = 1, intercept = 1)+
    theme_minimal()+
    geom_text(x= -0.45,  y= 0.1, col= "black" , label = label_auc)+
    ggtitle(label = paste("validation gene sets in :", node))+
    scale_x_reverse()

  return(list(plot= roc.plot,
              auc= list("dis"= round(x$set_dis$auc,3),
                        "phe"= round(x$set_phe$auc,3)
              )
  )
  )
}
# gene_results= g.hf2
# set = unlist(sets)
#this function returns AUROC, AUPRC, as well as enrichment results:
validate_results2= function(set, gene_results){
  library(fgsea)
  library(PRROC)

  #gene_results= gene_results %>% filter(value >0)


  roc.df= gene_results %>%
    arrange(desc(value)) %>%
    mutate(validation_genes = ifelse(name %in% set, 1, 0))

  coverage= prop.table(table(unique(set) %in% gene_results$name))

  message(paste0("your validation set coverage is ",
                 round(coverage["TRUE"]*100,2),
                 "%"
                 )
          )



  pr= pr.curve(scores.class0 = roc.df%>% filter(validation_genes== 1) %>% pull(value),
               scores.class1 = roc.df%>% filter(validation_genes== 0) %>% pull(value),
               sorted = T,
               curve =T)

  roc= roc.curve(scores.class0 = roc.df%>% filter(validation_genes== 1) %>% pull(value),
                 scores.class1 = roc.df%>% filter(validation_genes== 0) %>% pull(value),
                 curve =T)

  # stats= gene_results$value
  # names(stats)= gene_results$name
  # fg.res= fgsea::fgseaSimple(list(set= set), stats, nperm= 500, scoreType = "pos")
  # p.fgsea= plotEnrichment(set,stats,gseaParam = 1)

  return(list(roc= roc,
              pr= pr
              )
         )


}

## partial AUROC
validate_results_partial= function(set, gene_results, FPR_cutoff= 0.015){
  library(pROC)
  library(PRROC)



  # get the roc data frame
  roc.df= gene_results %>%
    arrange(desc(value)) %>%
    mutate(validation_genes = ifelse(name %in% set, 1, 0)) %>%
    mutate(rank= rank(dplyr::desc(value)))

  print("number of genes with FPR_cutoff")
  print(dim(roc.df)[1] *FPR_cutoff)


  #### pAUROC
  roc4 <- roc(roc.df$validation_genes,
              roc.df$value,
              percent=F,
              # arguments for auc
              partial.auc=c(0.0,FPR_cutoff),
              partial.auc.correct=FALSE,
              partial.auc.focus="sens",
              # arguments for ci
              ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
              # arguments for plot
              plot=FALSE,
              auc.polygon=TRUE,
              max.auc.polygon=TRUE,
              grid=TRUE,
              print.auc=TRUE, show.thres=TRUE)

  pAUROC= as.numeric(roc4$auc)

  #### pAUROC with corrected
  roc.c <- roc(roc.df$validation_genes,
              roc.df$value,
              percent=F,
              # arguments for auc
              partial.auc=c(0.0,FPR_cutoff),
              partial.auc.correct=T,
              partial.auc.focus="sens",
              # arguments for ci
              ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
              # arguments for plot
              plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
              print.auc=TRUE, show.thres=TRUE)

  ### median rank stat:
  medi.rank= roc.df %>%
    filter(validation_genes== 1) %>%
    mutate(med.rank= median(rank))

  median_rank_stat= medi.rank$med.rank[1] / max(roc.df$rank)

  ## pR.

  pr= pr.curve(scores.class0 = roc.df%>% filter(validation_genes== 1) %>% pull(value),
               scores.class1 = roc.df%>% filter(validation_genes== 0) %>% pull(value),
               curve =T)

  ## full AUROC

  roc= roc.curve(scores.class0 = roc.df%>% filter(validation_genes== 1) %>% pull(value),
                 scores.class1 = roc.df%>% filter(validation_genes== 0) %>% pull(value),
                 curve =T)

  ## coverage of the validation set:

  coverage= prop.table(table(unique(set) %in% gene_results$name))

  message(paste0("your validation set coverage is ",
                 round(coverage["TRUE"]*100,2),
                 "%"
  )
  )


  return(list(pAUROC= pAUROC,
              pAUROC.object= roc4,
              AUROC = roc$auc,
              AUROC.object= roc,
              pr= pr,
              median_rank_stat= median_rank_stat,
              pAUROC.object.corrected= roc.c
  )
  )


}

## plot different AUCs for a given genesets in relation to cut offs used to create networks

plot_cutoff_robustness_multiRW= function(
                                 gamma.seq= seq(0.1,0.9,0.1),
                                 weights= seq(0.0, 0.55, 0.05),
                                 set,
                                 ...){

  auc_r= map(gamma.seq, function(g){
    print(g)
    g.test= run_multi_RW( r =g,
                             weight_dg = 0.1,
                             weight_gg= 0.1,
                         weight_dd= 0.1,
                         prefiltered = T,
                          ...)
    g.hf.val.t= validate_results2(set, g.test)
    return(g.hf.val.t$roc$auc)

  })

  auc_dg= map(weights, function(wdg){
    print(wdg)
    g.test= run_multi_RW( r =0.7,
                               weight_dg = wdg,
                               weight_gg= 0.1,
                               weight_dd= 0.1,
                               prefiltered = T,
                               ...)

    g.hf.val.t= validate_results2(set, g.test)
    return(g.hf.val.t$roc$auc)

  })

  auc_gg= map(weights, function(wdg){
    print(wdg)
    g.test= run_multi_RW( r =0.7,
                          weight_dg = 0.1,
                          weight_gg= wdg,
                          weight_dd= 0.1,
                          prefiltered = T,
                          ...)

    g.hf.val.t= validate_results2(set, g.test)
    return(g.hf.val.t$roc$auc)

  })

  weights2=  seq(0, 0.3, 0.05)
  auc_dd= map(weights2, function(wdg){
    print(wdg)
    g.test= run_multi_RW( r =0.7,
                          weight_dg = 0.1,
                          weight_gg= 0.1,
                          weight_dd= wdg,
                          prefiltered = T,
                          ...)

    g.hf.val.t= validate_results2(set, g.test)
    return(g.hf.val.t$roc$auc)

  })

  res.aucs= list(auc_r, auc_dg, auc_gg, auc_dd)

  plot.aucs= list(
    enframe(unlist(res.aucs[[1]]), name = "r") %>% mutate(r= gamma.seq),
    enframe(unlist(res.aucs[[2]]), name = "dg") %>% mutate(dg= weights),
    enframe(unlist(res.aucs[[3]]), name = "gg") %>% mutate(gg= weights),
    enframe(unlist(res.aucs[[4]]), name = "dd") %>% mutate(dd= weights2)
  )

  names(plot.aucs) = c("r", "dg", "gg", "dd")
  return(plot.aucs)
}


  plot_cutoff_robustness= function(.net,
                                 diseaseoi,
                                 gamma.seq= seq(0.1,0.9,0.1),
                                 weights= seq(0.0, 0.6, 0.05),
                                 set,
                                 ...){

  aucs= map(gamma.seq, function(g){
    print(g)
    g.hf2= get_genes_for_net(.net= .net,
                             gamma =g,
                             eps= 0.001,
                             diseaseoi = diseaseoi,
                             weighted = T,
                             weight_dg = 0.3,
                             weight_gg= 0.2,
                             ...)
    g.hf.val.t= validate_results2(set, g.hf2)
    return(g.hf.val.t$roc$auc)

  })




  #plot( gamma.seq,unlist(aucs))

  #dg weights

  auc_dg= map(weights, function(dg){
    print(dg)
    g.hf2= get_genes_for_net(.net= .net,
                             gamma =0.7,
                             eps= 0.001,
                             diseaseoi = diseaseoi,
                             weighted = T,
                             weight_dg =  dg,
                             ...)

    g.hf.val.t= validate_results2(set, g.hf2)
    return(g.hf.val.t$roc$auc)

  })

  #plot(weights,unlist(auc_dg))

  # gg weights

  auc_gg= map(weights, function(gg){
    print(gg)
    g.hf2= get_genes_for_net(.net= .net,
                             gamma =0.2,
                             eps= 0.001,
                             diseaseoi = diseaseoi,
                             weighted = T,
                             weight_dg =0.20 ,
                             weight_gg = gg,
                             weight_dd= 0.1,
                             ...)
    g.hf.val.t= validate_results2(set, g.hf2)
    return(g.hf.val.t$roc$auc)

  })

  # plot( weights,unlist(auc_gg))


  # dd weights

  weights2=  seq(0, 0.3, 0.05)
  auc_dd= map(weights2, function(gg){
    print(gg)
    g.hf2= get_genes_for_net(.net= .net,
                             .gamma =0.7,
                             eps= 0.001,
                             diseaseoi = diseaseoi,
                             weighted = T,
                             weight_dd = gg,
                             ...)
    g.hf.val.t= validate_results2(set, g.hf2)
    return(g.hf.val.t$roc$auc)

  })

  par(mfrow= c(2,2))
  plot( gamma.seq,unlist(aucs))
  plot( weights2,unlist(auc_dd))
  plot( weights,unlist(auc_gg))
  plot(weights,unlist(auc_dg))



}


## process the mouse consortium phenotype mouse data:

process_IMPC= function(){
  impc= read_tsv(file = "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/IMPC_Cardiovascular_System.tsv")
  gene_trans= readRDS(file= "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/gene_translate.rds")
  df= impc %>%
    left_join(gene_trans %>% dplyr::rename(Gene= MGI.symbol)) %>%
    group_by(Phenotype) %>% filter(!is.na(Gene.name))

  gene_list= split(df$Gene.name,df$Phenotype)
  names(gene_list) = str_replace_all(names(gene_list), " ", "_")
  return(gene_list)
}
## pass two generankings and enrich impc sets and return fgsea results as plot
enrichIMPC= function(generank1, generank2, name1, name2){

  impc= process_IMPC()

  gene_results1= generank1 %>% filter(value>0)
  stats= generank1$value
  names(stats)= generank1$name
  fg.res.1= fgsea::fgseaSimple(impc, stats, nperm= 1000, minSize = 5)

  gene_results2= generank2 %>% filter(value>0)
  stats2= generank2$value
  names(stats2)= generank2$name
  fg.res.2= fgsea::fgseaSimple(impc, stats2, nperm= 1000, minSize = 5)

  fg.res.1$group = name1
  fg.res.2$group = name2
  fg.comb= rbind(fg.res.1, fg.res.2)

  fg.comb$stars <- cut(fg.comb$pval, breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))

  fg.comb%>% ggplot(., aes(x= pathway, y= group, fill =NES ))+
    geom_tile()+
    scale_fill_gradient(low= "white", high = "red")+
    geom_text(aes(label=stars), color="black", size=5)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle= 45, hjust = 1))

}

# function to extract validation results for plotting
plot.val= function(val.res){
  data.frame(auroc= val.res$roc$auc,
             #auprc= val.res$pr$auc.integral,
             gsea.nes=  val.res$fgsea.res$NES,
             gsea.p = val.res$fgsea.res$pval
  )

}


#function to calculate degree of a set of nodes within a network and perform stats test

analyze_node_degrees= function(nodes,
                               g){
  i.g= graph_from_data_frame(g)
  null= degree(i.g, V(i.g))
  oi= degree(i.g, V(i.g)[V(i.g)$name %in% nodes])
  print(paste0("mean degree net= ",mean(null)))
  print(paste0("mean degree set= ",mean(oi)))
  print(paste0("median degree net= ",median(null)))
  print(paste0("median degree set= ",median(oi)))
  print(pnorm(q= mean(oi), mean = mean(null),sd = sd(null),lower.tail = F
          ))
  #mean(oi)
}

# compare rankings:  ------------------------------------------------------

jaccard_rankings= function(rank1, rank2, top= 180){
  set_hfref= rank1 %>% top_n(top) %>% pull(name)
  set_hfpef= rank2 %>% top_n(top) %>% pull(name)

  length(intersect(set_hfref, set_hfpef)) / length(union(set_hfref, set_hfpef))
}

setdiff_rankings= function(rank1, rank2, top= 200){
  set_hfref= rank1 %>% top_n(top,wt =value ) %>% pull(name)
  set_hfpef= rank2 %>% top_n(top, wt= value) %>% pull(name)
  #setdiff(set_hfref, set_hfpef)
  c(set_hfpef[!set_hfpef %in% set_hfref], set_hfref[!set_hfref %in% set_hfpef])
  #length(intersect(set_hfref, set_hfpef)) / length(union(set_hfref, set_hfpef))
}

get_rankings= function(rank1, rank2, name1= "rank1", name2= "rank2"){

  ranks1= rank1 %>% arrange(desc(value)) %>% mutate(value.r1= rank(-value, ties.method = "random"))
  ranks= rank2 %>%
    mutate(value.r2= rank(-value, ties.method = "random")) %>%
    full_join(ranks1,by = "name" ) %>%
    mutate(rank.diff= value.r2-value.r1) %>%
    arrange(desc(abs(rank.diff)))
}
