# this script performs a continguency table representation of diagnostic data
# 1. part : create an intersect matrix
# 2. part : use intersect matrix to calculate the rest of continguency table
# 3. part : use continguency tables to perform fisher's exact test, along  with other stats (RR, OR..)
# 4. part : perform partial correlation  (code adopted from Cory Brunson, https://bitbucket.org/corybrunson/comorbidity/src/master/R/)

# 5. simple test on small cohort

# libraries and data ------------------------------------------------------

## libs

library(igraph)
library(tidyverse)

## data
#directory= "T:/fsa04/MED2-HF-Comorbidities/data/RWH_March2020/"

##

# function to calculate frequencies for a cohort (pids)
#' @param pids, the patient ids that define the cohort
#' @param icd, the icd code table, where every row represents a diagnosis for a patient at a time
#' @param feature, the feature whose frequency should be calculated (e.g. PheCode, Icd10-code etc.)

disease_frequencies = function(pids, icd, feature= "PheCode"){

  #x= as.name(feature)

  # filter data to contain only desired patients and the selected diseas classifier:
  icd_pids = icd %>%
    dplyr::filter(pid %in% pids) %>%
    dplyr::select(pid, {{feature}})%>%
    drop_na

  # remove duplicated rows (count a patient only once)
  icd_pids= icd_pids[!duplicated(icd_pids),]



  # calculate frequencies with table-function (this is only based on the icd3 code)
  disease_freq= as_tibble(as.data.frame(sort(table(icd_pids %>% pull({{feature}})))))
  colnames(disease_freq) = c(feature,"freq")

  #     add relative frequency.
  disease_freq= disease_freq %>%
    mutate(rel_freq = freq / length(unique(pids))) %>%
    arrange(desc(freq))
}

#function to quickly check for "not in"
'%ni%' <- Negate('%in%')

cols.nice= c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51")

# PART 1 ------------------------------------------------------------------

## funciton that will return a matrix where a cell indicates the frequency of the occurence
## of disease i (row) and disease j (column) in a specific cohort
#' @param pid, a vector of patient ids that define the cohort
#' @param data, a table containing each diagnostic entry for a patient. Each row is one unique entry and is defined by time and a pid.
#' @param phecodes, the disease codes of interest as phecodes
#' @return union matrix A

get.intersect= function(pids,phecodes,data){
## prep
  # filter dataset for pids.
  data= data %>%
    filter(pid %in% pids,
           PheCode %in% phecodes)

  # sort phecodes
  phecodes= sort(phecodes)

  data %>% distinct(pid, PheCode) %>% group_by(pid, PheCode) %>% count(pid)

  #create empty union matrix
  A = data.frame(matrix(ncol= length(phecodes), nrow = length(phecodes)))
  colnames(A) = phecodes
  rownames(A) = phecodes

## calculate A
  # group data by patients and count every disease only once per patient
  data_grouped = data %>%
    group_by(pid) %>%
    distinct(pid,PheCode) %>%
    drop_na

  # loop over every diagnosis i and calculate its interesections with every other disease
  # and store information in matrix A
  for (i in rownames(A)){
    #create that contains all patients with diagnosis i and their disease codes
    cop = data_grouped %>%
      filter(any(PheCode==i))  %>%
      #filter(PheCode != i) %>%
      arrange(pid)

    #count how many times each disease appears (in patients with disease i)
    counts = enframe(table(cop$PheCode)) %>%
      mutate(value = as.numeric(value))

    # create the row of A for disease i
    # counts does not contain all diseases, only those with an intersect
    # therfore add the missing diseases with the value 0
    missing = colnames(A)[colnames(A) %ni% counts$name]

    if(length(missing)>0){
      Z= data.frame(name=missing, value =0)
      y = rbind(counts, Z)
    }else{
      y= counts
    }

    y= y %>%
      arrange(name)

    # this condition was created to make sure the order is the same when adding the vector for disease i to the matrix
    if (all(y$name == colnames(A))){
      A[i,]= y$value
      }else{
      break()
    }
  }

  if(isSymmetric(as.matrix(A))){
    message("Intersect Matrix is symmetric!")
      }else{
        message("careful, intersect matrix is not symmetric. Debug please")
  }


  return(A)

}


# PART2 ------------------------------------------------------
## function that will return a list of continguency tables.
#' @param  A, intersect matrix from part 1
#' @param pid, a vector of patient ids that define the cohort
#' @param data, a table containing each diagnostic entry for a patient. Each row is one unique entry and is defined by time and a pid.
#' @param phecodes, the disease codes of interest as phecodes (same as rows and cols in A)

get.continguency.columns= function(A, pids, phecodes, data, disease_ontology){

  ##prep

  # deleting half of A, as A is symmetric
  # This saves computing time, as every disease pair is only calculated once
  if(isSymmetric(as.matrix(A))){
    A[lower.tri(A)] <- NA
  }else{
    message("why is your intersect matrix not Symmetric?")
  }

  # tidy format of A
  A.tidy= A %>% rownames_to_column("disease1") %>%
    as_tibble()%>%
    gather(-disease1, key= disease2, value= a) %>%
    drop_na() %>%
    filter(a >0) %>%
    arrange(disease1)

  # calculate total number of patients for every disease
  total_freq = disease_frequencies(pids,data, feature = disease_ontology) %>%
    mutate(PheCode = as.character(PheCode))

  ## provide measures for continguency table

  # A.tidy contains every pairwise combination of all diseases.
  # This intersect is called a, now we add values for the other 4 cells of the continguency table,
  # based on the frequency count

  disease1= total_freq %>% rename(disease1= PheCode)
  disease2= total_freq %>% rename(disease2= PheCode)
  N= length(pids) #total number of patients

  # a= patients with disease 1 and 2
  # b= patient with disease 1 but not disease 2,
  # c= patients with disease 2 but not 1,
  # d= patients without disease 1 or 2
  A.tidy2= A.tidy %>%
    filter(disease1 != disease2) %>%
    left_join(disease1, by= "disease1") %>%
    rename(b= freq) %>%
    mutate(b= b-a) %>%    # total number of disease 1 minus the union(a) yields number of patients with dis 1 but not dis2
    left_join(disease2, by= "disease2") %>%
    rename(c= freq) %>%
    mutate(c= c-a) %>% # total number of disease 2 minus the union(a) yields number of patients with dis 2 but not dis1
    mutate(d= N-b-c-a) %>% # all patients minus all other cells yields number of patients with neither dis1 nor dis2
    dplyr::select(disease1,disease2,a,b,c,d)

  return(A.tidy2)
}



# PART3 Fisher exact test -------------------------------------------------

## Function that will run a fisher.exact test on every continguency table, and correct pvalues afterwards.
# Output is a df, suitable for network representation.
#' @param A.tidy,

get.pairwise.metric = function(A.tidy){
  require(generics)
  require(psych)
  library(qdapTools)
  require(DescTools)
  ## first part: run fisher exact test on pairwise
  A.res= A.tidy %>%
    mutate(a= a+1,
           b= b+1,
           c= c+1,
           d= d+1)%>%
    rowwise() %>%
    mutate( fisher.p = fisher.test(matrix(c(a,b,c,d), nrow= 2, byrow = T))$p.value,
            odds.ratio = fisher.test(matrix(c(a,b,c,d), nrow= 2, byrow = T))$estimate,
            corr.phi = psych::phi(matrix(c(a,b,c,d), nrow= 2, byrow = T), digits= 5),
            corr.tet = psych::tetrachoric(matrix(c(a,b,c,d), nrow= 2, byrow = T))$rho,
            somersD= ((a*d)- (b*c)) / min((a*d)+(b*c)+ (a*b)+ (c*d), (a*d)+(b*c)+((a*c)+ (b*d)))
           )%>%
    ungroup()%>%
    mutate( fisher.p.adj= p.adjust(fisher.p, method = "fdr"))

    return(A.res)

}


# PART 4 adding partial correlation ---------------------------------------

## This function will calculate the multivariate partial correlation.

get.part.corr = function(pids, phecodes, data){
  #1. part -
  #partial correlation will not be calculated based on the continguency table as this is a multi-variate approach
  #Thus, transform tidy data to a modelling matrix with each row a patient and each colmn a disease

  # filter data for desired patients and phecodes
  data = data %>%
    filter(pid %in% pids,
           PheCode %in% phecodes)%>%
    distinct(pid, PheCode) # important to select only a unique disease entry per patient

  # perform data transformation to disease x disease matrix

  df= split(as.character(data$PheCode), data$pid) %>%
    lapply(function(x) strsplit(x, "-")) %>%
    mtabulate()%>%
    matrix2df("pid")
  df= as.data.frame(df) %>%
    column_to_rownames("pid")

  #2.part

  ## convert to matrix format
  df= df %>%
    mutate_all(as.numeric) %>%
    as.matrix

  #calculate tet correlation, based on a ML estimation
  #t1= sirt::tetrachoric2(df)


  t1 <- sirt::tetrachoric2(dat = df, maxit= 1000,
                            cor.smooth = F)
  ?sirt::tetrachoric2()

  #t1$rho[ "394.2", "395.1"]=t1$rho["395.1","394.2"]= 0.8
  table(rowSums(is.na(t1$rho)) <1,colSums(is.na(t1$rho)) <1)
  #remove NAs , if any
  m= t1$rho[rowSums(is.na(t1$rho)) <1,colSums(is.na(t1$rho)) <1 ]

  eigs <- eigen(m)

  # if (min(eigs$values) < median(abs(m)) * 2^(-6)) {
  #    warning("Sample correlation matrix is far from positive definite.")
  # }

  get_confint= function(t,n = length(pids), z= 3){
    v= dim(t)[1] #number of variables (diseases)

    rho_partial_Z <- 0.5 * log((1 + t) / (1 - t))
    se_p <- matrix(1 / sqrt(n - 3 - (v - 2)), nrow = v, ncol = v)
    rho_partial_Z_lower = rho_partial_Z - z / sqrt(n - 3 - (v - 2))
    rho_partial_Z_upper = rho_partial_Z + z / sqrt(n - 3 - (v - 2))

    lower= as.data.frame(as.matrix(rho_partial_Z_lower)) %>%
      rownames_to_column("disease1") %>%
      gather("disease2", "rho_part_lower", -disease1) %>%
      as_tibble

    upper= as.data.frame(rho_partial_Z_upper) %>%
      rownames_to_column("disease1") %>%
      gather("disease2", "rho_part_upper", -disease1) %>%
      as_tibble

    left_join(lower, upper, by= c("disease1", "disease2"))
  }

  # calculate partial corrs:
  # with shrinkage ( dim >> vol)
  t2= corpcor::pcor.shrink(m)

  # without shrinkage (vol >> dim )
  t3 <- psych::partial.r(m)
  t3.s= cor.shrink(t3) #not used


  # transform to link format
  t1_df= as.data.frame(t1$rho) %>%
    rownames_to_column("disease1") %>%
    gather("disease2", "corr.tet.sirt", -disease1) %>%
    as_tibble

  t2_df= as.data.frame(t2[,]) %>%
    rownames_to_column("disease1") %>%
    gather("disease2", "pcorr.corpor", -disease1) %>%
    as_tibble

  t3_df= as.data.frame(as.matrix(t3)) %>%
    rownames_to_column("disease1") %>%
    gather("disease2", "pcorr.psych", -disease1) %>%
    as_tibble

  t3s_df = as.data.frame(as.data.frame(t3.s[,]) ) %>%
    rownames_to_column("disease1") %>%
    gather("disease2", "pcorr.psych.shrink", -disease1) %>%
    as_tibble


  t_df= t1_df %>%
    left_join(t2_df, by= c("disease1", "disease2")) %>%
    left_join(t3_df, by= c("disease1", "disease2")) %>%
    left_join(t3s_df, by= c("disease1", "disease2")) %>%
    left_join(get_confint(as.data.frame(t2[,]))) %>%
    left_join(get_confint(as.data.frame(t1$rho))%>%
                rename(tet_lower= rho_part_lower ,
                       tet_upper= rho_part_upper)) %>%
    left_join(get_confint(as.data.frame(t3))%>%
                rename(rho_part_lower_psych= rho_part_lower ,
                       rho_part_upper_psych= rho_part_upper)) %>% ## add confident intervalls (for pcor by psych pacakge(t3))
    #rowwise() %>%
    filter(!disease1 == disease2) # remove disease pairs
  apply(t_df, 2, range)
  return(t_df)

  }

# PART 5 wrapper func  ------------------------------------------------------------

## this function wraps the above funcitons and outputs a merged link table

create.links = function(pids, phecodes, data, disease_ontology= "PheCode"){

  #dependencies:
  require(sirt)
  require(corpcor)
  require(qdapTools)
  require(generics)
  require(psych)
  require(tidyverse)
  require(igraph)

  A= get.intersect(pids, phecodes, data)
  message("1/5")
  #print(head(A))
  #2
  A.tidy= get.continguency.columns(A, pids, phecodes, data, disease_ontology )
  message("2/5")
  #print(head(continguency.list))
  #3 stats
  links= get.pairwise.metric(A.tidy = A.tidy)
  message("3/5")

  #4 part. corr
  cors= get.part.corr(pids, phecodes, data)
  message("4/5")
  print(head(cors))
  # merge link results
  links= links %>%
    left_join(cors, by= c("disease1", "disease2")) %>%
    left_join(data %>%
                rename(disease1 = PheCode,
                              dis1_phenotype= Phenotype) %>%
                distinct(disease1, dis1_phenotype),
              by= "disease1") %>%
    left_join(data %>%
                rename(disease2 = PheCode,
                       dis2_phenotype= Phenotype) %>%
                distinct(disease2, dis2_phenotype),
              by= "disease2")

  message("5/5")

  return(links)
}
#
#
# pairwise_matrix <- function(d) {
#   vars <- if (is.factor(d$v1)) {
#     levels(d$v1)
#   } else {
#     d %>%
#       select(v1, v2) %>%
#       unlist() %>% unique()
#   }
#   m <- d %>%
#     select(v1, v2, rho) %>%
#     bind_rows(select(., v1 = v2, v2 = v1, rho)) %>%
#     arrange(v2) %>%
#     spread(v2, rho) %>%
#     arrange(v1) %>%
#     select(-v1) %>%
#     as.matrix() %>% unname()
#   diag(m) <- 1
#   rownames(m) <- colnames(m) <- vars
#   m
# }
#
# partial_data <- function(d, z = 1.96, use.shrinkage = FALSE) {
#   n <- unique(d$n)
#   vars <- if (is.factor(d$v1)) {
#     levels(d$v1)
#   } else {
#     d %>%
#       select(v1, v2) %>%
#       unlist() %>% unique()
#   }
#   v <- colnames()
#   # (pairwise) correlation matrix
#   R_s <- d %>%
#     select(v1, v2, rho) %>%
#     bind_rows(select(., v1 = v2, v2 = v1, rho)) %>%
#     arrange(v2) %>%
#     spread(v2, rho) %>%
#     arrange(v1) %>%
#     select(-v1) %>%
#     as.matrix() %>% unname()
#   diag(df) <- 1
#   R_s= m
#   # warn if the sample correlation matrix is not positive-definite
#   eigs <- eigen(R_s)
#   if (min(eigs$values) < median(abs(R_s)) * 2^(-6)) {
#     warning("Sample correlation matrix is far from positive definite.")
#   }
#   # partial correlation matrix
#   if (use.shrinkage) {
#     R_p <- corpcor::cor2pcor(corpcor::cor.shrink(R_s))
#   } else {
#     #R_p <- psych::partial.r(R_s)
#     R_p <- -(solve(R_s))
#     diag(R_p) <- 1 / (1 - psych::smc(R_s))
#     R_p <- cov2cor(R_p)
#   }
#   Z_p <- .5 * log((1 + R_p) / (1 - R_p))
#   dp <- tibble(
#     v1 = factor(rep(vars, times = ncol(R_p)), levels = vars),
#     v2 = factor(rep(vars, each = nrow(R_p)), levels = vars),
#     rho_partial = as.vector(R_p),
#     rho_partial_Z = as.vector(Z_p)
#   )
#   # augment with confidence intervals
#   d <- d %>%
#     left_join(dp, by = c("v1", "v2")) %>%
#     mutate(
#       rho_partial_Z_lower = rho_partial_Z - z / sqrt(n - 3 - (v - 2)),
#       rho_partial_Z_upper = rho_partial_Z + z / sqrt(n - 3 - (v - 2)),
#       rho_partial_lower = (exp(2 * rho_partial_Z_lower) - 1) /
#         (exp(2 * rho_partial_Z_lower) + 1),
#       rho_partial_upper = (exp(2 * rho_partial_Z_upper) - 1) /
#         (exp(2 * rho_partial_Z_upper) + 1)
#     ) %>%
#     select(-contains("_Z_"))
#   # return data
#   d
# }
