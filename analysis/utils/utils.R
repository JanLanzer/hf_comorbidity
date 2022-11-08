##UTILS



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

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector= col_vector[c(-4,-12)]
hfpef_cols= c("#EF946C", "#785474")
up_dn_cols= c("#0353A4", "#C41E3D")
#function to characterize a given patient cohort:


get_summary_table= function(
  data,
  pids.list){

  require(comorbidity)
  require(cowplot)
  require(ggpubr)
  require(lubridate)

  ########################### DF joining ###########################

  ## create data frame with all additional information,

  # filter for pids of interest:
  df= data# %>% dplyr::filter(pid %in% pids)

  ## add cohort category
  df$patient_cohort= "none"

  for(i in names(pids.list)){
    df = df %>% mutate(patient_cohort= ifelse(pid %in% pids.list[[i]],
                                              i, patient_cohort))
  }

  df = df %>% filter(pid != "none")

  #add hfpef label
  # if(!exists("pids.list")){
  #   pids.list = pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/pidslist_v2021.rds")
  # }
  pid.df= read.csv( "T:/fsa04/MED2-HF-Comorbidities/data/RWH_September2022/raw/levinson_comorbidities_pids_2022-10-06.csv", sep = ";") %>%
      as_tibble


  #add age calculation
  pid.df= pid.df  %>%
    mutate(birthday= lubridate::as_date(birthday),
           entry_date = lubridate::as_date(diag_date_min ),
           ICDint = lubridate::interval(birthday, entry_date),
           age.at.icd = ICDint/dyears(1)) %>%
    select(pid, sex, age.at.icd)

  # clean sex
  pid.df= pid.df%>%
    mutate(sex= ifelse(sex=="m, f", NA,sex ),
           sex= ifelse(grepl("m", sex), "m", sex),
           sex= ifelse(grepl("f", sex), "f", sex))


  df = df %>% select(-sex) %>% left_join(pid.df,by= "pid")

  #phenotypic data
  if(!exists("pheno_df")){
    pheno_df= readRDS("T://fsa04/MED2-HF-Comorbidities/data/processed_data/full_clinic_info.rds")
  }

  df= df %>%
    left_join(pheno_df %>%
                select(-sex)%>%
                dplyr::rename(pid= patient_id),by= "pid")

  HF_allcause= c("I11.0", "I13.0", "I13.2", "I25.5", "I42.0", "I42.5", "I42.8", "I42.9", "I50.0", "I50.1", "I50.9")

  df= df %>% filter(!is.na(PheCode))

  ## calculate comorbidity indexes:

  charlson = comorbidity::comorbidity(df, id= "pid", code ="entry_value", score= "charlson", assign0 = F)
  elixhauser = comorbidity::comorbidity(df, id= "pid", code ="entry_value", score= "elixhauser", assign0 = F)
  df= df %>% left_join(charlson %>% select(pid, wscore, score) %>%
                         rename(charlson_wscore= wscore,
                                charlson_score= score)%>%
                         mutate(pid = as.numeric(pid))) %>%
    left_join(elixhauser %>% select(pid, wscore_ahrq,wscore_vw, score) %>%
                rename(elixhauser_wscore= wscore_vw ,
                       elixhauser_score= score)%>%
                mutate(pid = as.numeric(pid)))

  ## add lipid blood values
  if(!exists("full_lipids")){
    full_lipids= readRDS( file.path("T:/fsa04/MED2-HF-Comorbidities/data/processed_data/sept2022/median_lipids_20221006.rds"))

  }
  df = df %>% left_join(full_lipids)

  if(!exists("hba1c")){
    hba1c= readRDS( file.path("T:/fsa04/MED2-HF-Comorbidities/data/processed_data/sept2022/hba1c_20221006.rds"))

  }

  df= df%>% left_join(hba1c %>%
                        group_by(pid) %>%
                        summarise(median.hba1c = median(entry_value))
                      )
  ##gfr
  if(!exists("gfrcg")){
    gfrcg= readRDS( file.path("T:/fsa04/MED2-HF-Comorbidities/data/processed_data/sept2022/gfrcg_20221006.rds"))

  }

  df= df%>% left_join(gfrcg %>% group_by(pid ) %>% summarise(median.gfrcg = median(entry_value)))

  ## add pharma info:
  if(!exists("pharma")){
    pharma= readRDS( file.path("T:/fsa04/MED2-HF-Comorbidities/data/processed_data/", "Medication_df.rds"))

  }

  df = df %>% left_join(pharma %>% dplyr::rename(pid= patient_id), by= "pid")

  ## add endpoints (HTX, ICD and Intubation) from OPS data:
  if(!exists("ops_endpoint")){
    ops_endpoint = readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_output/pids_endpoints.rds")

  }

  df= df %>% mutate(intu= ifelse(pid %in% ops_endpoint$intubation, "yes", "no"),
                    htx= ifelse(pid %in% ops_endpoint$htx, "yes", "no"),
                    defi= ifelse(pid %in% ops_endpoint$defi, "yes", "no"),
                    pci= ifelse(pid %in% ops_endpoint$pci, "yes", "no"))

  ## we have icd_counts, now we add also phecode counts

  df_phe= df %>% group_by(pid) %>% distinct(pid, PheCode) %>% summarise(PheCode_count= n_distinct(PheCode))

  df= df %>% left_join(df_phe)

  # x= "edd.max"
  #
  # df
  # df.t= df %>% distinct(pid, !!as.symbol(x))
  # missing= sum(is.na(df.t$edd.max))
  # missing_percent= missing/dim(df.t)[1] *100
  # median2 = median(df.t$edd.max,na.rm = T)
  # IQR2= IQR(df.t$edd.max,na.rm = T)
  #
  #
  #library(gt)
  #library(glue)
  #library(gtsummary)
  #df= df %>% mutate(sex= ifelse(sex== "u", NA, sex))
  return(df)
}

plot_summary_table=function(df, col.order= c()){

  table.df= df %>%
    filter(patient_cohort != "none") %>%
    distinct(pid,
             #age_at_HF,
             median.BMI,
             mean.sys,
             mean.dias,
             median.LDL,
             median.HDL,
             median.Chol,
             median.Trigs,
             intu,
             htx,
             defi,
             PheCode_count,
             icd10gm_count,
             charlson_score,
             elixhauser_wscore,
             nyha.max,
             median.bnp,
             patient_cohort,
             edd.max,
             ea.min,
             ea.max,
             ee.max,
             ef.min,
             median.gfrcg,
             median.hba1c,
             age.at.mean,
             sex)

  table.df %>% select(-pid) %>%
    select(patient_cohort,
           sex,
           age.at.mean,
           median.BMI,
           mean.sys,
           mean.dias,
           median.LDL,
           median.HDL,
           median.Trigs,
           median.Chol,
           median.gfrcg,
           median.hba1c,
           median.bnp,
           PheCode_count,
           icd10gm_count,
           charlson_score,
           elixhauser_wscore,
           intu,
           htx,
           defi,
           ef.min,
           edd.max,
           ea.min,
           ea.max,
           ee.max

    )%>%
    tbl_summary(by = patient_cohort) %>%
    add_p()%>%
    add_overall()%>%
    add_n() #%>%
    #modify_header(label ~ "**Variable**") %>%
    #modify_spanning_header(c("stat_1", "stat_2") ~ "**Main Cohort**")

}


#get_summary_table(data, pids.list = pids.hf_ct)


unify_axis= function(x){
  x+
    theme(axis.text = element_text(size= 11, color = "black"),
          axis.title= element_text(size= 9, color= "black"),
          legend.title = element_text(size= 11),
          legend.text = element_text(size= 9))
}
