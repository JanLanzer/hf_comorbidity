library(dbparser)

## parse data from XML and save it to memory
drug_b= read_drugbank_xml_db("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/drug_bank.xml")


## load drugs data
drugs <- dbparser::drugs()

targets= dbparser::targets()
## load drug groups data
drug_groups <- drug_groups()
dbparser::tar
## load drug targets actions data
drug_targets_actions <- targets_actions()

drugs$pathway
drugs$snp_effects
drugs$drug_classification
drugs$general_information
saveRDS(list("general"= drugs,
             "groups"= drug_groups,
             "target_action" = drug_targets_actions,
             "targets"= targets),
        "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/drug_bank.rds" )

drugbank= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/data_other/Databases/drug_bank.rds" )
drugbank
unique(drug_targets_actions$action)
drugs$interactions

unique(targets$id)
drugs$reactions_enzymes

glimpse(targets)
