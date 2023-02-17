## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2023-02-06
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2023
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## assess risk factor distribution for revision
## ---------------------------

library(tidyverse)

directory= "T:/fsa04/MED2-HF-Comorbidities/data/processed_data/sept2022"
risks= readRDS(file.path(directory, "risk.factors.2023.rds"))


pids.list= readRDS("T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/hf_cohort_data/cohort_pids/hf_types_pids2022.rds")
map(pids.list,length)

risk= risks %>%
  filter(pid %in% c(pids.list$hfpef,pids.list$hfmref,  pids.list$hfref))

double.pid= risk %>%
  dplyr::group_by(pid, entry_date, entity) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) %>%
  print(n=100)%>%
  pull(pid)

risk2= risk%>%
  filter(!pid %in% double.pid)%>%
  pivot_wider(names_from = entity, values_from = entry_value)

map(colnames(risk2)[-c(1,2)], function(x){ unique(risk2[[x]])})

glimpse(risk2)

# only use those
ggplot(risk2, aes(x= demografie.risikofaktor.nikotinanamnese))+
  geom_histogram(stat="count")+
  coord_flip()


risk3= risk2 %>%
  rename(smoker=demografie.risikofaktor.nikotinanamnese )%>%
  mutate(smoker= tolower(smoker))%>%
  mutate(smoker= ifelse(is.na(smoker), "unbekannt", smoker),
         smoker= ifelse(smoker=="nie", "nein", smoker) ,
         smoker= ifelse(smoker=="n.a.", "unbekannt", smoker),
         smoker= ifelse(smoker=="n.b.", "unbekannt", smoker),
         smoker= ifelse(smoker=="raucher", "ja", smoker),
         smoker= ifelse(grepl("fr", smoker), "frueher", smoker),
         smoker= ifelse(grepl("gelegentlich", smoker), "gelegentlich", smoker)
         )%>%
  filter(smoker %in% c("unbekannt", "ja", "gelegentlich", "nein","frueher"))%>%
  mutate(hf.type = ifelse(pid %in% pids.list$hfref,
                          "hfref",
                          ifelse(pid %in% pids.list$hfpef,
                                 "hfpef",
                                 ifelse(pid %in% pids.list$hfmref,
                                        "hfmref",
                                        "unlabeled")
                          )))

risk3%>%
  group_by(pid)%>%
  summarise(count())

risk4 = risk3 %>% filter(smoker != "unbekannt")

risk4 = risk4 %>%group_by(pid)%>%
  add_count(smoker) %>%
  select(pid, smoker, hf.type, n)%>%
  print(n=200)%>%
  mutate(maxn = max(n))%>%
  filter(n==maxn)%>%
  distinct()

prop.table(table( c(pids.list$hfpef, pids.list$hfref, pids.list$hfmref) %in% risk4$pid))


library(ComplexHeatmap)

x<- prop.table(table(risk4$smoker, risk4$hf.type), margin = 2)
Heatmap(x, name= "% column")


ggplot(risk3, aes(x= smoker, fill = hf.type))+
  geom_histogram(stat="count")+
  coord_flip()

table(unique(risks$pid) %in% pids.list$hfref)


