## ---------------------------
##
## Purpose of script:
##
## Author: Jan D. Lanzer
##
## Date Created: 2022-09-01
##
## Copyright MIT License Copyright (c) Jan Lanzer, 2022
##
## Email: jan.lanzer@bioquant.uni-heidelberg.de
##
## ---------------------------
##
## Notes:
##
## Identifying correlation features of binary data
## ---------------------------


library(tidyverse)

.links= readRDS( file ="T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_table_hf_cohort_fil.rds")
link.data= readRDS( "T:/fsa04/MED2-HF-Comorbidities/lanzerjd/manuscript/data/networks/comorbidity/link_data.rds")



## asses whether diseases with large incident difference tend to have smaller correlations:

df= disease_frequencies(link.data$pid, link.data$data)%>% select(-rel_freq)

x= link.data$links%>%
  left_join(df %>% rename(disease1= PheCode,
                          dis1.f= freq))%>%
  left_join(df %>% rename(disease2= PheCode,
                          dis2.f= freq))%>%
  mutate(delta.f= abs(dis1.f-dis2.f))


x= x %>%rowwise()%>% mutate(RR= (a*d)/(dis1.f*dis2.f),
                sd.= (1/a) + (1/(dis1.f*dis2.f)) - (1/d) -(1/(d*d)),
                RR.l= RR-2.56*sd.,
                RR.u= RR+2.56*sd.)

hist(x$RR)
ggplot(x, aes(x= RR.l, y= RR.u))+
  geom_point()


1/x$a + (1/(x$dis1.f*x$dis2.f)) - (1/x$d) -(1(x$d*x$d))


ggplot(x%>% filter(fisher.p.adj<0.05), aes(x= delta.f, y= corr.phi))+
   geom_point()

x %>% distinct(dis1.f, dis1_phenotype)%>%
  ggplot(., aes(x= reorder(dis1_phenotype,-dis1.f), y= dis1.f))+
  geom_col()+
  scale_y_log10()

ggplot(x, aes(x= RR, y= delta.f))+
  geom_point()+
  scale_y_log10()

ggplot(x, aes(x= corr.phi, y= delta.f))+
  geom_point()+
  scale_y_log10()

ggplot(x, aes(y= delta.f, x= corr.tet))+
  geom_point()+
  scale_y_log10()+
  geom_density_2d_filled(alpha = 0.5)

ggplot(x, aes(y= (RR), x= corr.phi))+
  geom_point()+
  scale_y_log10()

ggplot(x, aes(x= (RR), y= log10(odds.ratio)))+
  geom_point()


p1= ggplot(.links%>%
             mutate(edge= ifelse((corr.tet>0 & fisher.p.adj<0.05), "y", "n")), aes(x= pcorr.corpor, y= corr.tet))+
  geom_point(size= 0.1)

p1+ geom_density_2d_filled(alpha = 0.5)+
  geom_vline(xintercept = 0.01)+
  geom_hline(yintercept = 0.1)




# check whether hfepf and hfref patients difference, who is sicker? -------

pids.list$hfpef

data %>%mutate(hf = ifelse(pid   %in% pids.list$hfref,
                        "HFrEF",
                        ifelse(pid   %in% pids.list$hfpef,
                               "HFpEF",
                               ifelse(pid   %in% pids.list$hfmref,
                                      "HFmrEF",
                                      "unlabeled"))))%>%
  filter(pid %in% pids.list$hf_all)%>%
  distinct(entry_date,hf,  pid)%>%
  group_by(hf,pid) %>% count()%>%
  ggplot(., aes(x= hf, y= n))+
  geom_jitter()+
  geom_boxplot()+
  scale_y_log10()

##plot n per visit per patient per group
data.mod = data %>%
  mutate(hf = ifelse(pid   %in% pids.list$hfref,
                           "HFrEF",
                           ifelse(pid   %in% pids.list$hfpef,
                                  "HFpEF",
                                  ifelse(pid   %in% pids.list$hfmref,
                                         "HFmrEF",
                                         "unlabeled"))))

p1= data.mod %>%
  filter(pid %in% pids.list$hf_all)%>%
  distinct(entry_date,hf,  pid)%>%
  group_by(hf,pid) %>% count()%>%
  ggplot(., aes(x= hf, y= n))+
  geom_jitter()+
  geom_boxplot()+
  scale_y_log10()+
  labs(y= "n visits")

##plot n disease per patient per group

data.mod = data.mod %>% left_join(hf_meta, by = "pid")

p2=data.mod %>%
  filter(pid %in% pids.list$hf_all)%>%
  distinct(icd10gm_count,hf,  pid)%>%
  ggplot(., aes(x= hf, y= icd10gm_count))+
  geom_jitter()+
  geom_boxplot()+
  scale_y_log10()


p3= data.mod %>%
  filter(pid %in% pids.list$hf_all)%>%
  distinct(ops_count ,hf,  pid)%>%
  ggplot(., aes(x= hf, y= ops_count ))+
  geom_jitter()+
  geom_boxplot()+
  scale_y_log10()


## scale dis count per visit


  data.mod.sc= data.mod %>%
  filter(pid %in% pids.list$hf_all)%>%
  distinct(entry_date,hf,  pid)%>%
  group_by(hf,pid) %>%
  count() %>%
  left_join(hf_meta)%>%
  mutate(icd10.p.visit = icd10gm_count/n,
         ops.p.visit= ops_count/n)

p5= data.mod.sc %>%
  ggplot(., aes(x= hf, y= icd10.p.visit))+
  geom_jitter()+
  geom_boxplot()+
  scale_y_log10()

p4=data.mod.sc %>%
  ggplot(., aes(x= hf, y= ops.p.visit))+
  geom_jitter()+
  geom_boxplot()+
  scale_y_log10()


cowplot::plot_grid(p1, p2, p3, p4,p5)


# test run scale phi ------------------------------------------------------


