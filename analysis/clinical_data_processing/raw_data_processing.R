## loading data from research data warehouse containing requested clinical variables.

library(vctrs, lib.loc = .libPaths()[1])
library(tidyverse, lib.loc = .libPaths()[1])

library(skimr)

directory= "T:/fsa04/MED2-HF-Comorbidities/data/processed_data/"

raw= read_csv(file= "T:/fsa04/MED2-HF-Comorbidities/data/RWH_August_2020/raw_data/levinson_comorbidities_additional_data_2020-11-11.csv")

# get list of unique features in this dataset:
sort(unique(raw$entity))



# preprocessing steps -----------------------------------------------------

# remove time (only keep date)
raw = raw %>% select(-entry_time, -src)

# remove variables with less than 100 entries
low_freq= raw %>% count(entity) %>% filter(n<100)
raw = raw %>% filter(!entity %in% low_freq$entity)

# remove src column?
raw = raw %>% select(-src)

# patientId as numeric
raw =raw %>% mutate(patient_id = as.numeric(patient_id))


# functions ---------------------------------------------------------------
## general func to turn turn the lists into strings

turn_to_vec = function(column){
  x= lapply(lapply(column, FUN = function(x) paste0(x, collapse=" ,")), FUN = noquote) %>%
    as.character(x)  
}


# Medication --------------------------------------------------------------

## remove the listed format and create a df with one entry of a drug for each row:

spre = raw %>% 
  filter(grepl("^therapie", entity)) %>% 
  group_by(patient_id, entry_date) %>%
  pivot_wider(names_from= entity,values_from= entry_value) %>% 
  mutate(same_length= length(unlist(therapie.medikament)) ==length(unlist(therapie.medikament.dosierung)) ) 

# turn the lists into character string, separated by " ," 
spre$therapie.medikament = 
  lapply(lapply(spre$therapie.medikament, FUN = function(x) paste0(x, collapse=" ,")), FUN = noquote)
spre$therapie.medikament.dosierung = 
  lapply(lapply(spre$therapie.medikament.dosierung, FUN = function(x) paste0(x, collapse=" ,")), FUN = noquote)

# mutate as character: 
spre= spre %>% mutate(therapie.medikament = as.character(therapie.medikament),
                therapie.medikament.dosierung= as.character(therapie.medikament.dosierung)) 


# now create a row per entry, separate entries by " ,"
x1=   spre %>% separate_rows(therapie.medikament, sep = " ,")
x2=  spre %>% separate_rows( therapie.medikament.dosierung, sep =" ,")

# these two data frames cannot be joined because they differ in lenght. 
# This means that for a list of medication (e.g. x= a,b,c) an incomplete list of dosage exists (y= a,b) 
# In this case medication dosage information is discarded. 
x4 =x1 %>% 
  filter(same_length == F) %>% 
  mutate(therapie.medikament.dosierung = NA)

# If the length of both lists matches, an orderd assignment  (e.g x[1] <- y[1], x[2] <- y[2] ..)

x3= x1 %>% 
  filter(same_length ==T ) %>% 
  select(-therapie.medikament.dosierung)%>%
  add_column("therapie.medikament.dosierung"= x2 %>% 
               filter(same_length== T) %>% 
               pull(therapie.medikament.dosierung))

#create final dataframe = 
pharma = rbind(x3,x4) %>% select(-same_length) %>% ungroup %>%
  mutate(therapie.medikament.dosierung = tolower(therapie.medikament.dosierung),
         therapie.medikament= tolower(therapie.medikament))

saveRDS(pharma, file = file.path(directory, "Medication_df.rds"))

pharma %>% filter(patient_id %in% pids.list$hf_all) %>% distinct(patient_id, entry_date, therapie.medikament ) %>% count(therapie.medikament) %>%
  arrange(desc(n)) %>% print(n= 100)

# Demografics -------------------------------------------------------------
demo = raw %>% 
  filter(grepl("^demo", entity)) %>% 
  group_by(patient_id, entry_date) %>%
  pivot_wider(names_from= entity,values_from= entry_value)



#### weight, height and BMI
demo$demografie.gewicht = turn_to_vec(demo$demografie.gewicht)

weight =demo %>% select(patient_id, entry_date, demografie.gewicht)%>%
  separate_rows(demografie.gewicht, sep = " ,") %>%
    mutate(demografie.gewicht = as.numeric(demografie.gewicht)) %>% 
  drop_na

#filter spurious weights
weight = weight%>% mutate(demografie.gewicht = ifelse(demografie.gewicht>400, NA,demografie.gewicht))

# check weight change in a patient:
weight %>% 
  group_by(patient_id)%>%
  mutate(ra = abs(min(demografie.gewicht)- max(demografie.gewicht)))%>%
  ggplot(aes(x= ra))+
  geom_histogram()
                                      


#### Gr??e

demo$demografie.koerpergroesse = turn_to_vec(demo$demografie.koerpergroesse)

height =demo %>% separate_rows(demografie.koerpergroesse, sep = " ,") %>%
  mutate(demografie.koerpergroesse = as.numeric(demografie.koerpergroesse)) %>% 
  drop_na

hist(height$demografie.koerpergroesse)

#filter spurious weights
height = height%>%
  mutate(demografie.koerpergroesse = ifelse(demografie.koerpergroesse>250,
                                                             NA,
                                                             demografie.koerpergroesse),
         demografie.koerpergroesse = ifelse(demografie.koerpergroesse<20,
                                            NA,
                                            demografie.koerpergroesse),)

# take the mean of multiple measures (height doesn change that much as opposed to weight

height = height %>% group_by(patient_id)%>%
  summarise(m_height = mean(demografie.koerpergroesse))

# combine height and weight to calculate bmi: 

bmi = weight %>% left_join(height) %>% mutate(bmi = demografie.gewicht/ ((m_height/100)^2) )

bmi %>% filter(bmi> 100) %>% print(n=100)
bmi %>% ggplot(.,aes(x= bmi))+geom_density()


saveRDS(bmi, file = file.path(directory, "processed_data/BMI_df.rds"))

       

# EF ----------------------------------------------------------------------

####echo based EF cleansing: 
ef.echo = raw %>% 
  filter(entity %in% c("diagnostik.echo.lv_ef")) %>% 
  group_by(patient_id, entry_date) %>%
  pivot_wider(names_from= entity,values_from= entry_value)

#transform the lists into a string 
ef.echo$diagnostik.echo.lv_ef = turn_to_vec(ef.echo$diagnostik.echo.lv_ef)

#separate the string into different rows and convert to numeric
ef.echo2= ef.echo %>%
  separate_rows(diagnostik.echo.lv_ef, sep = " ,") %>%
  mutate(diagnostik.echo.lv_ef =as.numeric(gsub(",", ".",diagnostik.echo.lv_ef ) ))
  
## filtering wrong entries
ef.echo2 = ef.echo2 %>% filter(diagnostik.echo.lv_ef<90 & diagnostik.echo.lv_ef> 0)

# plot distribution
ggplot(ef.echo2, aes(x= diagnostik.echo.lv_ef))+
  geom_density()

# looks OK! Rather binned variable, not really continous.

#### echo based MRI cleansing: 
ef.mri= raw %>% 
  filter(entity %in% c("diagnostik.mrt.lv_ef")) %>% 
  group_by(patient_id, entry_date) %>%
  pivot_wider(names_from= entity,values_from= entry_value)

#transform the lists into a string 
ef.mri$diagnostik.mrt.lv_ef = turn_to_vec(ef.mri$diagnostik.mrt.lv_ef)

#separate the string into different rows and convert to numeric
ef.mri2= ef.mri %>%
  separate_rows(diagnostik.mrt.lv_ef, sep = " ,") %>%
  mutate(diagnostik.mrt.lv_ef =as.numeric(gsub(",", ".",diagnostik.mrt.lv_ef ) ))

## filtering wrong entries
ef.mri2 = ef.mri2 %>% filter(diagnostik.mrt.lv_ef<90 & diagnostik.mrt.lv_ef> 0)

# plot distribution
ggplot(ef.mri2, aes(x= diagnostik.mrt.lv_ef))+
  geom_density()


### combination of both data sets to a single EF measurement. 

ef.full = full_join(ef.echo2, ef.mri2)

#around 315 cases wher echo and mri is both apparent on the same day. 
ef.full2= ef.full %>% drop_na %>% mutate(diff= abs(diagnostik.echo.lv_ef- diagnostik.mrt.lv_ef))
#histogram of difference:
hist(ef.full2$diff)

# rules for merging: 
# MRI measurment is more accurate and preferred if both entries are apparaent.otherwise merge

ef.full= ef.full %>% 
  mutate(ef = ifelse(is.na(diagnostik.mrt.lv_ef), diagnostik.echo.lv_ef, diagnostik.mrt.lv_ef )) %>%
  select(-diagnostik.echo.lv_ef, -diagnostik.mrt.lv_ef)

hist(ef.full$ef)

saveRDS(ef.full, file = file.path(directory, "processed_data/EF_df.rds"))
ef= readRDS( file = file.path(directory, "processed_data/EF_df.rds"))




# echo.lv -----------------------------------------------------------------
### LV EDD
#list of echo paramters:
unique(raw$entity)[grep(pattern = "diagnostik.",unique(raw$entity) )]

lv_edd.echo = raw %>% 
  filter(entity %in% c("diagnostik.echo.lv_edd")) %>%
  mutate(entry_value = gsub(",", ".", entry_value), #
         entry_value = gsub("<", "", entry_value), # remove the < > signs and take the margin entry. (which should be an important threshold already)
         entry_value = gsub(">", "", entry_value), 
         entry_value= as.numeric(entry_value)) 

lv_edd.echo.max= lv_edd.echo %>% 
  group_by(patient_id) %>%
  filter(entry_value <200) %>%
  summarise("maxLv.EDD"=max(entry_value)) 

## LV diagnostik.echo.lv_ee_ruhe
lv_ee.echo = raw %>% 
  filter(entity %in% c("diagnostik.echo.lv_ee_ruhe")) %>%
  mutate(entry_value = gsub(",", ".", entry_value), #
         entry_value = gsub("<", "", entry_value), # remove the < > signs and take the margin entry. (which should be an important threshold already)
         entry_value = gsub(">", "", entry_value), 
         entry_value= as.numeric(entry_value)) %>%
  filter(entry_value <100) # rationale

lv_ee.echo.max= 
  lv_ee.echo %>% 
  group_by(patient_id) %>%
  summarise(entry_value=max(entry_value)) #%>%
  ggplot(.,aes(x= entry_value))+ geom_histogram()

## LV diagnostik e/a
lv_ea.echo = raw %>% 
  filter(entity %in% c("diagnostik.echo.lv_ea")) %>%
  mutate(entry_value = gsub(",", ".", entry_value), #
         entry_value = gsub("<", "", entry_value), # remove the < > signs and take the margin entry. (which should be an important threshold already)
         entry_value = gsub(">", "", entry_value), 
         entry_value= as.numeric(entry_value)) %>%
  filter(entry_value <100) # rationale

lv_ea.echo.min= lv_ea.echo %>% 
  group_by(patient_id) %>%
  summarise(entry_value=min(entry_value))#%>%
  ggplot(.,aes(x= entry_value))+ geom_histogram()

lv_ea.echo.min %>% filter(entry_value<2 & entry_value >0.9)

lv_ea.echo.max= lv_ea.echo %>% 
    group_by(patient_id) %>%
    summarise(entry_value=max(entry_value))#%>%
  ggplot(.,aes(x= entry_value))+ geom_histogram()
  

lv_ea.echo.max = lv_ea.echo.max %>% filter(entry_value<25)


#combined echo_df

echo_df= lv_edd.echo.max %>% rename(edd.max= maxLv.EDD)  %>%
  full_join(lv_ee.echo.max %>% rename(ee.max=entry_value))  %>% 
  full_join(lv_ea.echo.min %>% rename(ea.min=entry_value)) %>%
  full_join(lv_ea.echo.max %>% rename(ea.max=entry_value))%>%
  full_join(ef  %>% group_by(patient_id) %>% 
              summarise(ef.min= min(ef)))

echo_df
saveRDS(echo_df, file = file.path(directory, "echo_measures_most_extreme.rds") )

# nyha --------------------------------------------------------------------

nyha= raw %>% 
  filter(grepl("nyha", entity)) %>%
  mutate(entry_value = factor(entry_value, levels = c("I", "II","III", "IV"), ordered = T))

unique(nyha$entry_value) 

nyha.max = nyha %>% 
  group_by(patient_id) %>%
  summarise(max(entry_value))

saveRDS(nyha, file = file.path(directory, "NYHA.rds") )

# Lab value - BNP ---------------------------------------------------------

# 1. transform all values to numeric. 
labs= raw %>% 
  filter(grepl("labor", entity)) %>%
  mutate(entry_value = gsub(",", ".", entry_value), #
         entry_value = gsub("<", "", entry_value), # remove the < > signs and take the margin entry. (which should be an important threshold already)
         entry_value = gsub(">", "", entry_value), 
         entry_value= as.numeric(entry_value)) %>%
  drop_na()

#2. every single variable needs to get cleaned itsel. 
# BNP
labs.bnp = labs %>% filter(grepl("ntbnp", entity)) 

labs.bnp%>% 
  ggplot(.,aes(x= entry_value))+geom_histogram()

labs.bnp %>% 
  filter(entry_value < 40000) %>% 
  ggplot(.,aes(x= entry_value))+geom_histogram()

##excessively high bnp values are not uncommon (10k-30k), a filter is applied at 60k
labs.bnp= labs.bnp %>% 
  filter(entry_value < 100000)

labs.bnp %>% group_by(patient_id) %>%
  summarise(sd(entry_value)) %>% 
  ggplot(.,aes(x= `sd(entry_value)`))+geom_histogram()

pids.maxbnp = labs.bnp %>% group_by(patient_id) %>%
  summarise("max.bnp"= max(entry_value)) %>% 
  filter(max.bnp > 125) %>% 
  pull(patient_id)

saveRDS(labs.bnp, file = file.path(directory, "labs_bnp.rds") )

# blood pressure ----------------------------------------------------------
bp = raw %>% filter(grepl("blut", entity))# %>% pivot_wider(id_cols = c(entity,entry_value))

bp2 = bp %>% 
  mutate(entry_value2= as.numeric(entry_value),
         patient_id = as.numeric(patient_id))

bp2[is.na(bp2$entry_value2),] # only one entry couldnt be transformed to numeric.We will dop it


bp2 = bp2 %>% 
  group_by(patient_id) %>% 
  summarise(count = n()) %>% 
  left_join( bp2 %>%
    group_by(patient_id) %>% 
    filter(entity== "diagnostik.systolischer_blutdruck") %>%
    summarise(mean.sys = median(entry_value2)) )  %>% 
  left_join( bp2 %>%
    group_by(patient_id) %>% 
    filter(entity== "diagnostik.diastolischer_blutdruck") %>%
    summarise(mean.dias = median(entry_value2)))%>%
  drop_na()

bp2= bp2 %>% filter(mean.sys <330, mean.dias< 200)
ggplot(bp2 %>% filter(count<40), aes(x= count))+
  geom_histogram(bins= 100)

bp2 %>% pivot_longer(cols = c(mean.sys, mean.dias)) %>%
  ggplot(aes(x= value, fill = name))+
  #geom_histogram(bins= 100)+
  geom_density(alpha = 0.5)

saveRDS(bp2, file = file.path(directory, "diag_bp.rds") )

range(bp2$mean.sys)




# weight ------------------------------------------------------------------
weight = raw %>% 
  filter(grepl("gewicht", entity))# %>% pivot_wider(id_cols = c(entity,entry_value))

weight= weight %>% 
  mutate(entry_value = str_replace_all(entry_value, ",", "."),
         entry_value = str_replace_all(entry_value, "kg", ""),
         entry_value = str_replace_all(entry_value, "ca", ""),
         entry_value = str_replace_all(entry_value, ". ", ""),
         weight= as.numeric(entry_value)) %>% 
filter(weight<500 & weight > 20) 

weight %>% filter(is.na(weight))

weight %>% 
  #filter(weight<500 & weight > 20) %>%
  ggplot(aes(x= weight))+
  geom_histogram()

range(weight$weight, na.rm = T)

height= raw %>% 
  filter(grepl("groesse", entity))

height= height %>% 
  mutate(entry_value = str_replace_all(entry_value, ",", "."),
         entry_value = str_replace_all(entry_value, "cm", ""),
         entry_value = str_replace_all(entry_value, "q", ""),
         entry_value = str_replace_all(entry_value, "-", ""),
         height= as.numeric(entry_value),
         height= ifelse(height<5, height*100, height)) %>% 
  filter(height != 0 & height <300) ## remove those

height  %>% ggplot(aes(x= height))+
  geom_histogram()

range(height$height, na.rm = T)

small_peops = height %>% filter(height < 100) %>% pull(patient_id)

hf %>% filter(pid %in% small_peops) %>% arrange(desc(birthday))


# many of these people are old and not children.. 
# we will discard heights below 130

height %>% 
  filter(height > 120) %>% 
  group_by(patient_id) %>%
  count(entry_date) %>% 
  ggplot(aes(n)) + geom_histogram()

# usually people only have 1 height




# load lipids and BMI . processed by rebecca ----------------------------------

bmi= as_tibble(readRDS(file.path(directory, "height_weight_pairs_temp_20210621.rds"))) 
tri= as_tibble(readRDS(file.path(directory, "trigs_20210621.rds"))) 
hdl = as_tibble(readRDS(file.path(directory, "HDL_20210621.rds"))) 
ldl = as_tibble(readRDS(file.path(directory, "LDL_20210621.rds"))) 
tc = as_tibble(readRDS(file.path(directory, "Total_Cholesterol_20210621.rds"))) 

df= tri
col = "Trigs"
get_median = function(df, 
                      col){
  #myenc <- enquo(col)
  
  df2= df %>% 
    group_by(pid) %>%
    mutate("median.{col}" := median(!! sym(col), na.rm = T)) 
  #colm = paste0("median.", col)
  #print(colnames(df2))
  df2  %>% 
    select(starts_with("median")) %>% 
    distinct()
}

# get median per patient: 
bmi = get_median(bmi %>% rename(pid= patient_id),"BMI")
tri = get_median(tri, "Trigs")
hdl = get_median(hdl, "HDL")
ldl = get_median(ldl, "LDL")
tc = get_median(tc, "Chol")

# merge to df: 

full_lipids= full_join(bmi, tri) %>% full_join(hdl) %>% full_join(ldl) %>% full_join(tc)
saveRDS(full_lipids, file.path(directory, "median_lipids_20210621.rds"))


