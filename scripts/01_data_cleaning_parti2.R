#Loading libraries
library(tidyverse)
library(janitor)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

#probably when this is added to and Rmd then i should always load libraries bc its becoming a problem

# Section I: ----

# Cleaning and wrangling of skyline output
transition_list <- read_csv("data_raw/NFEX2_parti_rawdata.csv")

clean_areas <- transition_list %>%
  mutate(cmpd_type=ifelse(str_detect(compound_name, "A"), "IS", "Non-IS")) %>%
  mutate(area=as.numeric(area)) %>%
  mutate(samp_type=str_extract(filename, "Poo|Blk|Smp|Std")) %>%
  mutate(day=str_extract(filename, "T\\d+")) %>%
  mutate(day=as.numeric(str_remove(day, "T"))) %>%
  mutate(filename=fct_inorder(filename)) %>%
  arrange(compound_name, filename)


# #Internal standard area from skyline
# IS_list <- clean_areas %>%
#   filter(cmpd_type=="IS") %>%
#   select(filename=`filename`, compound_name_IS=`compound_name`, area_IS=area)
# 
# #Read in Internal Standards list from Ingalls Lab
# standards_list <- read_csv("data_raw/NFEXIngalls_Lab_Standards.csv") %>%
#   clean_names() %>%
#   select("compound_name", "hilic_mix", "concentration_u_m", "z") %>%
#   filter(z == 1) %>%
#   mutate(cmpd_type=ifelse(str_detect(compound_name, ", "), "IS", "IS")) %>% #how can i change it so its always IS?
#   semi_join(clean_areas, by = "compound_name") %>%
#   rename(mix = hilic_mix, conc_um = concentration_u_m)
# 
# #Writing them out
# write_csv(clean_areas, "intermediates/NFEX_diss_clean_areas.csv")
# write_csv(IS_list, "intermediates/NFEX_diss_IS_list.csv")
# write_csv(standards_list, "intermediates/NFEX_diss_standards_list.csv")