#Loading libraries
library(tidyverse)
library(janitor)
library(dplyr)
library(tidyr)
library(readr)


# Section I: ----

# cleaning and wrangling
transition_list <- read_csv("data_raw/HILIC_QE_POSITIVE_G5_Particulate_NFEX_Isotopes1.csv")

clean_areas <- transition_list %>%
  select(filename=`Replicate Name`, compound_name=`Precursor Ion Name`, area=Area) %>%
  mutate(cmpd_type=ifelse(str_detect(compound_name, "A"), "IS", "Non-IS")) %>%
  mutate(area=as.numeric(area)) %>%
  mutate(samp_type=str_extract(filename, "Poo|Blk|Smp|Std")) %>%
  mutate(day=str_extract(filename, "T\\d+")) %>%
  mutate(day=as.numeric(str_remove(day, "T"))) %>%
  mutate(filename=fct_inorder(filename)) %>%
  arrange(compound_name, filename)

#Standard List that only contains Arsenobetaine from Skyline
IS_list <- clean_areas %>%
  filter(cmpd_type=="IS") %>%
  select(filename=`filename`, compound_name_IS=`compound_name`, area_IS=area)

#Reading in the standards sheet gotten from the Ingalls Lab github
stds_list <- read_csv("data_raw/NFEXIngalls_Lab_Standards.csv") %>%
  clean_names() %>%
  select("compound_name", "hilic_mix", "concentration_u_m", "z") %>%
  filter(z == 1) %>%
  mutate(cmpd_type=ifelse(str_detect(compound_name, ", "), "IS", "Non-IS")) %>%
  semi_join(clean_areas, by = "compound_name") %>%
  filter (cmpd_type == "Non-IS") %>%
  rename(mix = hilic_mix, conc_um = concentration_u_m)



write_csv(clean_areas, "intermediates/NFEX_particulate_clean_areas.csv")
write_csv(IS_list, "intermediates/NFEX_particulate_IS_list.csv")
