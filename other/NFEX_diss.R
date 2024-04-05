#Loading libraries
library(tidyverse)
library(janitor)
library(dplyr)
library(tidyr)
library(readr)


# Section I: ----

# cleaning and wrangling
transition_list <- read_csv("data_raw/HILIC_QE_POSITIVE_G5_Dissolve_NFEX_Isotopes1.csv")


clean_areas <- transition_list %>%
  select(filename=`Replicate Name`, compound_name=`Precursor Ion Name`, area=Area) %>%
  mutate(cmpd_type=ifelse(str_detect(compound_name, "A"), "IS", "Non-IS")) %>%
  mutate(area=as.numeric(area)) %>%
  mutate(samp_type=str_extract(filename, "Poo|Blk|Smp|Std")) %>%
  mutate(day=str_extract(filename, "T\\d+")) %>%
  mutate(time=as.numeric(str_remove(day, "T"))) %>%
  mutate(filename=fct_inorder(filename)) %>%
  arrange(compound_name, filename) %>%
  select(filename, compound_name, area, cmpd_type, samp_type, time)

#filename, area, compound_name?

IS_list <- clean_areas %>%
  filter(cmpd_type=="IS") %>%
  select(filename=`filename`, compound_name_IS=`compound_name`, area_IS=area)

#check there are no Non-IS?

write_csv(clean_areas, "intermediates/NFEX_dissolve_clean_areas.csv")
write_csv(IS_list, "intermediates/NFEX_dissolve_IS_list.csv")


# Section II: ----

#BMIS

BMISed_areas <- clean_areas %>%
  filter(cmpd_type=="Non-IS") %>%
  select(filename, compound_name, area) %>%
  left_join(IS_list) %>%
  select(filename, compound_name, area, area_IS) %>%
  group_by(compound_name) %>%
  mutate(bmised_area=(area/area_IS)*mean(area_IS[31:60], na.rm = TRUE)) %>%
  select(filename, compound_name, bmised_area)

#are there any bmised_areas that lowers 40% but its okay if its not, just keep in mind?

write_csv(BMISed_areas, "intermediates/NFEX_dissolve_BMISed_areas.csv")
