# Section II: Applying Best-Matched Internal Standard (B-MIS).

#Read out the srcripts from data_cleaning and stuff
clean_areas <- read_csv("intermediates/NFEX_parti2_clean_areas.csv")
IS_list <- read_csv("intermediates/NFEX_parti_IS_list.csv")
standards_list <- read_csv("intermediates/NFEX_parti_standards_list.csv")

#IS_list and standards_list were just used from the first time this was done

#BMIS
BMISed_areas <- clean_areas %>%
  filter(cmpd_type == "Non-IS") %>%
  mutate(filename = gsub(".mzML", "", filename)) %>% #Interesting i had not used this until i found in stackoverflow
  select(filename, compound_name, area) %>%
  left_join(IS_list) %>%
  select(filename, compound_name, area, area_IS) %>%
  group_by(compound_name) %>%
  mutate(bmised_area=(area/area_IS)*mean(area_IS[31:60], na.rm = TRUE)) %>%
  select(filename, compound_name, bmised_area)

#WRite it out

write_csv(BMISed_areas, "intermediates/NFEX_parti2_BMISed_areas.csv")