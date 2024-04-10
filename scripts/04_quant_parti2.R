# Section III: Quantifying peak area to umol/vial when possible.

#we use actually use the raw results
clean_areas <- read_csv("intermediates/NFEX_parti2_clean_areas.csv")
standards_list <- read_csv("intermediates/NFEX_parti_standards_list.csv")
cmpd_cncs <- read_csv("intermediates/NFEX_parti2_cmpd_cncs.csv")
BMISed_areas <- read_csv("intermediates/NFEX_parti2_BMISed_areas.csv")

# Setting the stage with clean_areas and the concs that were calculate
updated_clean_areas <- clean_areas %>%
  left_join(cmpd_cncs, by = c("filename", "compound_name")) %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type", "day","conc_um")

#wrangling and cleaning now that it has been combined
cleanareas_and_stdslist <- updated_clean_areas %>%
  left_join(standards_list, by = "compound_name") %>%
  filter(compound_name %in% c("Glycine betaine 0(15N), 0(13C)", "Guanine 0(15N), 0(13C)", "Homarine 0(15N, 0(13C)", "L- Glutamine 0(15N), 0(13C)", "L-Glutamic acid 0(15N), 0(13C)")) %>%
  mutate(conc_um = coalesce(conc_um.x, conc_um.y)) %>%
  mutate(cmpd_type = coalesce(cmpd_type.x, cmpd_type.y)) %>%
  select(-conc_um.x, -conc_um.y, -cmpd_type.x, -cmpd_type.y, -z) %>%
  mutate(mix = case_when(
    compound_name %in% c("Glycine betaine 0(15N), 0(13C)", "Guanine 0(15N), 0(13C)", "Homarine 0(15N, 0(13C)") ~ "Mix2",
    compound_name %in% c("L- Glutamine 0(15N), 0(13C)", "L-Glutamic acid 0(15N), 0(13C)") ~ "Mix1",
    TRUE ~ "NoMix"
  )) %>%
  select("filename", "compound_name", "area", "mix", "conc_um")


#now start to calculate the RFs
all_rfs <-  cleanareas_and_stdslist  %>%
  filter(str_detect(filename, "InMatrix|H2OinMatrix")) %>%
  filter(str_detect(filename, paste0(mix, "InMatrix|H2OinMatrix"))) %>%
  mutate(std_type = str_extract(filename, "Mix|H2O")) %>%
  select(-mix) %>%
  group_by(compound_name, std_type, conc_um) %>%
  summarise(mean_area=mean(area)) %>%
  pivot_wider(names_from = std_type, values_from = mean_area) %>%
  mutate(area_from_4um_addition=Mix-H2O) %>%
  mutate(rf=area_from_4um_addition/conc_um) %>%
  mutate(rf_ratio = Mix / H2O) %>%
  select(compound_name, rf, rf_ratio)


#Rf and rf ratios are being combine to all the labeled ones too to then
final_concs <- BMISed_areas %>%
  filter(!str_detect(compound_name, "Std")) %>%
  left_join(all_rfs) %>%
  mutate(N = as.integer(str_extract(compound_name, "\\d+(?=\\(15N\\))")),
         C = as.integer(str_extract(compound_name, "\\d+(?=\\(13C\\))"))) %>%
  group_by(filename) %>%
  fill(rf, rf_ratio, .direction = "down") %>%
  ungroup() %>%
  mutate(uM_in_vial=bmised_area/rf) %>%
  select(compound_name, filename, uM_in_vial) %>%
  mutate(uM_in_env=uM_in_vial/5000) %>%
  mutate(nM=uM_in_env*1000) %>%
  mutate(hour=str_extract(filename, "T\\d+")) %>%
  select(compound_name, filename, nM, hour)




write_csv(final_concs, "intermediates/NFEX_parti_final_concs.csv")


