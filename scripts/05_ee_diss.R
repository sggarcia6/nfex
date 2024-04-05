# Section IV: Limit of Detection for Dissolved Samples

diss_final_concs <- read_csv("intermediates/NFEX_diss_final_concs.csv")
all_rfs <- read_csv("intermediates/NFEX_diss_all_rfs.csv")


# Load Extraction efficiency list from Josh
ee_list <- read_csv("data_raw/js_diss_ee.csv") %>%
  clean_names() %>%
  filter(compound %in% c("Glycine betaine", "Homarine", "L-glutamic acid", "L-glutamine", "Guanine")) %>%
  mutate(parent_compound_name = case_when(
    compound == "Glycine betaine" ~ "Glycine betaine",
    compound == "Homarine" ~ "Homarine",
    compound == "Guanine" ~ "Guanine",
    compound == "L-glutamine" ~ "L- Glutamine",
    compound == "L-glutamic acid" ~ "L-Glutamic acid")) %>%
  select("parent_compound_name","fraction","ee_percent","rsd_of_ee_percent", "r2","lod_n_m" )

# RFs
all_rfs <- read_csv("intermediates/NFEX_diss_all_rfs.csv") %>%
  mutate(parent_compound_name = str_remove(compound_name, " \\d.*")) %>%
  select("parent_compound_name", "rf", "rf_ratio")


## Determine Limit of Detection

lod_blks <- diss_final_concs %>%
  filter(str_detect(filename, "230831_Smp_LabBlk")) %>%
  group_by(compound_name) %>%
  filter(!is.na(nM)) %>%
  mutate(blank_avg = mean(nM, na.rm = TRUE),
         SB = sd(nM))  %>%
  mutate(parent_compound_name = str_remove(compound_name, " \\d.*")) %>%
  left_join(ee_list, by = "parent_compound_name") %>%
  left_join(all_rfs, by = "parent_compound_name") %>%
  group_by(compound_name) %>%
  slice(1) %>%
  mutate(nis = 15) %>%
  mutate(LOD = (blank_avg + ((1.761 * SB) / sqrt(nis))) * (1 / ee_percent)) %>%
  select(parent_compound_name, compound_name, LOD)





## Extraction Efficiency

#took into account extraction efficiency
diss_final_concs_w_ee <- diss_final_concs %>%
  mutate(parent_compound_name = str_remove(compound_name, " \\d.*")) %>%
  left_join(ee_list, by = "parent_compound_name") %>%
  mutate(final_nM = nM * ee_percent) %>%
  left_join(lod_blks, by = "compound_name") %>%
  mutate(cxc_final_nM = final_nM >= LOD) %>%
  filter(final_nM > LOD) %>%
  select("compound_name", "filename", "final_nM")




## Write out the CSV for all
write_csv(diss_final_concs_w_ee, "intermediates/NFEX_diss_full_final_concs_after_ee.csv")


