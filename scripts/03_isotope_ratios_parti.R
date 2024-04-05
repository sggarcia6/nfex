# Extra Section to find ratios of compounds in order to find concentrations

all_cmpd_ref_areas <- read_csv("intermediates/NFEX_parti_clean_areas.csv") %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type")

cmpd_list <- c("Glycine betaine", "Guanine", "L- Glutamine")

loop_output <- list()
for(cmpd_i in cmpd_list){
  cmpd_ratios_i <- all_cmpd_ref_areas %>%
    filter(cmpd_type == "Non-IS", grepl(cmpd_i, compound_name)) %>%
    mutate(ratio_to_reference = area / area[compound_name == paste(cmpd_i, "0(15N), 0(13C)")])
  loop_output[[cmpd_i]] <- cmpd_ratios_i
}
do.call(rbind, loop_output) %>%
# rbind(loop_output[[1]], loop_output[[2]], loop_output[[3]])
  mutate(parent_conc_um=case_when(TRUE~4),
         conc_um= ratio_to_reference*parent_conc_um)


loop_output <- sapply(cmpd_list, function(cmpd_i){
  all_cmpd_ref_areas %>%
    filter(cmpd_type == "Non-IS", grepl(cmpd_i, compound_name)) %>%
    mutate(ratio_to_reference = area / area[compound_name == paste(cmpd_i, "0(15N), 0(13C)")])
}, simplify = FALSE) %>%
  bind_rows() %>%
  mutate(parent_conc_um=case_when(TRUE~4),
         conc_um= ratio_to_reference*parent_conc_um)


#GBT ratios

  gbt_cmpd_ratios <- read_csv("intermediates/NFEX_parti_clean_areas.csv") %>%
    select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
    filter(cmpd_type == "Non-IS", grepl("Glycine betaine", compound_name)) %>%
    mutate(ratio_to_reference = area / area[compound_name == "Glycine betaine 0(15N), 0(13C)"]) %>%
    mutate(parent_conc_um=case_when(TRUE~4),
           conc_um= ratio_to_reference*parent_conc_um)

#Guanine ratios

  guanine_cmpd_ratios <- read_csv("intermediates/NFEX_parti_clean_areas.csv") %>%
    select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
    filter(cmpd_type == "Non-IS", grepl("Guanine", compound_name)) %>%
    mutate(ratio_to_reference = area / area[compound_name == "Guanine 0(15N), 0(13C)"]) %>%
    mutate(parent_conc_um=case_when(TRUE~4),
           conc_um= ratio_to_reference*parent_conc_um)

#Homarine ratios

  homarine_cmpd_ratios <- read_csv("intermediates/NFEX_parti_clean_areas.csv") %>%
    select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
    filter(cmpd_type == "Non-IS", grepl("Homarine", compound_name)) %>%
    mutate(ratio_to_reference = area / area[compound_name == "Homarine 0(15N, 0(13C)"]) %>%
    mutate(parent_conc_um=case_when(TRUE~4),
           conc_um= ratio_to_reference*parent_conc_um)


#L-Glutamic acid- ratios

  glutamic_cmpd_ratios <- read_csv("intermediates/NFEX_parti_clean_areas.csv") %>%
    select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
    filter(cmpd_type == "Non-IS", grepl("L-Glutamic acid", compound_name)) %>%
    mutate(ratio_to_reference = area / area[compound_name == "L-Glutamic acid 0(15N), 0(13C)"]) %>%
    mutate(parent_conc_um=case_when(TRUE~4),
           conc_um= ratio_to_reference*parent_conc_um)

#L-Glutamine ratios

  glutamine_cmpd_ratios <- read_csv("intermediates/NFEX_parti_clean_areas.csv") %>%
    select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
    filter(cmpd_type == "Non-IS", grepl("L- Glutamine", compound_name)) %>%
    mutate(ratio_to_reference = area / area[compound_name == "L- Glutamine 0(15N), 0(13C)"]) %>%
    mutate(parent_conc_um=case_when(TRUE~4),
           conc_um= ratio_to_reference*parent_conc_um)

#Bind them all together

  cmpd_cncs <- bind_rows(
    gbt_cmpd_ratios,
    guanine_cmpd_ratios,
    homarine_cmpd_ratios,
    glutamic_cmpd_ratios,
    glutamine_cmpd_ratios) %>%
    select("filename", "compound_name", "conc_um")



write_csv(cmpd_cncs, "intermediates/NFEX_parti_cmpd_cncs.csv")



