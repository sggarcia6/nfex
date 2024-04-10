# Extra Section to find ratios of compounds in order to find concentrations

all_cmpd_ref_areas <- read_csv("intermediates/NFEX_parti2_clean_areas.csv") %>%
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


#carnitine ratios

carnitine_cmpd_ratios <- read_csv("intermediates/NFEX_parti2_clean_areas.csv") %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
  filter(cmpd_type == "Non-IS", grepl("Carnitine", compound_name)) %>%
  mutate(ratio_to_reference = area / area[compound_name == "Carnitine 0(15N), 0(13C)"]) %>%
  mutate(parent_conc_um=case_when(TRUE~1),
         conc_um= ratio_to_reference*parent_conc_um)

#choline ratios

choline_cmpd_ratios <- read_csv("intermediates/NFEX_parti2_clean_areas.csv") %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
  filter(cmpd_type == "Non-IS", grepl("Choline", compound_name)) %>%
  mutate(ratio_to_reference = area / area[compound_name == "Choline 0(15N), 0(13C)"]) %>%
  mutate(parent_conc_um=case_when(TRUE~4),
         conc_um= ratio_to_reference*parent_conc_um)

#dmg ratios

dimethylglycine_cmpd_ratios <- read_csv("intermediates/NFEX_parti2_clean_areas.csv") %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
  filter(cmpd_type == "Non-IS", grepl("Dimethylglycine", compound_name)) %>%
  mutate(ratio_to_reference = area / area[compound_name == "Dimethylglycine 0(15N), 0(13C)"]) %>%
  mutate(parent_conc_um=case_when(TRUE~4),
         conc_um= ratio_to_reference*parent_conc_um)


#ectoine ratios

ectoine_cmpd_ratios <- read_csv("intermediates/NFEX_parti2_clean_areas.csv") %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
  filter(cmpd_type == "Non-IS", grepl("Ectoine", compound_name)) %>%
  mutate(ratio_to_reference = area / area[compound_name == "Ectoine 0(15N), 0(13C)"]) %>%
  mutate(parent_conc_um=case_when(TRUE~8),
         conc_um= ratio_to_reference*parent_conc_um)

#glycine ratios

glycine_cmpd_ratios <- read_csv("intermediates/NFEX_parti2_clean_areas.csv") %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
  filter(cmpd_type == "Non-IS", grepl("Glycine", compound_name)) %>%
  mutate(ratio_to_reference = area / area[compound_name == "Glycine 0(15N), 0(13C)"]) %>%
  mutate(parent_conc_um=case_when(TRUE~4),
         conc_um= ratio_to_reference*parent_conc_um)
#lysine
lysine_cmpd_ratios <- read_csv("intermediates/NFEX_parti2_clean_areas.csv") %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
  filter(cmpd_type == "Non-IS", grepl("Lysine", compound_name)) %>%
  mutate(ratio_to_reference = area / area[compound_name == "Lysine 0(15N), 0(13C)"]) %>%
  mutate(parent_conc_um=case_when(TRUE~8),
         conc_um= ratio_to_reference*parent_conc_um)

#proline ratios

proline_cmpd_ratios <- read_csv("intermediates/NFEX_parti2_clean_areas.csv") %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
  filter(cmpd_type == "Non-IS", grepl("Proline", compound_name)) %>%
  mutate(ratio_to_reference = area / area[compound_name == "Proline 0(15N), 0(13C)"]) %>%
  mutate(parent_conc_um=case_when(TRUE~4),
         conc_um= ratio_to_reference*parent_conc_um)
#proline betaine ratios
prolinebetaine_cmpd_ratios <- read_csv("intermediates/NFEX_parti2_clean_areas.csv") %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
  filter(cmpd_type == "Non-IS", grepl("Proline Betaine", compound_name)) %>%
  mutate(ratio_to_reference = area / area[compound_name == "Proline Betaine 0(15N), 0(13C)"]) %>%
  mutate(parent_conc_um=case_when(TRUE~8),
         conc_um= ratio_to_reference*parent_conc_um)

#TMAO ratios

tmao_cmpd_ratios <- read_csv("intermediates/NFEX_parti2_clean_areas.csv") %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
  filter(cmpd_type == "Non-IS", grepl("Trimethylamine N-Oxide", compound_name)) %>%
  mutate(ratio_to_reference = area / area[compound_name == "Trimethylamine N-Oxide 0(15N), 0(13C)"]) %>%
  mutate(parent_conc_um=case_when(TRUE~4),
         conc_um= ratio_to_reference*parent_conc_um)


#Bind them all together

cmpd_cncs <- bind_rows(
  carnitine_cmpd_ratios,
  choline_cmpd_ratios,
  dimethylglycine_cmpd_ratios,
  ectoine_cmpd_ratios,
  glycine_cmpd_ratios,
  lysine_cmpd_ratios,
  proline_cmpd_ratios,
  prolinebetaine_cmpd_ratios,
  tmao_cmpd_ratios) %>%
  select("filename", "compound_name", "conc_um")



write_csv(cmpd_cncs, "intermediates/NFEX_parti2_cmpd_cncs.csv")