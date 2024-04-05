

read_csv("intermediates/NFEX_parti_clean_areas.csv") %>%
  select("filename", "compound_name", "area", "cmpd_type", "samp_type") %>%
  filter(cmpd_type == "Non-IS", grepl("Glycine betaine", compound_name)) %>%
  mutate(ratio_to_reference = area / area[compound_name == "Glycine betaine 0(15N), 0(13C)"])
