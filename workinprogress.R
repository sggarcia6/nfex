

library(tidyverse)

# Load the data
diss_concs <- read_csv("intermediates/NFEX_diss_full_final_concs_after_ee.csv") %>%
  mutate(
    samp_type = case_when(
      str_detect(filename, "Smp") ~ "Smp",
      str_detect(filename, "Poo") ~ "Poo",
      str_detect(filename, "Blk") ~ "Blk",
      str_detect(filename, "Std") ~ "Std",
      TRUE ~ "Other"
    ),
    hour = str_extract(filename, "T\\d+"),
    hour = as.numeric(str_remove(hour, "T")),
    filename = str_replace(filename, "^[^_]*_(.*)_NFEX", "\\1"),
    diss_nM = final_nM) %>%
  filter(samp_type == "Smp"& !is.na(hour)) %>% 
  select(compound_name, filename, diss_nM)


parti_concs <- read_csv("intermediates/NFEX_parti_final_concs.csv") %>%
  mutate(
    samp_type = case_when(
      str_detect(filename, "Smp") ~ "Smp",
      str_detect(filename, "Poo") ~ "Poo",
      str_detect(filename, "Blk") ~ "Blk",
      str_detect(filename, "Std") ~ "Std",
      TRUE ~ "Other"
    ),
    hour = str_extract(filename, "T\\d+"),
    hour = as.numeric(str_remove(hour, "T")),
    filename = str_replace(filename, "^[^_]*_(.*)_NFEX", "\\1"),
    parti_nM = nM) %>%
  filter(samp_type == "Smp"& !is.na(hour)) %>%
  select(compound_name, filename, parti_nM, hour)


all_concs <- left_join(diss_concs, parti_concs, by = c("filename", "compound_name")) %>%
  mutate(treatment = str_extract(filename, "NO3|NH4|Urea"),
         treatment = if_else(is.na(treatment), "Control", treatment),
         trpl = str_extract(filename, "([ABC])$"))


tbd_all_concs <- left_join(diss_concs, parti_concs, by = c("filename", "compound_name")) %>%
  mutate(
    treatment = str_extract(filename, "NO3|NH4|Urea"),
    treatment = if_else(is.na(treatment), "Control", treatment),
    trpl = str_extract(filename, "([ABC])$")) %>%
  filter(compound_name %in% c("Glycine betaine 0(15N), 0(13C)", 
                              "Homarine 0(15N, 0(13C)", 
                              "L- Glutamine 0(15N), 0(13C)", 
                              "L-Glutamic acid 0(15N), 0(13C)", 
                              "Guanine 0(15N), 0(13C)")) %>%
  filter(treatment != "Control") %>%
  mutate(compound_name = str_replace(compound_name, " 0\\(15N\\),? 0\\(13C\\)", ""),
         compound_name = str_replace(compound_name, "Homarine 0\\(15N,? 0\\(13C\\)", "Homarine"))



# Define custom colors
custom_colors <- c("#581845", "#900C3F", "#ee9090", "#90ee90", "#9090ee")

metabolite_colors <- c(
  "Glycine betaine" = custom_colors[5],
  "Homarine" = custom_colors[2],
  "Guanine" = custom_colors[1],
  "L- Glutamine" = custom_colors[4],
  "L-Glutamic acid" = custom_colors[3])



tbd_all_concs  %>%
  # filter(compound_name=="Homarine") %>%
  # filter(trpl!="C") %>%
  ggplot(aes(x=parti_nM, y=diss_nM, color=compound_name,
             group=interaction(trpl, compound_name))) +
  geom_abline(slope = 1, linetype=2) +
  geom_abline(slope = 1, linetype=2, intercept = 2) +
  geom_abline(slope = 1, linetype=2, intercept = -2) +
  geom_path() +
  geom_point() +
  geom_point(color="black", aes(alpha=hour)) +
  facet_wrap(~treatment, ncol=1) +
  scale_y_log10(limits=c(0.01, 2e3), name="Dissolved concentration") +
  scale_x_log10(limits=c(0.01, 2e3), name="Particulate concentration") +
  scale_alpha_continuous(range = c(0, 0.5)) +
  coord_fixed()




# Your ggplot code with custom colors
nfex_partivsdiss <- tbd_all_concs %>%
  group_by(trpl, treatment, compound_name) %>%
  mutate(parti_nM_end = lead(parti_nM), diss_nM_end=lead(diss_nM)) %>%
  ggplot(aes(x=parti_nM, y=diss_nM, xend=parti_nM_end, yend=diss_nM_end, color=compound_name)) +
  geom_abline(slope = 1, linetype=2, color="grey") +
  geom_abline(slope = 1, linetype=2, intercept = 2, color="grey") +
  geom_abline(slope = 1, linetype=2, intercept = -2, color="grey") +
  geom_segment(arrow = arrow(type="closed", length = unit(0.02, "npc"))) +
  facet_wrap(~treatment, ncol = 2, scales="free_x") +
  scale_y_log10(limits=c(0.01, 2e3), name="Dissolved concentration",
                breaks=c(0.01, 1, 100), labels=c(".01", 1, 100)) +
  scale_x_log10(limits=c(0.01, 50), name="Particulate concentration",
                breaks=c(0.01, 0.1, 1, 10, 100), labels=c(".01", ".1", 1, 10, 100)) +
  scale_alpha_continuous(range = c(0, 0.5)) +
  scale_color_manual(values = metabolite_colors) +  # Use custom colors
  theme_bw() +
  theme(axis.text = element_text(size = 24),
        strip.text = element_text(size = 24),
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.3),
        axis.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24)) +
  labs(color="Metabolite")


print(nfex_partivsdiss)



#----

ggsave(filename = "NFEXPartivsDiss_plot.png",
       plot = nfex_partivsdiss,
       width = 10.25, height = 11.20, units = "in")


