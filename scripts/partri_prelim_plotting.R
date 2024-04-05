#Loading Libraries

library(ggplot2)
library(dplyr)

#Prelim plotting

parti_final_concs <- read_csv("intermediates/NFEX_parti_final_concs.csv") %>%
  filter(!str_detect(filename, "Std"),
         !str_detect(filename, "Poo"),
         !str_detect(filename, "Blk")) %>%
  mutate(day = str_extract(filename, "T\\d+"),
         day = as.numeric(str_remove(day, "T")),
         treatment = str_extract(filename, "NO3|NH4|Urea")) %>%
  mutate(treatment = ifelse(is.na(treatment), "Control", treatment)) %>%
  group_by(compound_name, day, treatment) %>%
  summarise(mean_nM = mean(nM)) %>%
  ungroup() %>%
  mutate(N = as.integer(str_extract(compound_name, "\\d+(?=\\(15N\\))")),
         C = as.integer(str_extract(compound_name, "\\d+(?=\\(13C\\))"))) %>%
  mutate(mean_nM = replace_na(mean_nM, 0),
         N = replace_na(N, 0),
         C = replace_na(C, 0)) %>%
  #select(-filename) %>%
  unique()

# Convert day into factor so it can be better labeled
parti_final_concs$day <- factor(parti_final_concs$day)

# Identifier for the compounds
parti_final_concs$compound_id <- as.integer(as.factor(parti_final_concs$compound_name))

# #Mass Isotope Distribution (MID) for NO3 ----
#
# GBT_NO3_final_concs_plot <- final_concs %>%
#   filter(treatment == "NO3", grepl("Glycine betaine", compound_name)) %>%
#   ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
#   geom_col(position = "fill", color = "black")
#
# Guanine_NO3_final_concs_plot <- final_concs %>%
#   filter(treatment == "NO3", grepl("Guanine", compound_name)) %>%
#   ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
#   geom_col(position = "fill", color = "black")
#
# Homarine_NO3_final_concs_plot <- final_concs %>%
#   filter(treatment == "NO3", grepl("Homarine", compound_name)) %>%
#   ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
#   geom_col(position = "fill", color = "black")
#
# Glutamic_NO3_final_concs_plot <- final_concs %>%
#   filter(treatment == "NO3", grepl("L-Glutamic acid", compound_name)) %>%
#   ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
#   geom_col(position = "fill", color = "black")
#
# Glutamine_NO3_final_concs_plot <- parti_final_concs %>%
#   filter(treatment == "NO3", grepl("L- Glutamine", compound_name)) %>%
#   ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
#   geom_col(position = "fill", color = "black")
#
# Glutamine_NO3_final_concs_plot


# size and settings
options(repr.plot.width = 10, repr.plot.height = 6)
theme_set(theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "bottom"))

# make groups seems to be repetitive
part_final_concs <- parti_final_concs %>%
  mutate(compound_group = case_when(
    grepl("Glycine betaine", compound_name) ~ "Glycine betaine",
    grepl("Guanine", compound_name) ~ "Guanine",
    grepl("Homarine", compound_name) ~ "Homarine",
    grepl("L-Glutamic acid", compound_name) ~ "Glutamic Acid",
    grepl("L- Glutamine", compound_name) ~ "Glutamine",
    TRUE ~ "Other"
  ))


# Filter and plot for all compounds grouped by compound_group
NO3_compound_group_plot <- part_final_concs %>%
  mutate(iso_label=paste0(N, "x15N, ", C, "x13C")) %>%
  filter(treatment == "NO3") %>%
  ggplot(aes(x = day, y = mean_nM, fill = iso_label)) +
  geom_col(position = "fill", color = "black", alpha = 0.9) +
  facet_wrap(~ compound_group, scales = "free_y", nrow = 1) +
  labs(title = "Particulate Metabolite Mass Isotope Distribution treated by NO3",
       x = "Hour",
       y = "Mean Concentration (nM)",
       fill = "Isotope Label") +
  theme(legend.key.size = unit(0.5, "cm")) +
  scale_fill_manual(values = c("#581845", "#900C3F", "#C70039", "#FF5733", "#FFC30F","#9290C3" ))

NO3_compound_group_plot



#Mass Isotope Distribution (MID) for Urea ----

GBT_Urea_final_concs <- final_concs %>%
  filter(treatment == "Urea", grepl("Glycine betaine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Guanine_Urea_final_concs_plot <- final_concs %>%
  filter(treatment == "Urea", grepl("Guanine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Homarine_Urea_final_concs_plot <- final_concs %>%
  filter(treatment == "Urea", grepl("Homarine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Glutamic_Urea_final_concs_plot <- final_concs %>%
  filter(treatment == "Urea", grepl("L-Glutamic acid", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Glutamine_Urea_final_concs_plot <- final_concs %>%
  filter(treatment == "Urea", grepl("L- Glutamine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")



# size and settings
options(repr.plot.width = 10, repr.plot.height = 6)
theme_set(theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "bottom"))

# make groups seems to be repetitive
final_concs <- final_concs %>%
  mutate(compound_group = case_when(
    grepl("Glycine betaine", compound_name) ~ "Glycine betaine",
    grepl("Guanine", compound_name) ~ "Guanine",
    grepl("Homarine", compound_name) ~ "Homarine",
    grepl("L-Glutamic acid", compound_name) ~ "Glutamic Acid",
    grepl("L- Glutamine", compound_name) ~ "Glutamine",
    TRUE ~ "Other"
  ))


# Filter and plot for all compounds grouped by compound_group
Urea_compound_group_plot <- final_concs %>%
  mutate(iso_label=paste0(N, "x15N, ", C, "x13C")) %>%
  filter(treatment == "Urea") %>%
  ggplot(aes(x = day, y = mean_nM, fill = iso_label)) +
  geom_col(position = "fill", color = "black") +
  facet_wrap(~ compound_group, scales = "free_y", nrow = 1) +
  labs(title = "Particulate Metabolite Mass Isotope Distribution treated by Urea",
       x = "Hour",
       y = "Mean Concentration (nM)",
       fill = "Isotope Label") +
  theme(legend.key.size = unit(0.5, "cm")) +
  scale_fill_manual(values = c("#581845", "#900C3F", "#C70039", "#FF5733", "#FFC30F","#9290C3" ))

Urea_compound_group_plot

#Mass Isotope Distribution (MID) for NH4 ----

GBT_NH4_final_concs <- final_concs %>%
  filter(treatment == "NH4", grepl("Glycine betaine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Guanine_NH4_final_concs_plot <- final_concs %>%
  filter(treatment == "NH4", grepl("Guanine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Homarine_NH4_final_concs_plot <- final_concs %>%
  filter(treatment == "NH4", grepl("Homarine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Glutamic_NH4_final_concs_plot <- final_concs %>%
  filter(treatment == "NH4", grepl("L-Glutamic acid", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Glutamine_NH4_final_concs_plot <- final_concs %>%
  filter(treatment == "NH4", grepl("L- Glutamine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")


# size and settings
options(repr.plot.width = 10, repr.plot.height = 6)
theme_set(theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "bottom"))

# make groups seems to be repetitive
final_concs <- final_concs %>%
  mutate(compound_group = case_when(
    grepl("Glycine betaine", compound_name) ~ "Glycine betaine",
    grepl("Guanine", compound_name) ~ "Guanine",
    grepl("Homarine", compound_name) ~ "Homarine",
    grepl("L-Glutamic acid", compound_name) ~ "Glutamic Acid",
    grepl("L- Glutamine", compound_name) ~ "Glutamine",
    TRUE ~ "Other"
  ))


# Filter and plot for all compounds grouped by compound_group
NH4_compound_group_plot <- final_concs %>%
  mutate(iso_label=paste0(N, "x15N, ", C, "x13C")) %>%
  filter(treatment == "NH4") %>%
  ggplot(aes(x = day, y = mean_nM, fill = iso_label)) +
  geom_col(position = "fill", color = "black") +
  facet_wrap(~ compound_group, scales = "free_y", nrow = 1) +
  labs(title = "Particulate Metabolite Mass Isotope Distribution treated by NH4",
       x = "Hour",
       y = "Mean Concentration (nM)",
       fill = "Isotope Label") +
  theme(legend.key.size = unit(0.5, "cm")) +
  scale_fill_manual(values = c("#581845", "#900C3F", "#C70039", "#FF5733", "#FFC30F","#9290C3" ))

NH4_compound_group_plot


title = "Particulate Metabolite Mean Concentration by Treatment and Hour"



#Mass Isotope Distribution (MID) for Control ----


GBT_Control_final_concs <- final_concs %>%
  filter(treatment == "Control", grepl("Glycine betaine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Guanine_Control_final_concs_plot <- final_concs %>%
  filter(treatment == "Control", grepl("Guanine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Homarine_Control_final_concs_plot <- final_concs %>%
  filter(treatment == "Control", grepl("Homarine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Glutamic_Control_final_concs_plot <- final_concs %>%
  filter(treatment == "Control", grepl("L-Glutamic acid", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")

Glutamine_Control_final_concs_plot <- final_concs %>%
  filter(treatment == "Control", grepl("L- Glutamine", compound_name)) %>%
  ggplot(aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_col(position = "fill", color = "black")




# Convert day into factor so it can be better labeled
final_concs$day <- factor(final_concs$day)

# Identifier for the compounds
final_concs$compound_id <- as.integer(as.factor(final_concs$compound_name))




#Bargraph ----


# Bar graph for each compound's absolute concentrations using facets
ggplot(final_concs, aes(x = day, y = mean_nM, fill = compound_name)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Particulate Metabolite Mean Concentration by Treatment and Hour",
       x = "Hour",
       y = "Mean Concentration (nM)",
       fill = "Treatment") +
  theme_minimal() +
  scale_fill_manual(values = c("NO3" = "#E05793", "NH4" = "#FE927D", "Urea" = "#FECF5D", "Control" = "#D2B9DE"),
                    labels = c("NO3" = "NO3", "NH4" = "NH4", "Urea" = "Urea", "Control" = "Control")) +
  facet_grid(treatment~compound_name, scales = "free_y", nrow = 5) +  # Facet by compound_name with 2 rows
  guides(fill = guide_legend(title = "Treatment", ncol = 2))  # Adjust legend appearance





