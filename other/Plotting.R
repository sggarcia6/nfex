#loading libraries

library(tidyverse)
library(ggplot2)
library(janitor)

# Putting the data together ----
diss_concs <- read_csv("intermediates/full_diss_final_concs.csv") %>%
  select(compound_name, filename, diss_nM=final_nM) %>%
  mutate(organism=str_extract(filename, "1005|1335")) %>%
  mutate(day=str_extract(filename, "D\\d+")) %>%
  mutate(tripl=str_extract(filename, "[A-C]$")) %>%
  filter(str_detect(filename, "Smp")) %>%
  filter(!str_detect(filename, "Blk"))
parti_concs <- read_csv("intermediates/full_parti_final_concs.csv") %>%
  select(compound_name, filename, part_nM=nM) %>%
  mutate(organism=str_extract(filename, "1005|1335")) %>%
  mutate(day=str_extract(filename, "D\\d+")) %>%
  mutate(tripl=str_extract(filename, "[A-C]$")) %>%
  filter(str_detect(filename, "Smp"))
rfus <- readxl::read_excel("data_raw/cultures_crr_rfu-fcm.xlsx", sheet = 2) %>%
  select(filename, sample_id, cpm=`cells/mL`, avg_RFU) %>%
  mutate(organism=str_extract(filename, "1005|1335")) %>%
  mutate(day=str_extract(filename, "D\\d+")) %>%
  mutate(tripl=str_extract(filename, "[A-C]$")) %>%
  filter(!sample_id%in%sprintf("D%02d", 4:12))
# Combined all the data into one place.
comb_df <- diss_concs %>%
  left_join(part_concs, by=c("compound_name", "organism", "day", "tripl"), suffix=c("_diss", "_part")) %>%
  left_join(rfus, by=c("organism", "day", "tripl")) %>%
  mutate(day=factor(day, levels=sprintf("D%01d", 0:11))) #%>% this is check
  filter(compound_name =="Homarine")

#Plotting  Dissolved
##Define the compound-to-group mapping
  compound_groups <- c(
    "(3-Carboxypropyl)trimethylammonium" = "Amino Acids",
    "4-Aminobutyric acid" = "Amino Acids",
    "5-Hydroxyectoine" = "Other",
    "5-Methylcytosine" = "Other",
    "Adenine" = "Nucleic Acids",
    "Adenosine" = "Nucleic Acids",
    "Betonicine" = "Betaines",
    "Butyrylcarnitine" = "Carnitines",
    "Citrulline" = "Amino Acids",
    "Creatine" = "Amino Acids",
    "Cytidine" = "Nucleic Acids",
    "Cytosine" = "Nucleic Acids",
    "Dimethylglycine" = "Betaines",
    "Ectoine" = "Amino Acids",
    "Glucosamine" = "Amino Sugars",
    "Glutamylphenylalanine" = "Amino Acids",
    "Glycine betaine" = "Betaines",
    "Guanine" = "Nucleic Acids",
    "Guanosine" = "Nucleic Acids",
    "Homarine" = "Betaines",
    "Hordenine" = "Amino Acids",
    "Hydroxyisoleucine" = "Amino Acids",
    "Hydroxyproline" = "Amino Acids",
    "L-Homoserine" = "Amino Acids",
    "Melamine" = "Other",
    "Muramic acid" = "Amino Sugars",
    "N6-Methyladenine" = "Nucleic Acids",
    "O-Methylmalonyl-L-carnitine" = "Carnitines",
    "Proline betaine" = "Betaines",
    "Sarcosine" = "Amino Acids",
    "Trigonelline" = "Amino Sugars",
    "beta-Alanine" = "Amino Acids",
    "beta-Alaninebetaine" = "Betaines",
    "beta-Glutamic acid" = "Amino Acids")

plot_diss_concs <- read_csv("intermediates/full_diss_final_concs.csv") %>%
    select(compound_name, filename, diss_nM = final_nM) %>%
    mutate(organism = str_extract(filename, "1005|1335")) %>%
    mutate(day = str_extract(filename, "D\\d+")) %>%
    mutate(tripl = str_extract(filename, "[A-C]$")) %>%
    filter(str_detect(filename, "Smp")) %>%
    filter(!str_detect(filename, "Blk")) %>%
    group_by(compound_name, organism, day) %>%
    mutate(avg_nM = mean(diss_nM, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(category = case_when(
      compound_name %in% names(compound_groups) ~ compound_groups[compound_name],
      TRUE ~ "Other"
    ))

#plotting the categories of compounds in each organism
category_plot <- ggplot(plot_diss_concs, aes(x = category, y = avg_nM, fill = category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") +
  labs(x = "Category", y = "Average avg_nM") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  guides(fill = guide_legend(title = "Category")) +
  ggtitle("Average avg_nM by Category") +
  facet_wrap(~ organism, ncol = 2, scales = "free_y")


day_order <- c("D0", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11")

# Created a plot with custom day order and dots
avg_nMbycategory <- ggplot(plot_diss_concs, aes(x = factor(day, levels = day_order), y = avg_nM, color = compound_name, group = 1)) +
  geom_point(size = 3) +
  ggh4x::facet_grid2(organism ~ category, scales = "free", independent = "all") +
  scale_x_discrete(labels = day_order) +  # Set the custom day order
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom")
avg_nMbycategory
#focusing on just the betaines

betaines_diss_concs <- plot_diss_concs %>%
  filter(category == "Betaines") %>%
  filter(compound_name != "Betonicine")

my_colors <- c( "pink", "orange")

# Create a plot for "Betaines" trends
betaines_plot <- ggplot(betaines_diss_concs, aes(x = factor(day, levels = day_order), y = avg_nM, color = organism)) +
  geom_line(size = 1.5) +  # Increase line size for better visibility
  geom_point(size = 4) +   # Increase point size for better visibility
  facet_wrap(organism ~ compound_name, scales = "free_y") +
  scale_x_discrete(labels = day_order) +
  scale_color_manual(values = my_colors) +  # Use custom colors
  labs(x = "Day", y = "Average avg_nM", color = "Organism") +
  ggtitle("Trends of Betaines by Organism and Compound") +
  theme_light()  # Use a lighter theme for the background

#COMBINED (normalized by frank) ----
diss_comb_df <- diss_concs %>%
  left_join(rfus, by = c("organism", "day", "tripl")) %>%
  mutate(day = factor(day, levels = sprintf("D%01d", 0:11))) %>%
  mutate(normalized_diss_nM = diss_nM / cpm) %>%
  group_by(compound_name, organism, day) %>%
  mutate(normalized_avg_nM = mean(normalized_diss_nM, na.rm = TRUE)) %/%
  mutate(category = case_when(
    compound_name %in% names(compound_groups) ~ compound_groups[compound_name],
    TRUE ~ "Other"
  ))

# Created a plot with custom day order and dots
comb_diss <- ggplot(diss_comb_df, aes(x = factor(day, levels = day_order), y = normalized_avg_nM, color = compound_name, group = 1)) +
  geom_point(size = 3) +
  ggh4x::facet_grid2(organism ~ category, scales = "free", independent = "all") +
  scale_x_discrete(labels = day_order) +  # Set the custom day order
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom")














#----
comb_df %>%
  ggplot(aes(label=filename_diss)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x=part_nM, y=diss_nM, color=day, shape=tripl), size=4) +
  ggh4x::facet_grid2(organism~compound_name, scales="free", independent = "all") +
  # scale_y_log10() + scale_x_log10() +
  scale_shape_identity() +
  theme()

rfu_gp <- rfus %>%
  mutate(day=factor(day, levels=sprintf("D%01d", 0:11))) %>%
  arrange(day) %>%
  mutate(filename=fct_inorder(filename)) %>%
  distinct(avg_RFU, filename, organism, cpm, day) %>%
  pivot_longer(c(avg_RFU, cpm)) %>%
  ggplot() +
  geom_col(aes(x=filename, y=value, fill=day)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_grid(name~organism, scales="free") +
  theme(legend.position = "top")


rfus %>%
  mutate(day = factor(day, levels = sprintf("D%01d", 0:11))) %>%
  arrange(day) %>%
  mutate(filename = fct_inorder(filename)) %>%
  distinct(RFU = avg_RFU, filename, organism, Cells_per_mL = cpm, day) %>%
  pivot_longer(c(RFU, Cells_per_mL)) %>%
  ggplot() +
  geom_col(aes(x = filename, y = value, fill = day)) +
  labs(
    x = "Filename",
    y = "",
    fill = "Day",
    title = "RFU and Cells/mL by Filename and Day"  # Add a title
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  ) +
  facet_grid(name ~ organism, scales = "free")


metab_gp <- comb_df %>%
  ggplot(aes(label=filename_diss)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x=part_nM, y=diss_nM, color=day, shape=tripl), size=4) +
  ggh4x::facet_grid2(organism~compound_name, scales="free", independent = "all") +
  # scale_y_log10() + scale_x_log10() +
  scale_shape_identity() +
  theme(legend.position="top")

gp_together <- egg::ggarrange(rfu_gp, metab_gp, widths = c(0.3, 0.7))
ggsave("~/../Desktop/gp_together.pdf", plot = gp_together, device = "pdf",
       width = 18, height=8, units = "in")
