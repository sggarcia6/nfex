library(tidyverse)
library(ggplot2)
library(janitor)

#Creating the filepath for file
filepath <- "intermediates/full_diss_final_concs.csv"


#Having an original file of the file received directly from Skyline
cncs_raw <- read.csv(filepath, check.names = TRUE) %>%
  clean_names()

# wiggling the data
cncs <- cncs_raw %>%
  filter(grepl("^230710_Smp_1", filename)) %>%
  mutate(species=str_extract(filename, "1335|1005")) %>%
  # filter(species=="1335") %>%
  mutate(group = sub("_[A-Z]$", "", filename)) %>%
  group_by(group, compound_name, species) %>%
  summarise(avg_nM = mean(final_n_m, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(day = as.integer(sub(".*D(\\d+)_.*", "\\1", group)))

#Plotting the data (log scale)
ggplot(cncs, aes(x = day, y = avg_nM)) +
  geom_line(size = 1.2, aes(color = compound_name), show.legend = FALSE) +
  geom_point(size = 3, aes(color = compound_name)) +
  labs(x = "Day", y = "Average nM") +
  ggtitle("Average nM for Unique Combinations of Compounds") +
  # facet_wrap(~ compound_name, scales = "free_y", ncol = 2) +
  facet_grid(~species, scales = "free_y") +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom")

#plotting data (normal scale)
ggplot(cncs, aes(x = day, y = avg_nM)) +
  geom_line(size = 1.2, aes(color = compound_name), show.legend = FALSE) +
  geom_point(size = 3, aes(color = compound_name)) +
  facet_wrap(species~compound_name, scales = "free_y") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom")



# Filter the data to include only 'L-Glutamic acid'
compound_to_plot <- "L-Glutamic acid"
filtered_data <- betaine_cncs %>%
  filter(compound_name == compound_to_plot)

# Calculate the slopes and intercepts for 'L-Glutamic acid'
slopes_intercepts <- filtered_data %>%
  summarise(slope = coef(lm(avg_nM ~ day))[2], intercept = coef(lm(avg_nM ~ day))[1])

# Create a plot for 'L-Glutamic acid' with a slope line
ggplot(filtered_data, aes(x = day, y = avg_nM)) +
  geom_line(size = 1.2, aes(color = compound_name), show.legend = FALSE) +
  geom_point(size = 3, aes(color = compound_name)) +
  labs(x = "Day", y = "Average nM") +
  ggtitle("Average nM for L-Glutamic acid") +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3)) +  # Format y-axis labels
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom") +
  geom_abline(data = slopes_intercepts, aes(intercept = intercept, slope = slope), color = "black", linetype = 2)



#plotting data but in barplot mode to see variation
betaine_cncs_raw %>%
  filter(grepl("^230710_Smp_1", filename)) %>%
  mutate(species=str_extract(filename, "1335|1005")) %>%
  mutate(day = as.integer(sub(".*D(\\d+)_.*", "\\1", filename))) %>%
  mutate(tripl=str_extract(filename, "[A-C]$")) %>%
  ggplot(aes(x = day, y = n_m)) +
  geom_hline(yintercept = 0) +
  # geom_boxplot(aes(group=day)) +
  # geom_point(aes(color=tripl)) +
  geom_col(aes(fill=tripl), position = "dodge") +
  facet_grid(compound_name~species, scales = "free_y") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom") %>%
  plotly::ggplotly()


# plotting

# Create a bar plot that divide by compound
ggplot(betaine_cncs, aes(x = day, y = avg_nM, fill = compound_name)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Day", y = "Average nM") +
  ggtitle("Concentrations of Compounds by Day") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "top") +
  facet_wrap(~compound_name, scales = "free_y")

#another plot but looks like katherine's


# calculating the total concentration for each compound within each day and diatom
total_concentration <- cncs %>%
  group_by(day, species, compound_name) %>%
  summarise(total_avg_nM = sum(avg_nM))

# custom color palette for compounds from internet

custom_palette <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#fdbf6f", "#b2df8a")

# a prettier bar plot for total concentration
ggplot(total_concentration, aes(x = as.factor(day), y = total_avg_nM, fill = compound_name)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Day", y = "Total Average nM") +
  ggtitle("Total Concentration of Compounds by Day and Diatom") +
  scale_fill_manual(values = custom_palette) +  # Custom color palette
  facet_grid(~species, scales = "free_y", labeller = labeller(species = c("1005" = "T. Oceanica", "1335" = "T. Pseudonana"))) +
  theme_minimal() +  # Minimal theme
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.spacing = unit(0.5, "lines")) +
  guides(fill = guide_legend(ncol = 2)) %>%
  plotly::ggplotly()



#IDK
# Replace total_concentration with your new dataset name
new_total_concentration <- cncs %>%
  group_by(day, species, compound_name) %>%
  summarise(total_avg_nM = sum(avg_nM))


# Get unique compound names and generate colors for each
unique_compounds <- unique(new_total_concentration$compound_name)
num_compounds <- length(unique_compounds)
color_palette <- viridis(num_compounds)

# Create a named vector to map compound names to colors
compound_colors <- setNames(color_palette, unique_compounds)

# Use the named vector for fill colors
ggplot(new_total_concentration, aes(x = as.factor(day), y = total_avg_nM, fill = compound_name)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Day", y = "Total Average nM") +
  ggtitle("Total Concentration of Compounds by Day and Diatom") +
  scale_fill_manual(values = compound_colors) +  # Use the generated colors
  facet_grid(~species, scales = "free_y", labeller = labeller(species = c("1005" = "T. Oceanica", "1335" = "T. Pseudonana"))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.spacing = unit(0.5, "lines")) +
  guides(fill = guide_legend(ncol = 2))














plotly::ggplotly()
