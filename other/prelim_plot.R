#PLOTTING

#preliminary plotting
final_samplepeaks <- BMISed_areas %>%
  filter(str_starts(filename, "230810_Smp_G5_NFEX_T")) %>%
  mutate(day=str_extract(filename, "T\\d+")) %>%
  mutate(time=as.numeric(str_remove(day, "T"))) %>%
  select(filename, compound_name, bmised_area, time)






# Assuming 'final_samplepeaks' has columns 'filename', 'compound_name', 'bmised_area', and 'time'

# Plotting
ggplot(final_samplepeaks, aes(x = compound_name, y = bmised_area, fill = factor(time))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Changes in bmised_area Over Time",
       x = "Compound Name",
       y = "BMISed Area",
       fill = "Time (hrs)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Filter data for compounds containing "Glycine betaine" in their names
glycine_betaine_data <- final_samplepeaks %>%
  filter(str_detect(compound_name, "Glycine betaine"))

# Plotting for compounds containing "Glycine betaine"
ggplot(glycine_betaine_data, aes(x = compound_name, y = bmised_area, fill = factor(time))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Changes in Compounds Containing Glycine Betaine Over Time",
       x = "Time (hrs)",
       y = "BMISed Area",
       fill = "Time (hrs)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Filter data for compounds containing "Glycine betaine" in their names
homarine_data <- final_samplepeaks %>%
  filter(str_detect(compound_name, "Homarine"))

# Plotting for compounds containing "Glycine betaine"
ggplot(homarine_data, aes(x = compound_name, y = bmised_area, fill = factor(time))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Changes in Compounds Containing Homarine Over Time",
       x = "Time (hrs)",
       y = "BMISed Area",
       fill = "Time (hrs)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Filter data for compounds containing "Glycine betaine" in their names
guanine_data <- final_samplepeaks %>%
  filter(str_detect(compound_name, "Guanine"))

# Plotting for compounds containing "Glycine betaine"
ggplot(guanine_data, aes(x = compound_name, y = bmised_area, fill = factor(time))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Changes in Compounds Containing Guanine Over Time",
       x = "Time (hrs)",
       y = "BMISed Area",
       fill = "Time (hrs)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Filter data for compounds containing "Glycine betaine" in their names
glutamine_data <- final_samplepeaks %>%
  filter(str_detect(compound_name, "L- Glutamine"))

# Plotting for compounds containing "Glycine betaine"
ggplot(glutamine_data, aes(x = compound_name, y = bmised_area, fill = factor(time))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Changes in Compounds Containing Glutamine Over Time",
       x = "Time (hrs)",
       y = "BMISed Area",
       fill = "Time (hrs)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

glutamic_data <- final_samplepeaks %>%
  filter(str_detect(compound_name, "L-Glutamic acid"))

# Plotting for compounds containing "Glycine betaine"
ggplot(glutamic_data, aes(x = compound_name, y = bmised_area, fill = factor(time))) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Changes in Compounds Containing Glutamic acid Over Time",
       x = "Time (hrs)",
       y = "BMISed Area",
       fill = "Time (hrs)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





#Big showcase
specific_hours <- c(0, 6, 12, 24)
ggplot(final_samplepeaks, aes(x = , y = bmised_area, color = compound_name)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 1) +
  labs(title = "Particulate BMISed Areas by Compound and Time",
       x = "Time (hours)",
       y = "BMISed Area") +
  theme_minimal() +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_x_continuous(breaks = specific_hours) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = FALSE)  # Remove the legend for the color aesthetic








