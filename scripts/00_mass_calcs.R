# Define the compounds and their corresponding masses
compounds <- c("Choline", "Carnitine", "Proline", "Proline Betaine", "Lysine",
               "Trimethylamine N-Oxide", "Glycine", "Dimethylglycine",  "Ectoine")
masses <- c(104.107539070-1.007276, 161.10519334, 115.063328530, 143.094628657, 146.105527694,
            75.068413911, 75.032028402, 103.063328530, 142.074227566) + 1.007276

# Define the constants
N15 <- 0.997035
C13 <- 1.003355
O18 <- 2.004244

# Function to calculate isotopes and return a dataframe
calculate_isotopes_df <- function(compound, mass, N15, C13) {
  compound_0_15N_0_13C <- mass
  compound_1_15N_0_13C <- compound_0_15N_0_13C + N15
  compound_2_15N_0_13C <- if (compound %in% c("Ectoine", "Lysine")) compound_1_15N_0_13C + N15 else NA
  compound_0_15N_1_13C <- compound_0_15N_0_13C + C13
  compound_1_15N_1_13C <- compound_0_15N_0_13C + C13 + N15
  compound_2_15N_1_13C <- if (compound %in% c("Ectoine", "Lysine")) compound_1_15N_1_13C + N15 else NA

  results <- data.frame(
    Compound_name = c(paste0(compound, " 0(15N), 0(13C)"), paste0(compound, " 1(15N), 0(13C)"), paste0(compound, " 2(15N), 0(13C)"), paste0(compound, " 0(15N), 1(13C)"), paste0(compound, " 1(15N), 1(13C)"), paste0(compound, " 2(15N), 1(13C)")),
    Value = c(compound_0_15N_0_13C, compound_1_15N_0_13C, compound_2_15N_0_13C, compound_0_15N_1_13C, compound_1_15N_1_13C, compound_2_15N_1_13C),
    Isotope_Change = c("No change", "15N added", ifelse(compound %in% c("Ectoine", "Lysine"), "15N added", "No change"), "13C added", "15N and 13C added", ifelse(compound %in% c("Ectoine", "Lysine"), "15N and 13C added", "No change"))
  )

  return(results)
}


# Call the function for each compound and combine the results into one dataframe
all_results <- do.call(rbind, Map(calculate_isotopes_df, compounds, masses, list(N15), list(C13))) %>%
  na.omit(all_results)

# Print the combined dataframe
print(all_results)


write_csv(all_results,"intermediates/RAMS_mass_calcs.csv")
