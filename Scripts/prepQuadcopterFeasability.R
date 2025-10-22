library(here)
library(dsims)
library(dssd)
library(scales)

# Define the path to the folder
folder_path <- here("Output", "Simulation")

# List all files starting with 'QC_Sys_design_'
files <- list.files(path = folder_path, pattern = "density_sim-A.*\\.RData$", full.names = TRUE)

# Extract the required numbers from each file name and load the file to extract the summary value
stats <- lapply(files, function(file) {
  # Regular expression to extract the desired components
  matches <- regmatches(file, regexec("QC_([^_]+)_design_([^\\-]+)-density_sim-A(\\d+)", file))

  # Convert to data frame for clarity
  simConfigs <- do.call(rbind, lapply(matches, function(x) data.frame(Overlap = x[2], NSubplots = x[3], Abundance = x[4])))
  # if(is.na(simConfigs[1,1])){simConfigs[1,1] <- "gimbal"}

  # Load the file
  load(file)

  # Extract the summary value
  if (is.na(QC_Sys_sim_density)) {
    mean_n <- NA
    true_D <- NA
    area <- NA
    subplots <- NA
    sd_of_means <- NA
    effortKM <- NA
  } else {
    summary_data <- summary(QC_Sys_sim_density, description.summary = FALSE)
    mean_n <- round(summary_data@individuals$summary$mean.n, 0)
    true_D <- round(summary_data@individuals$D$Truth * 1e9, 3)
    area <- round(summary_data@individuals$summary$mean.Cover.Area / 1e6, 0)
    subplots <- length(summary_data@design.summary$design.type)
    sd_of_means <- round(summary_data@individuals$D$sd.of.means / summary_data@individuals$D$mean.Estimate, 3)
    effortKM <- round((summary_data@individuals$summary$mean.Effort) / 1000, 2)
  }
  data.frame(subplots = subplots, mean_n = mean_n, area = area, trueDensity = true_D, CV = sd_of_means, NSubplots = simConfigs[1, 2], Abundance = simConfigs[1, 3], EffortKM = effortKM, Overlap = simConfigs[1, 1])
})

# Combine the list of data frames into a single data frame
sum_stats_df <- do.call(rbind, stats)

# Convert all columns in extracted_df to numeric
sum_stats_df[, 1:8] <- lapply(sum_stats_df[, 1:8], as.numeric)


# Scale the area column so that 0 remains 0 and other values are scaled between 0 and 1
max_area <- max(sum_stats_df$area, na.rm = TRUE)
sum_stats_df$area_p <- round(sum_stats_df$area / max_area, 2)

# Print the resulting data frame to verify the changes
print(sum_stats_df)


# Write the dataframe to a CSV file
output_path <- here("Output", "Simulation", "quadcopterFeasability.csv")
write.csv(extracted_df, file = output_path, row.names = FALSE)
