library(here)
library(dsims)
library(dssd)
library(scales)

# Define the path to the folder
folder_path <- here("Output", "Simulation")

# List all files starting with 'QC_Sys_design_'
files <- list.files(path = folder_path, pattern = "^QC_Sys_design_", full.names = TRUE)

# Extract the required numbers from each file name and load the file to extract the summary value
extracted_numbers <- lapply(files, function(file) {
  d_number <- sub(".*-D(\\d+\\.\\d+|\\d+).*", "\\1", file)

  # Load the file
  load(file)

  # Extract the summary value
  summary_data <- summary(QC_Sys_sim_density, description.summary = FALSE)
  mean_n <- round(summary_data@individuals$summary$mean.n, 0)
  true_D <- round(summary_data@individuals$D$Truth * 1e9, 3)
  area <- round(summary_data@individuals$summary$mean.Cover.Area / 1e6, 0)
  subplots <- length(summary_data@design.summary$design.type)
  sd_of_means = round(summary_data@individuals$D$sd.of.means / summary_data@individuals$D$mean.Estimate, 3)

  data.frame(subplots = subplots, mean_n = mean_n, area = area, trueDensity = true_D, d_number = d_number, CV = sd_of_means)
})

# Combine the list of data frames into a single data frame
extracted_df <- do.call(rbind, extracted_numbers)

# Convert all columns in extracted_df to numeric
extracted_df[] <- lapply(extracted_df, as.numeric)

# Remove rows with density > 8
filtered_df <- subset(extracted_df, d_number < 9)

# Scale the area column so that 0 remains 0 and other values are scaled between 0 and 1
max_area <- max(filtered_df$area)
filtered_df$area_p <- round(filtered_df$area / max_area, 2)

# Print the resulting data frame to verify the changes
print(filtered_df)


# Write the dataframe to a CSV file
output_path <- here("Output", "Simulation", "quadcopterFeasability.csv")
write.csv(filtered_df, file = output_path, row.names = FALSE)
