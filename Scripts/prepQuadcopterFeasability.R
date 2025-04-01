library(here)
library(dsims)
library(dssd)

# Define the path to the folder
folder_path <- here("Output", "Simulation")

# List all files starting with 'QC_Sys_design_'
files <- list.files(path = folder_path, pattern = "^QC_Sys_design_", full.names = TRUE)

# Extract the required numbers from each file name and load the file to extract the summary value
extracted_numbers <- lapply(files, function(file) {
  design_number <- sub(".*design_(\\d{2}).*", "\\1", file)
  d_number <- sub(".*-D(\\d+\\.\\d+|\\d+).*", "\\1", file)

  # Load the file
  load(file)

  # Extract the summary value
  mean_n <- summary(QC_Sys_sim_density, description.summary = FALSE)@individuals$
    summary$
    mean.n
  print(QC_Sys_sim_density@mean_cover_area)

  data.frame(subplots = design_number, density = d_number, mean_n = mean_n)
})

# Combine the list of data frames into a single data frame
extracted_df <- do.call(rbind, extracted_numbers)

# Convert all columns in extracted_df to numeric
extracted_df[] <- lapply(extracted_df, as.numeric)

# Remove rows with density > 16
filtered_df <- subset(extracted_df, density < 9)

# Print the resulting data frame
print(filtered_df)

# Write the dataframe to a CSV file
output_path <- here("Output", "Simulation", "quadcopterFeasability.csv")
write.csv(filtered_df, file = output_path, row.names = FALSE)

design_comparison_df