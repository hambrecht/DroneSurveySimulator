# load libraries
library(here)
library(knitr)
library(sf)
library(dplyr)
library(GGally)
library(dssd)
# library(reshape2)

# Create a named vector with old and new group names
group_names <- c(
  "FW_Sys_2C" = "FW-Sys-2C",
  "FW_ZZ_2C" = "FW-ZZ-2C",
  "FW_Sys_G" = "FW-Sys-G",
  "FW_ZZ_G" = "FW-ZZ-G",
  "QC_Sys_NADIR" = "QC-Sys-NADIR"
  "QC_Sys" = "QC-Sys"
  "H_SG" = "H-SG"
)

# List all files containing 'cover' in folder Output/Simulation
files <- list.files(path = here("Output", "Simulation"), pattern = "cover-WMU", full.names = TRUE)
wmu_path_list <- list.files(path = here("Output", "PrepData"), pattern = "prepared", full.names = TRUE)
# wmu_path_list <- wmu_path_list[-4]

# Extract the last three digits before the file extension for each file
file_ids <- sapply(files, function(file) {
    match <- regmatches(file, regexpr("\\d{3}(?=\\.\\w+$)", file, perl = TRUE))
    return(match)
})


# Load the first file and list its contents to verify object names
# Combine names with the first item in file_ids

for (i in seq_along(files)) {
    before_load <- ls()
    load(files[i])
    load(wmu_path_list[i])
  rm(segdata, distdata, obsdata)
    after_load <- ls()
    new_objs <- setdiff(after_load, before_load)
    suffix <- file_ids[[i]]
    for (obj in new_objs) {
        assign(paste0(obj, "_", suffix), get(obj))
        rm(list = obj)
    }
}
ls()


# Merge all objects into one dataframe and add `i` as a column to identify the origins
merged_df <- data.frame()

for (i in 1:length(file_ids)) {
    obj_name <- paste0("design_comparison_df_", file_ids[i])
    if (exists(obj_name)) {
        temp_df <- get(obj_name)
      # temp_df$Simulation <- paste0(temp_df$Simulation, "_", file_ids[i])
      temp_df$WMU <- file_ids[i]
        merged_df <- bind_rows(merged_df, temp_df)
    }
}


# Sort merged_df by the `Simulation` column
merged_df <- merged_df %>% arrange(Simulation)

# # Remove all rows with 'Short-Zig' in the 'Simulation' column
# merged_df <- merged_df %>%
#   filter(Simulation != "Short-Zig")

# Convert Mean_Cover_Area from m² to km² and round to 2 decimal places
merged_df$Mean_Cover_Area <- round(merged_df$Mean_Cover_Area / 1e6, 2)

# Convert specified columns from m to km and round to 2 decimal places
columns_to_convert <- c("Mean_Line_Length", "Mean_Trackline_Length", "Mean_Cyclic_Trackline_Length", 
                        "Mean_On_Effort", "Mean_Off_Effort", "Mean_Return_to_Home", "Mean_Off_Effort_Return")

merged_df[columns_to_convert] <- round(merged_df[columns_to_convert] / 1e3, 2)

kable(merged_df)



# Reorder columns and sort by WMU
merged_df <- merged_df %>%
  select(WMU, everything()) %>%
  arrange(WMU)

kable(merged_df)

# Identify list columns
list_columns <- sapply(merged_df, is.list)
# Convert list columns to character
merged_df[list_columns] <- lapply(merged_df[list_columns], as.character)

merged_df

# Save the transposed dataframe as a CSV file
output_path <- here("Output", "Simulation", "cover_overview.csv")
write.csv(merged_df, file = output_path, row.names = TRUE)


### Move to Python `vis_cover.ipynb`

# Assuming merged_df is your data frame and it has columns WMU, Simulation, and Mean_Line_Length
matrix_df <- dcast(merged_df, WMU ~ Simulation, value.var = "Mean_Line_Length")

# Print the matrix
kable(matrix_df)


#####



# Print the merged dataframe
## design stats
## For details see: https://examples.distancesampling.org/dssd-getting-started/GettingStarted-distill.html#appendix-trackline-and-cyclic-trackline-lengths
## Line lenght = on effort line length
## The trackline length is the sum of the lengths of the transects plus the off-effort transit time required to complete the survey from the beginning of the first transect to the end of the last transect. The off-effort transit distance is calculated as the crow flies and may be longer in reality if transit is required around lakes, islands or coastlines etc.
## The cyclic trackline length is the trackline length plus the off-effort transit distance required to return from the end of the last transect to the beginning of the first transect.


summary_by_simulation <- merged_df %>%
  group_by(Simulation) %>%
  summarise(across(where(is.numeric), \(x) round(mean(x, na.rm = TRUE), 1))) %>%
  select(1, 12:15)
# mutate(Simulation = recode(Simulation, !!!sim_names))


# Print the summarized table
kable(summary_by_simulation)





# Initialize an empty dataframe to store the results
coverage_stats <- data.frame()

# Loop through all objects in the environment
for (obj_name in ls()) {
  # Check if the object name contains 'design' but not 'comparison'
  if (grepl("design", obj_name) && !grepl("comparison", obj_name)) {
    # Get the coverage score
    coverage_score <- get.coverage(get(obj_name))


    # Calculate summary statistics manually
    min_val <- min(coverage_score, na.rm = TRUE)
    first_qu <- quantile(coverage_score, 0.25, na.rm = TRUE)
    median_val <- median(coverage_score, na.rm = TRUE)
    mean_val <- mean(coverage_score, na.rm = TRUE)
    third_qu <- quantile(coverage_score, 0.75, na.rm = TRUE)
    max_val <- max(coverage_score, na.rm = TRUE)

    # Create a dataframe with the statistics
    temp_df <- data.frame(
      Simulation = obj_name,
      Min = min_val,
      First_Qu = first_qu,
      Median = median_val,
      Mean = mean_val,
      Third_Qu = third_qu,
      Max = max_val
    )

    # Bind the results to the coverage_stats dataframe
    coverage_stats <- rbind(coverage_stats, temp_df)
  }
}


# Group the results by the names that match those in summary_by_simulation$Simulation
coverage_stats <- coverage_stats %>%
  mutate(Group = case_when(
    grepl(paste(substr(summary_by_simulation$Simulation, 1, nchar(summary_by_simulation$Simulation) - 11), collapse = "|"), substr(Simulation, 1, nchar(Simulation) - 11)) ~ substr(Simulation, 1, nchar(Simulation) - 11),
    TRUE ~ "Other"
  ))

# Create a box plot

# Update the Group column with new names
coverage_stats <- coverage_stats %>%
  mutate(Group = recode(Group, !!!group_names))

# Manually set the order of the groups
# coverage_stats$Group <- factor(coverage_stats$Group, levels = c(
#   "FW-Sys",
#   "FW-ZZ",
#   "QC-Sys",
#   "H-SG"
# ))
kable(coverage_stats)
# Save the transposed dataframe as a CSV file
output_path <- here("Output", "Simulation", "coverage_score.csv")
write.csv(coverage_stats, file = output_path, row.names = TRUE)