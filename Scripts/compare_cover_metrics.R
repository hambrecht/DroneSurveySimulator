library(here)
library(knitr)
library(sf)
library(dplyr)
library(GGally)
library(dssd)

# Create a named vector with old and new group names
group_names <- c(
  "H_SG" = "H-SG",
  "Sys" = "Sys",
  "Rnd" = "Rnd",
  "ZZ" = "ZZ",
  "ZZC" = "ZZC",
  "FW_Sys" = "FW-Sys",
  "FW_ZZ" = "FW-ZZ",
  "QC_Sys" = "QC-Sys"
)

# List all files containing 'cover' in folder Output/Simulation
files <- list.files(path = here("Output", "Simulation"), pattern = "cover", full.names = TRUE)
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
  if (exists("short_zigzag_design")) rm(short_zigzag_design)
  if (exists("short_zigzag_transects")) rm(short_zigzag_transects)
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

# Remove all rows with 'Short-Zig' in the 'Simulation' column
merged_df <- merged_df %>%
  filter(Simulation != "Short-Zig")

kable(merged_df)

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

cover_score <- get.coverage(FW_Sys_design_528)
# hist(cover_score)
summary(cover_score)


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
coverage_stats$Group <- factor(coverage_stats$Group, levels = c(
  "H-SG",
  "Sys",
  "Rnd",
  "ZZ",
  "ZZC",
  "FW-Sys",
  "FW-ZZ",
  "QC-Sys"
))

# Define custom colors for each group
group_colors <- c(
  "H-SG" = "#b2df8a",   # Pink
  "Sys" = "##b2df8a",    # Light Blue
  "Rnd" = "#b2df8a",    # Green
  "ZZ" = "#b2df8a",     # Light Green
  "ZZC" = "#b2df8a",    # Blue
  "FW-Sys" = "#1f78b4", # Orange
  "FW-ZZ" = "#1f78b4",  # Purple
  "QC-Sys" = "#a6cee3"  # Dark Purple
)

ggplot(coverage_stats, aes(x = Group, y = Mean)) +
  geom_boxplot() +
  scale_fill_manual(values = group_colors) +
  labs(x = "Design", y = "Mean coverage score") +
  theme_minimal()
