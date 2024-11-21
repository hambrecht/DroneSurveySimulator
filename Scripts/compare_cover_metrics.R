library(here)
library(knitr)
library(sf)
library(dplyr)
library(GGally)

# List all files containing 'cover' in folder Output/Simulation
files <- list.files(path = here("Output", "Simulation"), pattern = "cover", full.names = TRUE)
wmu_path_list <- list.files(path = here("Output", "PrepData"), pattern = "prepared", full.names = TRUE)
wmu_path_list <- wmu_path_list[-4]

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

# Assuming the WMU data is already loaded and renamed
# Create an empty dataframe to store the results
wmu_metrics <- data.frame()

# Loop through each renamed WMU object and calculate the metrics
for (i in seq_along(file_ids)) {
    wmu_id <- file_ids[[i]]
    wmu_data <- get(paste0("wmu_", wmu_id))
    
    # Calculate area, perimeter, and shape index
    wmu_data <- wmu_data %>%
        mutate(
            area = st_area(geometry),
            perimeter = st_length(st_cast(geometry, "MULTILINESTRING")),
            shapeindex = perimeter / (2 * sqrt(pi * area))
        )
    
    # Summarize the metrics
    wmu_summary <- wmu_data %>%
        summarise(
            wmu_id = wmu_id,
            total_area = sum(as.numeric(area)),
            total_perimeter = sum(as.numeric(perimeter)),
            mean_shapeindex = mean(as.numeric(shapeindex))
        ) %>%
        as.data.frame()
    
    # Bind the summary to the results dataframe
    wmu_metrics <- bind_rows(wmu_metrics, wmu_summary[1:4])
}

# Print the resulting dataframe
print(wmu_metrics)
wmu_metrics$total_area[1]

# Merge all objects into one dataframe and add `i` as a column to identify the origins
merged_df <- data.frame()

for (i in 1:length(file_ids)) {
    obj_name <- paste0("design_comparison_df_", file_ids[i])
    if (exists(obj_name)) {
        temp_df <- get(obj_name)
        temp_df$Simulation <- paste0(temp_df$Simulation, "_", file_ids[i])
        temp_df$total_area = wmu_metrics$total_area[i]
        temp_df$total_perimeter = wmu_metrics$total_perimeter[i]
        temp_df$mean_shapeindex = wmu_metrics$mean_shapeindex[i]
        merged_df <- bind_rows(merged_df, temp_df)
    }
}

# Sort merged_df by the `Simulation` column
merged_df <- merged_df %>% arrange(Simulation)

# Print the merged dataframe
kable(merged_df)





# Select the relevant columns
data <- merged_df %>%
  select(mean_shapeindex, ends_with("Percentage"), Design)

# Create the scatter plot matrix
ggpairs(data, aes(color = Design), 
        columns = 1:(ncol(data)-1), 
        upper = list(continuous = "points"),
        lower = list(continuous = "points"),
        diag = list(continuous = "densityDiag"))

# Install and load the necessary packages
library(ggplot2)

# Create the scatter plot with regression lines
ggplot(merged_df, aes(x = total_perimeter, y = Off_Effort_Return_Percentage, color = Design)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Influence of mean_shapeindex on Off_Effort_Return_Percentage",
       x = "Mean Shape Index",
       y = "Off Effort Return Percentage") +
  theme_minimal()

# Install and load the necessary packages
install.packages("broom")
library(dplyr)
library(broom)

# Perform linear regression with interaction term
interaction_model <- lm(Off_Effort_Return_Percentage ~ total_perimeter * Design, data = merged_df)

# Summarize the regression results
summary(interaction_model)

# Tidy the results for better readability
tidy_results <- tidy(interaction_model)
print(tidy_results)
