library(here)
library(ggplot2)
library(RColorBrewer)

# Ensure the output directory exists
dir.create(here("Data", "Plot"), showWarnings = FALSE)

# Define custom functions

#' Check for Required Columns in a Dataframe
#'
#' Checks if all specified columns are present in a dataframe. The check is case-insensitive.
#'
#' @param df A dataframe to check.
#' @param required_cols A character vector of column names that must be present in the dataframe.
#'
#' @return Throws an error if any required columns are missing.
#'
#' @examples
#' df <- data.frame(A = 1:5, b = 6:10, D = 11:15, E = 16:20)
#' required_columns <- c("A", "B", "C")
#' check_columns_present(df, required_columns)
#'
#' @seealso \link[base]{stop}
#' @export
check_columns_present <- function(df, required_cols) {
  # Convert column names to lowercase
  actual_cols <- tolower(colnames(df))
  required_cols_lower <- tolower(required_cols)
  
  # Identify missing columns
  missing_cols <- setdiff(required_cols_lower, actual_cols)
  
  # Stop execution if there are missing columns
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
}

# Function to extract the sf object from a custom class
extract_sf_from_region <- function(region_obj) {
  # Assuming the sf object is stored in a slot named 'region'
  if ("sf" %in% class(region_obj@region)) {
    return(region_obj@region)
  } else {
    stop("The 'region' slot does not contain an sf object.")
  }
}


# Define a colour palette
colour_palette <- brewer.pal(5, "RdBu")

# define list with WUM numbers
wmu_number_list <- c('501','503', '512', '517', '528')

# Define required columns for dataframes
segdata_required <- c("longitude", "latitude", "Transect.Label", "Sample.Label", "x", "y", "Effort")
distdata_required <- c("object", "size", "longitude", "latitude", "x", "y", "Effort")
obsdata_required <- c("object", "Sample.Label", "size", "distance", "Effort")

# loop through WUM numbers
for(wmu_number in wmu_number_list){
  # Load processed data
  input_path <- here("Output", "PrepData", paste0("prepared",wmu_number,".RData"))
  load(file = input_path)
  
  # Check presence of required columns
  check_columns_present(segdata, segdata_required)
  check_columns_present(distdata, distdata_required)
  check_columns_present(obsdata, obsdata_required)
  
  # Assign individal names to data frames
  assign(paste0("segdata",wmu_number), segdata)
  assign(paste0("distdata",wmu_number), distdata)
  assign(paste0("obsdata",wmu_number), obsdata)
  rm(segdata, distdata, obsdata)
}

## Exploratory data analysis
# Segmentation data
for(wmu_number in wmu_number_list){
  print(paste0("Segmentation data for WUM number: ", wmu_number))
  print(summary(get(paste0("segdata",wmu_number))))
}

# Observation data
for(wmu_number in wmu_number_list){
  print(paste0("Observation data for WUM number: ", wmu_number))
  print(summary(get(paste0("obsdata",wmu_number))))
}

# Distance data
for(wmu_number in wmu_number_list){
  print(paste0("Distance data for WUM number: ", wmu_number))
  print(summary(get(paste0("distdata",wmu_number))))
}

## Histogram
# Define the list of variables to plot
variables <- c("distance", "canopy_height", "canopy_cover", "agb", "vol")

# Loop through the WMU numbers and plot/save histograms
for (wmu_number in wmu_number_list) {
  # Retrieve the dataset for the current WMU number
  distdata <- get(paste0("distdata", wmu_number))
  
  # List of variables to plot
  variables <- c("distance", "canopy_height", "canopy_cover", "agb", "vol")
  
  # Loop through each variable
  for (variable in variables) {
    # Create histogram plot
    p <- ggplot(distdata, aes_string(x = variable)) +
      geom_histogram(bins = 30, fill = colour_palette[5], colour = colour_palette[3]) +
      labs(title = paste(wmu_number, variable), x = variable, y = "Frequency") +
      theme_minimal() +
      theme(text = element_text(size = 12),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 16, face = "bold"))
    
    print(p)
    # Save the plot
    output_file <- here("Output", "Plots", paste0("histogram_", wmu_number, "_", variable, ".png"))
    ggsave(output_file, plot = p, width = 8, height = 6)
  }
}
