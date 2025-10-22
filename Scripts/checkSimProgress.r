# Load necessary libraries
library(here)
library(dsims)
library(knitr)
library(RColorBrewer)
library(dplyr)
library(sf)
library(units)

# List all objects containing 'QC_Sys_design_'
input_path <- here("Output", "Simulation", paste0("designsQC.RData"))
load(input_path)


ABUNDANCE_LIST <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90)

loaded_objects <- ls(pattern = "^QC_")

# Initialize counters and storage
existing_files_count <- 0
total_files_count <- 0
missing_files <- c()
missing_designs <- c()

for (design_name in loaded_objects) {
  design <- get(design_name)
  
  for (ABUNDANCE in ABUNDANCE_LIST) {
    total_files_count <- total_files_count + 1
    output_path <- here("Output", "Simulation", paste0(design_name, "-density_sim-A", ABUNDANCE, "_adjust.RData"))
    
    if (file.exists(output_path)) {
      existing_files_count <- existing_files_count + 1
    #   print(paste("File exists:", output_path))
    } else {
      missing_files <- c(missing_files, output_path)
      missing_designs <- c(missing_designs, design_name)
      print(paste("Missing:", design_name))
    }
  }
}

# Summary output
print(paste("Total potential output files:", total_files_count))
print(paste("Existing files:", existing_files_count))
print(paste("Missing files:", total_files_count - existing_files_count))

# Optional detailed lists
cat("\nUnique designs with missing files:\n")
print(unique(missing_designs))

# cat("\nList of missing file paths:\n")
# print(missing_files)
