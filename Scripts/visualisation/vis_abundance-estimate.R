library(here)
library(sf)
library(dssd)
library(dsims)



# load processed data
# define list with WMU numbers
wmu_number_list <- c('501','503', '512', '517', '528')

# loop through WMU numbers
for(wmu_number in wmu_number_list){
  # Load processed data
  input_path <- here("Output", "Density", paste0("density",wmu_number,".RData"))
  load(file = input_path)

  # Assign individual names to data frames
  assign(paste0("density_",wmu_number), density)
  assign(paste0("total_abundance_",wmu_number), total_abundance)
  assign(paste0("region_",wmu_number), region)
  rm(density, total_abundance, region)
}
## Abundance Estimation
# Plot density surface
density_list <- paste0("density_",wmu_number)
total_abundance_list <- paste0("total_abundance_",wmu_number)


for (wmu_number in wmu_number_list){
  density <- get(paste0("density_",wmu_number))
  total_abundance <- get(paste0("total_abundance_",wmu_number))
  plot(density)
}

plot(density_501)

