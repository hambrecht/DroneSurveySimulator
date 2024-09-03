# load processed data
# define list with WUM numbers
wmu_number_list <- c('501','503', '512', '517', '528')

# loop through WUM numbers
for(wmu_number in wmu_number_list){
  # Load processed data
  input_path <- here("Output", "Density", paste0("density",wum_number,".RData"))
  load(file = input_path)

  # Assign individal names to data frames
  assign(paste0("density_",wmu_number), density)
  assign(paste0("total_abundance_",wmu_number), total_abundance)
  assign(paste0("region_",wmu_number), region)
  rm(density, total_abundance, region)
}
## Abundance Estimation
# Plot density surface
for (wmu_number in wmu_number_list){
  plot(paste0("density_",wmu_number))
  print(paste0("total_abundance_",wmu_number))
}




