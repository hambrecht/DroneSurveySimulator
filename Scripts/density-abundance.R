library(here)
library(dsims)
library(dsm) # for dsm_var_gam function
library(ggplot2)

# Define grid cell size (in meters)
GRID_SIZE <- 500
wmu_number_list <- c('501', '503', '512', '517', '528')
wmu_number <- wmu_number_list[2]
# save DSM list
input_path <- here("Output", "DSM", paste0("dsm",wmu_number,".RData"))
load(file = input_path)


## Density Surface Template

# Create a density object for coordinate retrieval and abundance prediction
# Define survey region
region <- make.region(
  region.name = "study area",
  shape = wmu
)

# Create density surface
density <- dsims::make.density(
  region = region,
  x.space = GRID_SIZE
)
# Extract coordinates
coords <- sf::st_drop_geometry(density@density.surface[[1]][, c("x", "y")])
coords$offset <- 0.5^2

## Abundance Estimation
# Select the best model based on the output of vis-detfc-dsm script. Often the model with the lowest AIC value is selected.
# in this case the model with hazard rate key and no adjustment is selected and 99.5% trunctaion distance is used.
# Get indices of entries with 'hr_null' in the name
# indices_hr_null <- grep("hr_null", dsm_list)
# # Filter out indices that also have 'trunc_0.6' in the name
# indices_final <- indices_hr_null[!grepl("trunc_0.6|cc", dsm_list[indices_hr_null])]
# indices_final <- indices_hr_null[!grepl("trunc_0.55|cc", dsm_list[indices_hr_null])]
indices_hn_cos_null <- grep("dsm_xy_detfc_hn_cos_trunc_0.6_bin_0.05", dsm_list)
names(dsm_list)

sel_model <- dsm_list[indices_hn_cos_null][[1]]

# par(mfrow = c(1, 2))
# sel_model <- dsm_list[6][[1]]
summary(sel_model)
vis.gam(sel_model, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = 500)
dsm_xy_pred <- predict(sel_model, coords, coords$offset)
model_variance <- dsm_var_gam(sel_model, coords, rep(sel_model$offset[1], length(coords[[1]])))


dsm_xy_pred <- dsm_xy_pred * 0.1809381
# Calculate total abundance over the survey area
total_abundance <- sum(dsm_xy_pred)
sum(dsm_xy_pred)
mean(dsm_xy_pred)*4

# Update density object with predicted values
density@density.surface[[1]]$density <- dsm_xy_pred



# Save processed data
output_path <- here("Output", "Density", paste0("density",wmu_number,".RData"))
save(density, total_abundance, region, wmu, model_variance, model_name, file = output_path)

# test
plot(density)
