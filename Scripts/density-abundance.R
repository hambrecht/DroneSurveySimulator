library(here)
library(dsims)
library(dsm) # for dsm_var_gam function

# Define grid cell size (in meters)
GRID_SIZE <- 500
wmu_number_list <- c('501', '503', '512', '517', '528')
wmu_number <- wmu_number_list[4]
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

## Abundance Estimation
# Select the best model based on the output of vis-detfc-dsm script. Often the model with the lowest AIC value is selected.
# in this case the model with hazard rate key and no adjustment is selected and 99.5% trunctaion distance is used.
# Get indices of entries with 'hr_null' in the name
indices_hr_null <- grep("hr_null", dsm_list)
# Filter out indices that also have 'trunc_0.6' in the name
indices_final <- indices_hr_null[!grepl("trunc_0.6|cc", dsm_list[indices_hr_null])]
# indices_final <- indices_hr_null[!grepl("trunc_0.55|cc", dsm_list[indices_hr_null])]
model_name <- names(dsm_list[indices_final])
print(paste0("The model that is use is: ",model_name))
sel_model <- dsm_list[indices_final][[1]]
par(mfrow = c(1, 2))
sel_model <- dsm_list[5][[1]]
summary(sel_model)
vis.gam(sel_model, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = 500)
dsm_xy_pred <- predict(sel_model, coords, rep(sel_model$offset[1], length(coords[[1]])))
model_variance <- dsm_var_gam(sel_model, coords, rep(sel_model$offset[1], length(coords[[1]])))

# Calculate total abundance over the survey area
total_abundance <- sum(dsm_xy_pred)

# Update density object with predicted values
density@density.surface[[1]]$density <- dsm_xy_pred



# Save processed data
output_path <- here("Output", "Density", paste0("density",wmu_number,".RData"))
save(density, total_abundance, region, wmu, model_variance, model_name, file = output_path)

# test
plot(density)
