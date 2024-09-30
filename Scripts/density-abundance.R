library(here)
library(dsims)

# Define grid cell size (in meters)
GRID_SIZE <- 500
wmu_number_list <- c('501', '503', '512', '517', '528')
wmu_number <- wmu_number_list[1]
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

# Predict abundance using the null model
# Find the index of the sublist with name "Second"
index <- which(sapply(dsm_list, function(x) x$name) == "dsm_xy_detfc_hn_cos_trunc_0.43_bin_0.1")

# Access the sublist

sel_model <- dsm_list[[index]]
dsm_xy_pred <- predict(sel_model, coords, rep(sel_model$offset[1], length(coords[[1]])))
dsm_var_gam(sel_model, coords, rep(sel_model$offset[1], length(coords[[1]])))

# Calculate total abundance over the survey area
total_abundance <- sum(dsm_xy_pred)

# Update density object with predicted values
density@density.surface[[1]]$density <- dsm_xy_pred

# Plot density surface
plot(density)
plot(density@density.surface[[1]]['density'])


# Save processed data
output_path <- here("Output", "Density", paste0("density",wmu_number,".RData"))
save(density, total_abundance, region, file = output_path)