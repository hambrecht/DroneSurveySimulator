


# Re-producable example
# load libraries
library(dsm)
library(dsims)
library(Distance)
library(sf)
library(units)

# load example data
data(mexdolphins)

# load prediction grid
pred_grid <- sf::st_read("prediction-grid.shp")
survey_area <- sf::st_read("survey-area.shp")

# Extract geometries
geometries <- st_geometry(pred_grid)

# Initialize vectors to store widths and heights
widths <- numeric(length(geometries))
heights <- numeric(length(geometries))

# Calculate the bounding box for each polygon and compute width and height
for (i in seq_along(geometries)) {
  bbox <- st_bbox(geometries[[i]])
  widths[i] <- bbox$xmax - bbox$xmin
  heights[i] <- bbox$ymax - bbox$ymin
}

# Calculate the average width and height
average_width <- mean(widths)
average_height <- mean(heights)

# Convert SpatVector to sf object
pred_grid_sf <- st_as_sf(pred_grid)

# Estimating the detection function
detfc.hr.null <- ds(distdata, max(distdata$distance), key = "hr", adjustment = NULL)

# Fitting DSM
dsm.xy <- dsm(count ~ s(x, y), detfc.hr.null, segdata, obsdata, method = "REML")

# Abundance estimation
dsm.xy.pred <- predict(dsm.xy, pred_grid_sf, 0)

# Assign prediction to grid
pred_grid_sf$predictions <- dsm.xy.pred

# Calculate centroids
centroids <- st_centroid(pred_grid_sf)
# Extract x and y coordinates from centroids
centroids_df <- st_coordinates(centroids)

# Create density surface
density_surface <- data.frame(x = centroids_df[, 1], y = centroids_df[, 2], density = pred_grid_sf$predictions)

# Create list, as required by make.density() function
density_surface_ls <- list(density_surface)

# Create the survey region
region <- make.region(
    region.name = "survey area",
    shape = survey_area
)

## Create the density with pre-defined density surface
pre_ds_density <- dsims::make.density(
    region = region,
    x.space = average_width,
    y.space = average_height,
    density.surface = density_surface_ls
)
# Error: All strata must have some cells with non-zero density. Check that you have correctly specified your density grid. Large grid spacing may also generate this error.

# Create density = 1
density_1 <- dsims::make.density(
    region = region,
    x.space = average_width,
    y.space = average_height
)

# Get density surface from density = 1
reverse_density_surface <- get.densities(density_1, coords = TRUE)

# Create the density from density surface retrieved from 1 density
rev_density <- density <- dsims::make.density(
    region = region,
    x.space = average_width,
    y.space = average_height,
    density.surface = reverse_density_surface
)
# Error: All strata must have some cells with non-zero density. Check that you have correctly specified your density grid. Large grid spacing may also generate this error.


# Use reverse density surface for prediction
dsm.xy.pred <- predict(dsm.xy, reverse_density_surface, 0)
reverse_density_surface$predictions <- dsm.xy.pred
rev_rev_density <- dsims::make.density(
    region = region,
    x.space = average_width,
    y.space = average_height,
    density.surface = reverse_density_surface
)
# Error: All strata must have some cells with non-zero density. Check that you have correctly specified your density grid. Large grid spacing may also generate this error.
