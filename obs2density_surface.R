# Density surface modelling
# based on: https://distancesampling.org/R/vignettes/mexico-analysis.html

# Load necessary libraries for analysis and visualization
library(dsm)
library(ggplot2)
library(terra)
library(Distance)
library(sf)
library(dplyr)
library(tidyr)
library(units)
library(purrr)
# install.packages("C:/Users/lhambrec/Sync/1_Projects/Simulator/dsims_1.0.4.tar.gz", repos = NULL, type = "source")
library(dsims)

# Define constants
GRID_SIZE <- 500

# Function to load and transform spatial data
load_spatial_data <- function(moose_path, transects_path, sbfi_path, crs) {
  moose <- sf::st_read(moose_path) %>%
    sf::st_transform(crs = crs) %>%
    select(Latitude, Longitude, date, name)

  transects <- sf::st_read(transects_path, layer = "tracks") %>%
    sf::st_transform(crs = crs) %>%
    select(1)

  sbfi <- sf::st_read(sbfi_path) %>%
    sf::st_transform(crs = crs)

  list(moose = moose, transects = transects, sbfi = sbfi)
}

# Function to split transects into segments
split_into_segments <- function(linestring) {
  total_length <- st_length(linestring)
  total_length <- if (!inherits(total_length, "units")) set_units(total_length, "m") else total_length
  stopifnot(as.numeric(total_length) > 0)

  total_length_km <- set_units(total_length, "km")
  num_segments <- ceiling(as.numeric(total_length_km))
  stopifnot(num_segments > 0)

  equal_segment_length <- total_length_km / num_segments
  segment_lengths <- rep(equal_segment_length, num_segments)
  segment_lengths[num_segments] <- total_length_km - sum(segment_lengths[1:(num_segments - 1)])

  points <- sf::st_line_sample(linestring, sample = seq(0, 1, length.out = num_segments + 1)) %>%
    st_cast("POINT")
  stopifnot(length(points) >= 2)

  segments <- map2(points[-length(points)], points[-1], ~ {
    segment <- st_sfc(st_linestring(x = st_coordinates(c(.x, .y))), crs = st_crs(linestring))
    st_sf(geometry = segment)
  })

  segments_sf <- do.call(rbind, segments)
  stopifnot(nrow(segments_sf) >= 1)

  segments_sf
}

# Load data
data_paths <- list(
  moose_path = "D:\\WMU\\survey_data\\501_moose_locations.shp",
  transects_path = "D:\\WMU\\survey_data\\WMU 501 (2018-2019)\\WMU501_transects_2018.gpx",
  sbfi_path = "D:\\WMU\\base_data\\CA_Forest_Satellite_Based_Inventory_2020\\clipped\\sbfi_501.shp"
)
data <- load_spatial_data(data_paths$moose_path, data_paths$transects_path, data_paths$sbfi_path, crs = 3400)
moose <- data$moose
transects <- data$transects
sbfi <- data$sbfi

# Rename columns in sbfi
org_sbfi_names <- colnames(sbfi)
colnames(sbfi) <- c("OBJECTID", "ID", "TILE", "AREA_HA", "PERIMETER_M", "JURISDICTION", "ECOZONE", "ECOPROVINCE", "ECOREGION", "MANAGEMENT", "LC_WATER", "LC_SNOW_ICE", "LC_ROCK_RUBBLE", "LC_EXPOSED_BARREN", "LC_BRYOIDS", "LC_SHRUBS", "LC_WETLAND", "LC_WETLAND_TREED", "LC_HERBS", "LC_CONIFEROUS", "LC_BROADLEAF", "LC_MIXEDWOOD", "LC_TREED", "LC_FAO_FOREST", "LC_WETLAND_VEGETATION", "DISTURB_FIRE_PERC", "DISTURB_FIRE_YEAR", "DISTURB_FIRE_MAGNITUDE_MIN", "DISTURB_FIRE_MAGNITUDE_MAX", "DISTURB_FIRE_MAGNITUDE_AVG", "DISTURB_FIRE_MAGNITUDE_SD", "DISTURB_FIRE_MAGNITUDE_MEDIAN", "DISTURB_HARVEST_PERC", "DISTURB_HARVEST_YEAR", "RECOVERY_FIRE_MIN", "RECOVERY_FIRE_MAX", "RECOVERY_FIRE_AVG", "RECOVERY_FIRE_SD", "RECOVERY_FIRE_MEDIAN", "RECOVERY_HARVEST_MIN", "RECOVERY_HARVEST_MAX", "RECOVERY_HARVEST_AVG", "RECOVERY_HARVEST_SD", "RECOVERY_HARVEST_MEDIAN", "AGE_MIN", "AGE_MAX", "AGE_AVG", "AGE_SD", "AGE_MEDIAN", "AGE_0_10", "AGE_10_20", "AGE_20_30", "AGE_30_40", "AGE_40_50", "AGE_50_60", "AGE_60_70", "AGE_70_80", "AGE_80_90", "AGE_90_100", "AGE_100_110", "AGE_110_120", "AGE_120_130", "AGE_130_140", "AGE_140_150", "AGE_GT_150", "STRUCTURE_CANOPY_HEIGHT_MIN", "STRUCTURE_CANOPY_HEIGHT_MAX", "STRUCTURE_CANOPY_HEIGHT_AVG", "STRUCTURE_CANOPY_HEIGHT_SD", "STRUCTURE_CANOPY_HEIGHT_MEDIAN", "STRUCTURE_CANOPY_COVER_MIN", "STRUCTURE_CANOPY_COVER_MAX", "STRUCTURE_CANOPY_COVER_AVG", "STRUCTURE_CANOPY_COVER_SD", "STRUCTURE_CANOPY_COVER_MEDIAN", "STRUCTURE_LOREYS_HEIGHT_MIN", "STRUCTURE_LOREYS_HEIGHT_MAX", "STRUCTURE_LOREYS_HEIGHT_AVG", "STRUCTURE_LOREYS_HEIGHT_SD", "STRUCTURE_LOREYS_HEIGHT_MEDIAN", "STRUCTURE_BASAL_AREA_MIN", "STRUCTURE_BASAL_AREA_MAX", "STRUCTURE_BASAL_AREA_AVG", "STRUCTURE_BASAL_AREA_SD", "STRUCTURE_BASAL_AREA_MEDIAN", "STRUCTURE_BASAL_AREA_TOTAL", "STRUCTURE_AGB_MIN", "STRUCTURE_AGB_MAX", "STRUCTURE_AGB_AVG", "STRUCTURE_AGB_SD", "STRUCTURE_AGB_MEDIAN", "STRUCTURE_AGB_TOTAL", "STRUCTURE_VOLUME_MIN", "STRUCTURE_VOLUME_MAX", "STRUCTURE_VOLUME_AVG", "STRUCTURE_VOLUME_SD", "STRUCTURE_VOLUME_MEDIAN", "STRUCTURE_VOLUME_TOTAL", "SPECIES_NUMBER", "SPECIES_1", "SPECIES_1_PERC", "SPECIES_2", "SPECIES_2_PERC", "SPECIES_3", "SPECIES_3_PERC", "SPECIES_4", "SPECIES_4_PERC", "SPECIES_5", "SPECIES_5_PERC", "SPECIES_CONIFEROUS_PERC", "SPECIES_CML_1", "SPECIES_CML_1_PERC", "SPECIES_CML_2", "SPECIES_CML_2_PERC", "SPECIES_CML_3", "SPECIES_CML_3_PERC", "SPECIES_CML_4", "SPECIES_CML_4_PERC", "SPECIES_CML_5", "SPECIES_CML_5_PERC", "SPECIES_CML_CONIFEROUS_PERC", "SPECIES_CML_ASSEMBLAGES", "SPECIES_CML_ASSEMBLAGES_PERC", "SYMB_LAND_BASE_LEVEL", "SYMB_LAND_COVER_LEVEL", "SYMB_VEGETATION_LEVEL", "SYMB_DISTURBANCE", "SYMB_RECOVERY", "SYMB_AGE", "Shape_Length", "Shape_Area", "layer", "path", "geometry")

# Perform spatial join to get the polygon attributes for each point
joined <- st_join(moose, sbfi)
result <- joined[, c("STRUCTURE_CANOPY_HEIGHT_MEDIAN", "STRUCTURE_CANOPY_COVER_MEDIAN")]

# Replace NA values with 0
result <- result %>% replace_na(list(STRUCTURE_CANOPY_HEIGHT_MEDIAN = 0, STRUCTURE_CANOPY_COVER_MEDIAN = 0))

# Add the results to the original points sf object
moose$canopy_height <- result$STRUCTURE_CANOPY_HEIGHT_MEDIAN
moose$canopy_cover <- result$STRUCTURE_CANOPY_COVER_MEDIAN


# Split transects into segments
transects_segments <- transects %>%
  st_cast("LINESTRING") %>%
  st_geometry() %>%
  map(split_into_segments) %>%
  bind_rows() %>%
  st_sf()

st_crs(transects_segments) <- st_crs(transects)
transects_segments$Sample.Label <- row_number(transects_segments)

# Calculate distances between moose points and transect segments
distances <- st_distance(moose, transects_segments, tolerance = set_units(600, "m"))

# Find closest segment for each moose point
closest_segments <- data.frame(
  moose_id = integer(nrow(moose)),
  Sample.Label = integer(nrow(moose)),
  distance = numeric(nrow(moose))
)

for (i in 1:nrow(moose)) {
  min_distance_index <- which.min(distances[i, ])
  closest_segments$moose_id[i] <- i
  closest_segments$Sample.Label[i] <- transects_segments$Sample.Label[min_distance_index]
  closest_segments$distance[i] <- distances[i, min_distance_index]
}

# Merge closest segments with moose data
moose <- merge(moose, closest_segments, by.x = "row.names", by.y = "moose_id")
colnames(moose)[1] <- "object"
colnames(moose)[5] <- "Transect.Label"
moose$size <- 1
moose$Effort <- 10
moose$distance <- moose$distance / 1000

# Prepare data for density surface modelling
segdata <- as.data.frame(sf::st_drop_geometry(moose[, c("Latitude", "Longitude", "Effort", "Transect.Label", "Sample.Label", "canopy_height", "canopy_cover")]))
distdata <- as.data.frame(sf::st_drop_geometry(moose[, c("object", "Latitude", "Longitude", "distance", "Effort", "size", "canopy_height", "canopy_cover")]))
distdata$detected <- 1
segdata$x <- distdata$x <- sf::st_coordinates(moose)[, 1]
segdata$y <- distdata$y <- sf::st_coordinates(moose)[, 2]
obsdata <- as.data.frame(sf::st_drop_geometry(moose[, c("object", "distance", "Effort", "Sample.Label", "size")]))

# Load WMU outline and create prediction grid
wmu <- sf::st_read("D:\\WMU\\base_data\\WMU\\wmu_501_3400.shp")
# Select only the OBJECTID column
wmu <- wmu[, "OBJECTID"]
# Create the survey region
region <- make.region(
  region.name = "study area",
  shape = wmu
)
plot(region)

## # Create the density surface
density <- dsims::make.density(
  region = region,
  x.space = GRID_SIZE,
)
coords <- sf::st_drop_geometry(density@density.surface[[1]][, c("x", "y")])

# Function to create plot objects for predictions
grid_plot_obj <- function(fill, name, sf_obj) {
  data <- data.frame(fill = fill)
  names(data) <- name
  combined_sf <- sf_obj %>%
    mutate(id = row_number()) %>%
    left_join(data, by = c("id" = "row.names"))

  ggplot(combined_sf) +
    geom_sf(aes_string(fill = name)) +
    theme_minimal()
}

# Exploratory data analysis and density surface modelling
# Null model
detfc.hr.null <- ds(distdata, max(distdata$distance), key = "hr", adjustment = NULL)
summary(detfc.hr.null)
par(mfrow = c(1, 2))
plot(detfc.hr.null, showpoints = FALSE, pl.den = 0, lwd = 2)
ddf.gof(detfc.hr.null$ddf)

# Canopy height model
detfc.hr.height <- ds(distdata, max(distdata$distance),
  formula = ~ as.factor(canopy_height),
  key = "hr", adjustment = NULL
)
summary(detfc.hr.height)
par(mfrow = c(1, 2))
plot(detfc.hr.height, showpoints = FALSE, pl.den = 0, lwd = 2)
ddf.gof(detfc.hr.height$ddf)

# Canopy cover model
detfc.hr.cover <- ds(distdata, max(distdata$distance),
  formula = ~ as.factor(canopy_cover),
  key = "hr", adjustment = NULL
)
summary(detfc.hr.cover)
par(mfrow = c(1, 2))
plot(detfc.hr.cover, showpoints = FALSE, pl.den = 0, lwd = 2)
ddf.gof(detfc.hr.cover$ddf)

par(mfrow = c(1, 1))


# Fit a DSM
dsm.xy <- dsm(count ~ s(x, y), detfc.hr.null, segdata, obsdata, method = "REML")
summary(dsm.xy)
vis.gam(dsm.xy, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = GRID_SIZE)
dsm.xy.height <- dsm(count ~ s(x, y, k = 10) + s(height, k = 20), detfc.hr.null, segdata, obsdata, method = "REML")
summary(dsm.xy.height)
plot(dsm.xy.height, select = 2)
vis.gam(dsm.xy.height, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = GRID_SIZE)
dsm.xy.cover <- dsm(count ~ s(x, y, k = 10) + s(cover, k = 20), detfc.hr.null, segdata, obsdata, method = "REML")
summary(dsm.xy.cover)
plot(dsm.xy.cover, select = 2)
vis.gam(dsm.xy.cover, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = GRID_SIZE)

# Spatial models when there are covariates in the detection function
dsm.est.xy_height <- dsm(abundance.est ~ s(x, y), detfc.hr.height, segdata, obsdata, method = "REML")
vis.gam(dsm.est.xy_height, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = GRID_SIZE)
dsm.est.xy_cover <- dsm(abundance.est ~ s(x, y), detfc.hr.cover, segdata, obsdata, method = "REML")
vis.gam(dsm.est.xy_cover, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = GRID_SIZE)


dsm.xy.tweedie <- dsm(count ~ s(x, y), detfc.hr.null, segdata, obsdata, family = tw(), method = "REML")
summary(dsm.xy.tweedie)

# model checking
gam.check(dsm.xy)
gam.check(dsm.xy.height)
gam.check(dsm.xy.cover)
gam.check(dsm.est.xy_height)
gam.check(dsm.est.xy_cover)

## randomised quantile residuals
rqgam.check(dsm.xy)
rqgam.check(dsm.xy.height)
rqgam.check(dsm.xy.cover)
rqgam.check(dsm.est.xy_height)
rqgam.check(dsm.est.xy_cover)

# gam.check(dsm.xy)
# rqgam_check(dsm.xy)
# vis_concurvity(dsm.xy)

# Check for autocorrelation
dsm_cor(dsm.xy, max.lag = 10, Segment.Label = "Sample.Label")
dsm_cor(dsm.xy.height, max.lag = 10, Segment.Label = "Sample.Label")
dsm_cor(dsm.xy.cover, max.lag = 10, Segment.Label = "Sample.Label")
dsm_cor(dsm.est.xy_height, max.lag = 10, Segment.Label = "Sample.Label")
dsm_cor(dsm.est.xy_cover, max.lag = 10, Segment.Label = "Sample.Label")

# Model selection
# make a data.frame to print out
mod_results <- data.frame(
  "Model name" = c(
    "`dsm.xy`", "`dsm.xy.height`", "`dsm.xy.cover`", "`dsm.est.xy_height`",
    "`dsm.est.xy_cover`"
  ),
  "Description" = c(
    "Bivariate smooth of location, quasipoisson",
    "Bivariate smooth of location, smooth of height, quasipoisson",
    "Bivariate smooth of location, smooth of cover, quasipoisson",
    "Bivariate smooth of location, Tweedie, height covariate in detection function",
    "Bivariate smooth of location, Tweedie, cover covariate in detection function"
  ),
  "Deviance explained" = c(unlist(lapply(
    list(
      dsm.xy,
      dsm.xy.depth,
      dsm.xy.tweedie,
      dsm.xy.tweedie.soap,
      dsm.est.xy
    ),
    function(x) {
      paste0(round(summary(x)$dev.expl * 100, 2), "%")
    }
  )), NA)
)

kable(mod_results, col.names = c("Model name", "Description", "Deviance explained"))

# Abundance estimation
dsm.xy.pred <- predict(dsm.xy, coords, dsm.xy$offset[1])
# Calculate total abundance over the survey area
sum(dsm.xy.pred)

# fill density object with predeicted values
density@density.surface[[1]]$density <- dsm.xy.pred

# plot density
plot(density)
