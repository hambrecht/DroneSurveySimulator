#' ---
#' title: "Covariate Extrapolation and Density Surface Modeling"
#' author: "Leonard Hambrecht"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---
#'
#' ## Overview
#' This script performs covariate extrapolation and density surface modeling for ecological data. 
#' It includes data preparation, raster manipulation, probability calculations, and distance sampling analysis.
#'
#' ## Libraries
#' The script uses the following libraries:
#' - `here`: For constructing file paths.
#' - `dsims`: For simulating distance sampling surveys.
#' - `raster`: For raster data manipulation.
#' - `sf`: For handling spatial data.
#' - `units`: For unit conversions in spatial data.
#' - `Distance`: For distance sampling analysis.
#' - `dplyr`: For data manipulation.
#' - `purrr`: For functional programming.
#' - `dsmextra`: For extrapolation metrics in ecological models.
#'
#' ## Functions
#' - `rescale_raster`: Rescales raster values between 0 and 1.
#' - `calculate_probability`: Calculates probabilities based on raster values or distances using hazard-rate decline functions.
#' - `split_into_segments`: Splits LINESTRING geometries into equal-length segments.
#' - `check_columns_present`: Ensures required columns are present in a dataframe.
#'
#' ## Data Loading
#' - Loads density and cover data for specified Wildlife Management Units (WMUs).
#' - Sets coordinate reference systems (CRS) for spatial data.
#'
#' ## Population Simulation
#' - Simulates a population based on density data using `dsims`.
#' - Generates example population and extracts individual locations.
#'
#' ## Raster Manipulation
#' - Creates raster layers with skewed random data.
#' - Applies smoothing filters and rescales raster layers.
#' - Stacks raster layers for further analysis.
#'
#' ## Covariate Extraction
#' - Extracts raster values for spatial points and adds them as attributes.
#' - Filters points based on calculated probabilities.
#'
#' ## Distance Sampling Preparation
#' - Clips points to the study region.
#' - Calculates distances between points and transect segments.
#' - Prepares data for distance sampling analysis, including segment, distance, and observation data.
#'
#' ## Extrapolation Analysis
#' - Computes extrapolation metrics using `dsmextra`.
#' - Calculates Gower's distances and nearby metrics.
#' - Maps extrapolation results and covariate comparisons.
#'
#' ## Density Surface Modeling
#' - Defines multiple density surface models with varying smoothness parameters and covariates.
#' - Evaluates model performance using deviance explained and diagnostic checks.
#'
#' ## Prediction
#' - Predicts density surfaces using the fitted models.
#' - Visualizes predicted density surfaces as rasters.
#'
#' ## Notes
#' - Ensure all required columns are present in dataframes before analysis.
#' - Adjust smoothing parameters (`k`) in density surface models for optimal performance.
#' - Validate raster data before prediction to avoid missing columns.
#'
#' ## Outputs
#' - Histograms of distances before and after filtering.
#' - Maps of extrapolation metrics and density surfaces.
#' - Predicted density values for the study area.
#'
#' ## References
#' - `dsims` and `dsm` documentation for population simulation and density surface modeling.
#' - `dsmextra` for extrapolation metrics and mapping.
#'
#' ## Disclaimer
#' This script assumes the availability of input data and appropriate CRS settings. Ensure all dependencies are installed and data paths are correct.

# Load libraries ====
library(here)
library(dsims)
library(raster)
library(sf)
library(units)   # Units for spatial data
library(Distance)
library(dsm)
library(dplyr)
library(purrr)   # Functional programming tools
library(dsmextra)     # Extrapolation toolkit for ecological models

# Define functions ====
# Function to rescale raster values between 0 and 1
rescale_raster <- function(raster_layer) {
  min_val <- raster::cellStats(raster_layer, stat = "min", na.rm = TRUE)
  max_val <- cellStats(raster_layer, stat = "max", na.rm = TRUE)
  (raster_layer - min_val) / (max_val - min_val)
}

# Define a combined function to calculate the probability
calculate_probability <- function(input, type = "value", lambda = 1, k = 1) {
  if (is.na(input)) {
    return(NA)  # Handle NA values
  }
  if (type == "value") {
    if (input == 0) {
      return(1)  # Keep all points with layer1 == 0
    } else if (input == 1) {
      return(0)  # Remove all points with layer1 == 1
    } else {
      return(1 - input)  # Gradually decrease probability as value increases
    }
  } else if (type == "distance") {
    # Hazard-rate decline function
    return(exp(-lambda * input^k))
  } else {
    stop("Invalid type. Use 'value' or 'distance'.")
  }
}

split_into_segments <- function(linestring) {
  total_length <- st_length(linestring)
  total_length <- if (!inherits(total_length, "units")) units::set_units(total_length, "m") else total_length

  num_segments <- ceiling(as.numeric(set_units(total_length, "km")))
  
  if (num_segments <= 0) {
    stop("The number of segments must be greater than 0. Check the input linestring.")
  }
  
  equal_segment_length <- set_units(total_length, "km") / num_segments
  segment_lengths <- rep(equal_segment_length, num_segments)
  segment_lengths[num_segments] <- set_units(total_length, "km") - sum(segment_lengths[1:(num_segments - 1)])

  points <- sf::st_line_sample(linestring, sample = seq(0, 1, length.out = num_segments + 1)) %>%
    st_cast("POINT")

  segments <- map2(points[-length(points)], points[-1], ~{
    segment <- st_sfc(st_linestring(x = st_coordinates(c(.x, .y))), crs = st_crs(linestring))
    st_sf(geometry = segment)
  })

  do.call(rbind, segments)
}

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


# Load data ====
wmu_number_list <- c("501", "503", "512", "517", "528")
wmu_number <- wmu_number_list[1]
input_path <- here("Output", "Density", paste0("density", wmu_number, ".RData"))
load(file = input_path)
input_path <- here("Output", "Simulation", paste0("cover-WMU", wmu_number, ".RData"))
load(file = input_path)

# Set variables ----
CRS <- st_crs(wmu)
CRSsp <- sp::CRS(CRS$wkt)

# Create example population based on density data ====
# Create population description
pop_desc <- make.population.description(
  region = region,
  density = density,
  N = total_abundance,
  fixed.N = T
)
detect_FWG <- make.detectability(
  key.function = "hn",
  scale.param = 170,
  # shape.param = 3,
  truncation = 260
)
plot(detect_FWG, pop_desc, legend = FALSE)
example_population <- generate.population(object = pop_desc, detectability = detect_FWG, region = region)

# Extract location of individuals
points_sf <- st_as_sf(example_population@population, coords = c("x", "y"), crs = CRS)
points_sp <- as(points_sf, "Spatial")

# Create covariate data ----
# Assuming 'region' is an sf object
region_extent <- st_bbox(region@region)
resolution <- 1000  # Define the resolution of the raster

# Create a raster template
raster_template <- raster(
  xmn = region_extent["xmin"],
  xmx = region_extent["xmax"],
  ymn = region_extent["ymin"],
  ymx = region_extent["ymax"],
  res = resolution,
  crs = CRS
)


# Generate three raster layers and fill them with skewed random data
raster_layer1 <- raster_template
values(raster_layer1) <- rbeta(ncell(raster_layer1), shape1 = 2, shape2 = 5)  # Skewed toward 0
raster_layer2 <- raster_template
values(raster_layer2) <- rbeta(ncell(raster_layer2), shape1 = 2, shape2 = 5)  # Skewed toward 0
raster_layer3 <- raster_template
values(raster_layer3) <- rbeta(ncell(raster_layer3), shape1 = 2, shape2 = 5)  # Skewed toward 0

# Define a smoothing filter (5x5 mean filter)
smoothing_filter <- matrix(1 / 9, nrow = 5, ncol = 5)

# Apply the smoothing filter to all raster layers in one line
raster_layer1_smooth <- focal(raster_layer1, w = smoothing_filter, fun = sum, na.rm = TRUE)
raster_layer2_smooth <- focal(raster_layer2, w = smoothing_filter, fun = sum, na.rm = TRUE)
raster_layer3_smooth <- focal(raster_layer3, w = smoothing_filter, fun = sum, na.rm = TRUE)

# Rescale the smoothed raster layers
# Apply the smoothing filter to all raster layers in one line
raster_layer1_rescaled <- rescale_raster(raster_layer1_smooth)
raster_layer2_rescaled <- rescale_raster(raster_layer2_smooth)
raster_layer3_rescaled <- rescale_raster(raster_layer3_smooth)

# create raster stack
raster_stack <- stack(raster_layer1_rescaled, raster_layer2_rescaled, raster_layer3_rescaled, crs = CRS)

# Extract values from raster layers for the points in points_sf
values_layer1 <- extract(raster_layer1_rescaled, points_sp)
values_layer2 <- extract(raster_layer2_rescaled, points_sp)
values_layer3 <- extract(raster_layer3_rescaled, points_sp)

# Add the extracted values as new columns to points_sf
points_sp$layer1 <- values_layer1
points_sp$layer2 <- values_layer2
points_sp$layer3 <- values_layer3


# Convert points_sp back to sf object
points_sf <- st_as_sf(points_sp)

# Histogram of distances before filtering
hist(points_sf$layer1, breaks = 20, col = "blue", main = "Histogram Before Filtering",
     xlab = "Layer1", ylab = "Frequency")

## Create point distrubtion ----
# Apply the function to calculate probabilities for each point
points_sf$keep_probability <- sapply(points_sf$layer1, calculate_probability, type = "value")

# Filter points based on the calculated probabilities
points_sf <- points_sf[runif(nrow(points_sf)) < points_sf$keep_probability,]

# Drop the temporary column
points_sf$keep_probability <- NULL

# Histogram of distances before filtering
hist(points_sf$layer1, breaks = 20, col = "blue", main = "Histogram After Filtering",
     xlab = "Layer1", ylab = "Frequency")

# Clip the points to the region ----
# Keep only the points within FW_Sys_design@region@region
points_within_region <- st_intersection(points_sf, FW_Sys_G_design@region@region)


points_transformed <- st_transform(points_within_region, crs = 4326)
coords <- st_coordinates(points_transformed)
points_within_region$longitude <- coords[, 1]
points_within_region$latitude <- coords[, 2]

# Prepare data for distance sampling analysis ====
# Transects are converted to LINESTRING geometries and split into segments using
# the `split_into_segments` function. The segments are then combined into a single
# sf object.
transects_segments <- FW_Sys_G_transects@samplers %>%
  st_cast("LINESTRING") %>%
  st_geometry() %>%
  map(split_into_segments) %>%
  bind_rows() %>%
  st_sf()

# Set CRS for the segments
st_crs(transects_segments) <- CRS

# Add sample labels to the segments
transects_segments$Sample.Label <- row_number(transects_segments)

# Calculate distances between moose points and transect segments
#
# Compute the distance between each moose point and each segment of the transects.
distances <- st_distance(points_within_region, transects_segments)

# Identify the closest segment for each moose point
#
# For each moose point, find the nearest transect segment and the distance to it.
closest_segments <- tibble(
  moose_id = seq_len(nrow(points_within_region)),
  Sample.Label = apply(distances, 1, which.min),
  distance = apply(distances, 1, min)
)


# Merge closest segment information with the moose data
#
# Append the closest segment and distance information to the moose dataset.
points_within_region <- points_within_region %>%
  mutate(
    object = row_number(),
    size = 1,
    Effort = 10,
    distance = closest_segments$distance / 1000
  ) %>%
  left_join(closest_segments, by = c("object" = "moose_id")) %>%
  rename(
    Transect.Label = ID,
    Sample.Label = Sample.Label,
    distance = distance.x
  ) %>%
  select(-distance.y)
# Filter out moose points with distance greater than 600 meters
#
# Remove moose points that are more than 260 meters from the nearest transect segment.
points_within_region <- points_within_region %>%
  filter(distance <= 0.26)

hist(points_within_region$distance, breaks = 20, col = "blue", main = "Histogram Before Filtering",
     xlab = "Distance", ylab = "Frequency")

points_with_filter <- points_within_region

# Apply the probability function to calculate probabilities for each point
# Check for NA values and handle them before applying the function
points_with_filter$distance_probability <- ifelse(
  is.na(points_with_filter$distance),
  NA,
  sapply(points_with_filter$distance, calculate_probability, type = "distance", lambda = 10, k = 1)
)

# Filter points based on the calculated probabilities
points_with_filter <- points_with_filter[runif(nrow(points_with_filter)) < points_with_filter$distance_probability,]

# Drop the temporary probability column
points_with_filter$distance_probability <- NULL

hist(points_with_filter$distance, breaks = 20, col = "blue", main = "Histogram After Filtering",
     xlab = "Distance", ylab = "Frequency")


# Extract x and y coordinates from the moose sf object
coords <- st_coordinates(points_with_filter)

# Prepare data for density surface modelling
#
# Convert moose data to a dataframe for modelling, including x and y coordinates.
segdata <- points_with_filter %>%
  select(
    Effort, longitude, latitude, Transect.Label, Sample.Label,
    layer1, layer2, layer3
  ) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  mutate(
    x = coords[, 1],
    y = coords[, 2]
  )

# Prepare distance data for density surface modelling
#
# Create a dataframe for distance data, including x and y coordinates.
distdata <- points_with_filter %>%
  select(
    object, longitude, latitude, distance, Effort, size,
    layer1, layer2, layer3
  ) %>%
  st_drop_geometry() %>%
  mutate(
    detected = 1,
    x = coords[, 1],
    y = coords[, 2]
  )

# Prepare observation data
#
# Create a dataframe with observation data.
obsdata <- points_with_filter %>%
  select(object, distance, Effort, Sample.Label, size) %>%
  st_drop_geometry() %>%
  mutate(
    detected = 1,
    x = coords[, 1],
    y = coords[, 2]
  ) %>%
  as.data.frame()

# Define required columns for dataframes
segdata_required <- c("longitude", "latitude", "Transect.Label", "Sample.Label", "x", "y", "Effort")
distdata_required <- c("object", "size", "longitude", "latitude", "x", "y", "Effort")
obsdata_required <- c("object", "Sample.Label", "size", "distance", "Effort")

# Check presence of required columns
check_columns_present(segdata, segdata_required)
check_columns_present(distdata, distdata_required)
check_columns_present(obsdata, obsdata_required)


# DSMExtra ====
#' documentation: https://densitymodelling.github.io/dsmextra/articles/dsmextra.html#quantifying-extrapolation-1
# Convert the raster stack to points (x, y, and values)
raster_points <- rasterToPoints(raster_stack)

# Convert the points to a dataframe
raster_df <- as.data.frame(raster_points, stringsAsFactors = FALSE)
# add the area of each cell, use later as offset
# Calculate the area of each raster cell in square kilometers
raster_df$area <- (resolution / 1000) * (resolution / 1000)

# Rename columns for clarity
colnames(raster_df) <- c("x", "y", "layer1", "layer2", "layer3", "area")

# Convert sfc_GEOMETRY to sf object
sf_object <- st_sf(geometry = FW_Sys_G_transects@samplers$geometry)
# Assign the CRS to the sf object
st_crs(sf_object) <- 3400

# Convert sf object to SpatialLines
spatial_lines <- sf::as_Spatial(sf_object)
# sp::proj4string(spatial_lines) <- st_crs(3400)

#'---------------------------------------------
# Create an outline of the study area boundaries
#'---------------------------------------------
study_area <- raster_df[, c("x", "y")]
study_area$value <- 1
study_area <- raster::rasterFromXYZ(study_area, crs = CRSsp)
study_area <- raster::rasterToPolygons(study_area, dissolve = TRUE)
# study_area <- sp::spTransform(study_area, CRSobj = sp::proj4string(transects))
study_area <- smoothr::smooth(study_area, method = "ksmooth", smoothness = 5)

#'---------------------------------------------
# Produce a simple plot
#'---------------------------------------------
plot(study_area, col = "lightskyblue1") # Region boundary
plot(spatial_lines, add = TRUE, col = "skyblue3") # Survey tracks
maps::map("world", fill = TRUE, col = "grey",
          xlim = range(obs$coords.x1),
          ylim = range(obs$coords.x2), add = TRUE)
pts <- obsdata # Sightings
coordinates(pts) <- ~x + y
axis(1); axis(2); box()
points(pts, pch = 16)


#'---------------------------------------------
# Define environmental covariates of interest
#'---------------------------------------------
predictive_covariates <- c("layer2", "layer3")

# Compute extrapolation for the study area
# This function calculates extrapolation metrics for the given samples and prediction grid.
# Inputs:
#   - samples: A dataframe containing sample data (e.g., transect segments with covariates).
#   - covariate.names: A vector of covariate names to be used in the extrapolation.
#   - prediction.grid: A dataframe representing the prediction grid with covariate values.
#   - coordinate.system: The coordinate reference system (CRS) for spatial data.
# Output:
#   - An object containing extrapolation metrics for the prediction grid.
moose_extrapolation <- compute_extrapolation(samples = segdata,
                                             covariate.names = predictive_covariates,
                                             prediction.grid = raster_df,
                                             coordinate.system = CRSsp)
# Removed duplicate or incomplete line to avoid confusion or errors

# Number of cells subject to univariate extrapolation (see below for definition)
raster_df %>%
  dplyr::filter(!dplyr::between(layer1, min(segdata$layer1), max(segdata$layer1)) |
                  !dplyr::between(layer2, min(segdata$layer2), max(segdata$layer2)) |
                  !dplyr::between(layer3, min(segdata$layer3), max(segdata$layer3))) %>%
  nrow()
str(moose_extrapolation, 2)

compare_covariates(extrapolation.type = "both",
                   extrapolation.object = moose_extrapolation,
                   n.covariates = NULL,
                   create.plots = TRUE,
                   display.percent = TRUE)

#'---------------------------------------------
# Calculate Gower's distances and %N
#'---------------------------------------------
moose_nearby <- compute_nearby(samples = segdata,
                               prediction.grid = raster_df,
                               coordinate.system = CRSsp,
                               covariate.names = predictive_covariates,
                               nearby = 1)

#'---------------------------------------------
# Rename coordinates and convert to SpatialPointsdf
#'---------------------------------------------
obs_sp <- obsdata %>%
  # dplyr::rename(., x = coords.x1, y = coords.x2) %>%
  sp::SpatialPointsDataFrame(coords = cbind(.$x, .$y), data = ., proj4string = CRSsp)
# sp::spTransform(., CRSobj = aftt_crs)


survey_extent <- methods::as(raster::extent(range(moose_extrapolation$prediction.grid$x), range(moose_extrapolation$prediction.grid$y)), "SpatialPolygons")

map_extrapolation(map.type = "extrapolation",
                  extrapolation.object = moose_extrapolation,
                  sightings = obs_sp,
                  tracks = spatial_lines)

map_extrapolation(map.type = "mic",
                  extrapolation.object = moose_extrapolation,
                  sightings = obs_sp,
                  tracks = spatial_lines)


map_extrapolation(map.type = "nearby",
                  extrapolation.object = moose_nearby,
                  sightings = obs_sp,
                  tracks = spatial_lines)



# everything in on of DSMextra
moose_analysis <- extrapolation_analysis(samples = segdata,
                                         covariate.names = predictive_covariates,
                                         prediction.grid = raster_df,
                                         coordinate.system = CRSsp,
                                         compare.covariates = TRUE,
                                         compare.extrapolation.type = "both",
                                         compare.n.covariates = NULL,
                                         compare.create.plots = TRUE,
                                         compare.display.percent = TRUE,
                                         nearby.compute = TRUE,
                                         nearby.nearby = 1,
                                         map.generate = TRUE,
                                         map.sightings = obs_sp,
                                         map.tracks = NULL)


## Distance Sampling Analysis ====
#' documentation: https://distancesampling.org/dsm/articles/lines_gomex/mexico-analysis.html
# Model detection function including covariates ----
detfc_hr <- Distance::ds(data = distdata, truncation = max(distdata$distance), transect = "line", key = "hr", formula=~as.factor(layer1))
detfc_hn <- Distance::ds(data = distdata, truncation = max(distdata$distance), transect = "line", key = "hn", formula=~as.factor(layer1))
detfc <- detfc_hn
gof_ds(detfc)
# Define density surface models with varying smoothness parameters (k) and basis types (bs)
# The k parameter controls the complexity of the smooth terms, ensuring sufficient flexibility.

# Model 0: Simple smooth over x and y without specifying k (default k is used)
dsm0 <- dsm::dsm(count ~ s(x, y), detfc, segdata, obsdata, method = "REML")

# Model 1: Simple smooth over x and y without specifying k (default k is used)
dsm1 <- dsm::dsm(abundance.est ~ s(x, y, k = 10), detfc, segdata, obsdata, method = "REML")

# Model 2: Smooth over x, y, and covariates layer2 and layer3 with specified k values
dsm2 <- dsm::dsm(abundance.est ~ s(x, y, k = 10) + s(layer2, k = 20) + s(layer3, k = 20), detfc, segdata, obsdata, method = "REML")

# Model 3: Separate smooths for x, y, and covariates with specified k values
dsm3 <- dsm::dsm(abundance.est ~ s(x, k = 10) +
  s(y, k = 10) +
  s(layer2, k = 20) +
  s(layer3, k = 20), detfc, segdata, obsdata, method = "REML")

# Model 4: Tensor product smooths for x, y, and covariates with basis type "ts"
dsm4 <- dsm::dsm(abundance.est ~ s(x, y, bs = "ts", k = 10) +
  s(layer2, bs = "ts", k = 20) +
  s(layer3, bs = "ts", k = 20), detfc, segdata, obsdata, method = "REML")

# Model 5: Separate tensor product smooths for x, y, and covariates with basis type "ts"
dsm5 <- dsm::dsm(abundance.est ~ s(x, bs = "ts", k = 10) +
  s(y, bs = "ts", k = 10) +
  s(layer2, bs = "ts", k = 20) +
  s(layer3, bs = "ts", k = 20), detfc, segdata, obsdata, method = "REML")


# The k parameter provided to s (and te) terms in dsm controls the complexity of the smooths in the model.
# By setting the k parameter we specify the largest complexity for that smooth term in the model; as long as this is high enough, we can be sure that there is enough flexibility.
# In the output from gam.check above, we can see that there is a “p-value” calculated for the size of the basis, this can be a good guide as to whether the basis size needs to be increased.
# The ?choose.k manual page from mgcv gives further guidance and technical details on this matter.

summary(dsm1)

gam.check(dsm1)
# Deviance explained
round(summary(dsm1)$dev.expl * 100, 2)
round(summary(dsm2)$dev.expl * 100, 2)
round(summary(dsm3)$dev.expl * 100, 2)
round(summary(dsm4)$dev.expl * 100, 2)
round(summary(dsm5)$dev.expl * 100, 2)

# Ensure required columns exist in raster_df
required_columns <- c("layer2", "layer3", "area")
missing_columns <- setdiff(required_columns, colnames(raster_df))
if (length(missing_columns) > 0) {
  stop("Missing required columns in raster_df: ", paste(missing_columns, collapse = ", "))
}

# Predict density surface using the dsm object

pred0 <- predict(dsm0, raster_df, off.set = raster_df$area)
pred1 <- predict(dsm1, raster_df, off.set = raster_df$area)
pred2 <- predict(dsm2, raster_df, off.set = raster_df$area)
pred3 <- predict(dsm3, raster_df, off.set = raster_df$area)
pred4 <- predict(dsm4, raster_df, off.set = raster_df$area)
pred5 <- predict(dsm5, raster_df, off.set = raster_df$area)
raster_df$pred_density0 <- pred0
raster_df$pred_density1 <- pred1
raster_df$pred_density2 <- pred2
raster_df$pred_density3 <- pred3
raster_df$pred_density4 <- pred4
raster_df$pred_density5 <- pred5
# Create a raster stack from the predicted density columns in raster_df
raster_obj <- stack(
  # rasterFromXYZ(raster_df[, c("x", "y", "pred_density0")]),
  rasterFromXYZ(raster_df[, c("x", "y", "pred_density1")]),
  rasterFromXYZ(raster_df[, c("x", "y", "pred_density2")]),
  rasterFromXYZ(raster_df[, c("x", "y", "pred_density3")]),
  rasterFromXYZ(raster_df[, c("x", "y", "pred_density4")]),
  rasterFromXYZ(raster_df[, c("x", "y", "pred_density5")])
)

# Plot the raster layers
# plot(raster_obj[[1]], main = "Predicted Density 0")
plot(raster_obj[[1]], main = "Predicted Density 1")
plot(raster_obj[[2]], main = "Predicted Density 2")
plot(raster_obj[[3]], main = "Predicted Density 3")
plot(raster_obj[[4]], main = "Predicted Density 4")
plot(raster_obj[[5]], main = "Predicted Density 5")
plot(density)
