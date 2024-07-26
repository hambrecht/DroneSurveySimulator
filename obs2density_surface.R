# Density surface modelling
# based on: https://distancesampling.org/R/vignettes/mexico-analysis.html

# Load necessary libraries for analysis and visualization
library(dsm)
library(ggplot2)
library(terra)
library(Distance)
library(sf)
library(dplyr)
library(units)
library(purrr)
# install.packages("C:/Users/lhambrec/Sync/1_Projects/Simulator/dsims_1.0.4.tar.gz", repos = NULL, type = "source")
library(dsims)

# Define constants
GRID_SIZE <- 500

# Function to load and transform spatial data
load_spatial_data <- function(moose_path, transects_path, crs) {
  moose <- sf::st_read(moose_path) %>%
    sf::st_transform(crs = crs) %>%
    select(Latitude, Longitude, date, name)

  transects <- sf::st_read(transects_path, layer = "tracks") %>%
    sf::st_transform(crs = crs) %>%
    select(1)

  list(moose = moose, transects = transects)
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
  transects_path = "D:\\WMU\\survey_data\\WMU 501 (2018-2019)\\WMU501_transects_2018.gpx"
)
data <- load_spatial_data(data_paths$moose_path, data_paths$transects_path, crs = 3400)
moose <- data$moose
transects <- data$transects

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
segdata <- as.data.frame(sf::st_drop_geometry(moose[, c("Latitude", "Longitude", "Effort", "Transect.Label", "Sample.Label")]))
distdata <- as.data.frame(sf::st_drop_geometry(moose[, c("object", "Latitude", "Longitude", "distance", "Effort", "size")]))
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
detfc.hr.null <- ds(distdata, max(distdata$distance), key = "hr", adjustment = NULL)
summary(detfc.hr.null)
par(mfrow = c(1, 2))
plot(detfc.hr.null, showpoints = FALSE, pl.den = 0, lwd = 2)
ddf.gof(detfc.hr.null$ddf)
par(mfrow = c(1, 1))

# Fit a DSM
dsm.xy <- dsm(count ~ s(x, y), detfc.hr.null, segdata, obsdata, method = "REML")
summary(dsm.xy)
vis.gam(dsm.xy, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = GRID_SIZE)
gam.check(dsm.xy)
rqgam_check(dsm.xy)
vis_concurvity(dsm.xy)

# Check for autocorrelation
dsm_cor(dsm.xy, max.lag = 10, Segment.Label = "Sample.Label")

# Estimate abundance
dsm.xy.pred <- predict(dsm.xy, coords, dsm.xy$offset[1])
# Calculate total abundance over the survey area
sum(dsm.xy.pred)

# fill density object with predeicted values
density@density.surface[[1]]$density <- dsm.xy.pred

# plot density
plot(density)
