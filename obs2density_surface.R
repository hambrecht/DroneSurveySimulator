# Density surface modelling
# based on: https://distancesampling.org/R/vignettes/mexico-analysis.html
library(dsm)
library(ggplot2)
library(terra)
library(Distance)
library(sf)
library(dplyr)
library(units)
library(purrr)

# plotting options
gg.opts <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank()
)

# make the results reproducible
set.seed(11123)

# load data
moose <- sf::st_read("D:\\WMU\\survey_data\\501_moose_locations.shp")
transects <- sf::st_read("D:\\WMU\\survey_data\\WMU 501 (2018-2019)\\WMU501_transects_2018.gpx", layer = "tracks") %>%
  sf::st_transform(crs = 3400)
head(moose)
head(transects)
transects <- transects %>%
  dplyr::select(1)
terra::plot(moose, max.plot = 1)
terra::plot(transects)
nrow(transects) # = 177



split_into_segments <- function(linestring) {
  # Ensure total_length is a units object
  total_length <- st_length(linestring)
  if (!inherits(total_length, "units")) {
    total_length <- set_units(total_length, "m") # Assuming meters, adjust as necessary
  }

  # Check if total_length is greater than 0
  if (as.numeric(total_length) <= 0) {
    stop("Linestring has a non-positive length")
  }

  # Convert total_length to kilometers
  total_length_km <- set_units(total_length, "km")

  # Calculate the number of segments, ensuring the result is compatible with units
  num_segments <- ceiling(as.numeric(total_length_km))

  # Check if num_segments is greater than 0
  if (num_segments <= 0) {
    stop("Number of segments calculated is not positive")
  }

  # Assuming equal division of segments
  equal_segment_length <- total_length_km / num_segments
  segment_lengths <- rep(equal_segment_length, num_segments)

  # Adjust the last segment length
  last_segment_length <- total_length_km - sum(segment_lengths[1:(num_segments - 1)])
  if (!inherits(last_segment_length, "units")) {
    last_segment_length <- set_units(last_segment_length, "km")
  }
  segment_lengths[num_segments] <- last_segment_length

  # Generate points along the linestring at the specified intervals
  points <- sf::st_line_sample(linestring, sample = seq(0, 1, length.out = num_segments + 1))

  # Convert MULTIPOINT to POINTs
  points <- st_cast(points, "POINT")

  # Check if points generation is successful
  if (length(points) < 2) {
    stop(paste("Failed to generate sufficient points for segments", id, sep = " "))
  }

  # Create segments between consecutive points
  segments <- map2(
    .x = points[-length(points)],
    .y = points[-1],
    .f = ~ {
      segment <- st_sfc(st_linestring(x = st_coordinates(c(.x, .y))), crs = st_crs(linestring))
      st_sf(geometry = segment)
    }
  )

  # Check if segments are created successfully
  if (length(segments) < 1) {
    stop("Failed to create segments")
  }

  # Combine all segments into a single sf object using do.call(rbind, ...)
  segments_sf <- do.call(rbind, segments)

  # Check if the final sf object is valid and not empty
  if (nrow(segments_sf) < 1) {
    stop("Final sf object is empty")
  }

  return(segments_sf)
}

# Assuming 'transects' is your sf object with MULTILINESTRING geometries

# Step 1: Explode MULTILINESTRINGs to LINESTRINGs
transects_linestrings <- transects %>%
  st_cast("LINESTRING")


# Step 2: Apply split_into_segments to each LINESTRING
transects_segments <- transects_linestrings %>%
  st_geometry() %>%
  map(split_into_segments) %>%
  bind_rows() %>% # Use bind_rows to combine all sf objects
  st_sf() # Ensure the result is an sf object



# Check the number of features now
nrow(transects_segments)
# rename columns
moose$Transect.Label <- moose$name
moose$object <- row_number(moose) # add unique ID to each observation based on row number
moose$Sample.Label <- paste(moose$name, moose$object, sep = "_")
moose$size <- 1

# TDDO split transects into segements

# add and convert units to km
moose$Effort <- 10
moose$distance <- moose$distance / 1000

# remove unnessecary colums
moose <- moose[, -which(names(moose) %in% c("distance_2", "aircraft", "date", "y_proj", "x_proj", "name", "layer", "ident"))]

# check moose
head(moose)

# create seperate df
segdata <- as.data.frame(sf::st_drop_geometry(moose[, which(names(moose) %in% c("Latitude", "Longitude", "Effort", "Transect.Label", "Sample.Label"))]))
distdata <- as.data.frame(sf::st_drop_geometry(moose[, which(names(moose) %in% c("object", "Latitude", "Longitude", "distance", "Effort", "size"))]))
distdata$detected <- 1
segdata$x <- distdata$x <- sf::st_coordinates(moose)[, 1]
segdata$y <- distdata$y <- sf::st_coordinates(moose)[, 2]
obsdata <- as.data.frame(sf::st_drop_geometry(moose[, which(names(moose) %in% c("object", "distance", "Effort", "Sample.Label", "size"))]))

head(segdata)
head(distdata)
head(obsdata)

# create a prediction grid
# base on https://examples.distancesampling.org/dsm-point/hare_point_transect_dsm-distill.html
# method from http://rfunctions.blogspot.co.uk/2014/12/how-to-create-grid-and-intersect-it.html

# load WMU outline shapefile
wmu <- sf::st_read("D:\\WMU\\base_data\\WMU\\wmu_501_3400.shp")

# Create an empty SpatRaster
grid <- rast(ext(wmu), resolution = 500) # Use ext instead of extent
crs(grid) <- crs(wmu) # Ensure the correct object (wmu) is referenced
gridpolygon <- as.polygons(grid)
wmu_vect <- vect(wmu)
pred_grid <- terra::intersect(wmu_vect, gridpolygon) # Ensure the correct object (wmu) is referenced
terra::plot(pred_grid)

# given the argument fill (the covariate vector to use as the fill) and a name,
# return a geom_polygon object
# fill must be in the same order as the polygon data
grid_plot_obj <- function(fill, name, sf_obj) {
  # Prepare the data
  data <- data.frame(fill = fill)
  names(data) <- name

  # Combine data with geometry from the sf object
  combined_sf <- sf_obj %>%
    mutate(id = row_number()) %>%
    left_join(data, by = c("id" = "row.names"))

  # Plot using ggplot2 and the combined sf object
  ggplot(combined_sf) +
    geom_sf(aes_string(fill = name)) +
    theme_minimal()
}

# Extract centroid coordinates of each polygon
centroids <- terra::centroids(pred_grid)
centroids_sf <- st_as_sf(centroids, coords = c("x", "y"), crs = 3400, agr = "constant")
# Calculate area of each polygon
areas <- terra::expanse(pred_grid)

# Combine coordinates and areas into a data frame
preddata <- data.frame(x = sf::st_coordinates(centroids_sf)[, 1], y = sf::st_coordinates(centroids_sf)[, 2], area = areas)

# Display the first few rows of the data frame
head(preddata)



# Explorartory data anaylsis
# Distance data
# histograms
hist(distdata$distance, main = "", xlab = "Distance (m)")

# Estimating the detection function
detfc.hr.null <- ds(distdata, max(distdata$distance), key = "hr", adjustment = NULL)
summary(detfc.hr.null)
par(mfrow = c(1, 2))
plot(detfc.hr.null, showpoints = FALSE, pl.den = 0, lwd = 2)
ddf.gof(detfc.hr.null$ddf)

# Fitting a DSM
dsm.xy <- dsm(count ~ s(x, y), detfc.hr.null, segdata, obsdata, method = "REML")
summary(dsm.xy)
vis.gam(dsm.xy, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = 100)
plot(data, max.plot = 1, add = TRUE)
gam.check(dsm.xy)
rqgam_check(dsm.xy)

# Autocorrelation
dsm_cor(dsm.xy, max.lag = 10, Segment.Label = "Sample.Label")

# Abundance estimation
dsm.xy.pred <- predict(dsm.xy, preddata, preddata$area)
p <- ggplot() +
  grid_plot_obj(dsm.xy.pred, "Abundance", pred_grid) +
  coord_equal() +
  gg.opts
p <- p + geom_path(aes(x = x, y = y), data = survey.area)
p <- p + labs(fill = "Abundance")
print(p)
