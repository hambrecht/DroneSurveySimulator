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
  # # Check if linestring is of type sfc_LINESTRING, convert if necessary
  # if (!inherits(linestring, "sfc_LINESTRING")) {
  #   if (inherits(linestring, "sfc") && any(st_geometry_type(linestring) == "LINESTRING")) {
  #     linestring <- st_cast(linestring, "LINESTRING")
  #   } else {
  #     stop("Input must be a LINESTRING geometry")
  #   }
  # }

  # Ensure total_length is a units object
  total_length <- st_length(linestring)
  if (!inherits(total_length, "units")) {
    total_length <- set_units(total_length, "m") # Assuming meters, adjust as necessary
  }

  # Convert total_length to kilometers
  total_length_km <- set_units(total_length, "km")

  # Calculate the number of segments, ensuring the result is compatible with units
  num_segments <- ceiling(as.numeric(total_length_km))

  # Initialize an empty list for segment lengths
  segment_lengths <- vector("list", num_segments)

  # Adjust the last segment length if necessary
  # This part seems to be missing the calculation of segment_lengths before adjusting the last segment
  # Assuming equal division of segments, here's a placeholder calculation
  equal_segment_length <- total_length_km / num_segments
  segment_lengths <- rep(equal_segment_length, num_segments)

  # Ensure the operation is compatible with units for the last segment length adjustment
  last_segment_length <- total_length_km - sum(segment_lengths[1:(num_segments - 1)])
  if (!inherits(last_segment_length, "units")) {
    last_segment_length <- set_units(last_segment_length, "km") # Adjust unit as necessary
  }
  segment_lengths[num_segments] <- last_segment_length

  # Generate points along the linestring at the specified intervals
  points <- sf::st_line_sample(linestring, sample = seq(0, 1, length.out = num_segments + 1))

  # Create segments between consecutive points
  segments <- map2(
    .x = points[-length(points)],
    .y = points[-1],
    .f = ~ st_sfc(st_linestring(x = st_coordinates(c(.x, .y))), crs = st_crs(linestring))
  )
  
  # Combine all segments into a single sf object
  do.call(rbind, segments)
}

# Assuming 'transects' is your sf object with MULTILINESTRING geometries

# Step 1: Explode MULTILINESTRINGs to LINESTRINGs
transects_linestrings <- transects %>%
  st_cast("LINESTRING")

terra::plot(transects_linestrings)
# Step 2: Apply split_into_segments to each LINESTRING
transects_segments <- transects_linestrings %>%
  st_geometry() %>%
  map(split_into_segments) %>%
  do.call(rbind, .) %>%
  st_sf()

# Check the number of features now
nrow(transects_segmentized)
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
