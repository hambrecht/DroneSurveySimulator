# Density surface modelling
# based on: https://distancesampling.org/R/vignettes/mexico-analysis.html
library(dsm)
library(ggplot2)
library(terra)
library(Distance)
library(sf)
library(dplyr)

# plotting options
gg.opts <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank()
)

# make the results reproducible
set.seed(11123)

# load data
data <- sf::st_read("D:\\WMU\\survey_data\\501_moose_locations.shp")
head(data)
plot(data, max.plot = 1)

# rename columns
data$Transect.Label <- data$name
data$Sample.Label <- paste(data$name, data$ident, sep = "_")
data$object <- data$ident
data$size <- 1

# add and convert units to km
data$Effort <- 10
data$distance <- data$distance / 1000

# remove unnessecary colums
data <- data[, -which(names(data) %in% c("distance_2", "aircraft", "date", "y_proj", "x_proj", "name", "layer", "ident"))]

# check data
head(data)

# create seperate df
segdata <- sf::st_drop_geometry(data[, which(names(data) %in% c("Latitude", "Longitude", "Effort", "Transect.Label", "Sample.Label"))])
distdata <- sf::st_drop_geometry(data[, which(names(data) %in% c("object", "Latitude", "Longitude", "distance", "Effort", "size"))])
distdata$detected <- 1
segdata$x <- distdata$x <- sf::st_coordinates(data)[, 1]
segdata$y <- distdata$y <- sf::st_coordinates(data)[, 2]
obsdata <- sf::st_drop_geometry(data[, which(names(data) %in% c("object", "distance", "Effort", "Sample.Label", "size"))])

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
pred_grid <- intersect(wmu_vect, gridpolygon) # Ensure the correct object (wmu) is referenced
plot(pred_grid)

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
