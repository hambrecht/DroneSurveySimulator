
### Example 1: Creating density and population for study_region and run simulation

# Load the required library
library(dsims)
library(sf)
library(dplyr)
library(here)


#' Create polygons when grid of centroids is provided
#'
#' @param center_points A set of center points for the polygons.
#' @param x_dim Width of each polygon.
#' @param y_dim Height of each polygon.
#' @return An sf object containing the created polygons.
create_sf_polygons <- function(center_points, x_dim, y_dim) {
  # Extract CRS from the center points
  crs <- st_crs(center_points)

  # Function to create a single polygon
  create_polygon <- function(center, x_dim, y_dim) {
    x_center <- st_coordinates(center)[1]
    y_center <- st_coordinates(center)[2]

    # Define the vertices of the rectangle
    vertices <- matrix(c(
      x_center - x_dim / 2, y_center - y_dim / 2,
      x_center + x_dim / 2, y_center - y_dim / 2,
      x_center + x_dim / 2, y_center + y_dim / 2,
      x_center - x_dim / 2, y_center + y_dim / 2,
      x_center - x_dim / 2, y_center - y_dim / 2
    ), ncol = 2, byrow = TRUE)

    # Create the polygon
    st_polygon(list(vertices))
  }

  # Create polygons for each center point
  polygons <- lapply(seq_len(nrow(center_points)), function(i) {
    create_polygon(center_points[i,], x_dim, y_dim)
  })

  # Generate IDs based on coordinates
  generate_id <- function(center) {
    coords <- st_coordinates(center)
    paste0(round(coords[1], 0), " ", round(coords[2], 0))
  }

  ids <- sapply(seq_len(nrow(center_points)), function(i) {
    generate_id(center_points[i,])
  })


  # Combine polygons into an sf object and add IDs
  sf_polygons <- st_sf(
    geometry = st_sfc(polygons),
    crs = crs,
    ID = ids
  )

  sf_polygons
}


# Define the study area
area_m <- matrix(c(0,0,500,0,500,500,0,500,0,0), ncol=2, byrow=TRUE)
area <- sf::st_polygon(list(area_m))

# Plot the study area and the polygons
plot(area)

# Create regions for the study area and subplots
region <- make.region(region.name = "study area", shape = area)


## Helicopter design
heli_design <- make.design(
  region = region,
  transect.type = "line",
  design = "segmentedgrid",
  spacing = 12, # segments seperated by 1.2km
  seg.length = 100, # segements of 10km
  design.angle = 0, # align transect with north south
  seg.threshold = 10, # any segments less than 10% of the segment length (i.e. 1km) will be discarded.
  edge.protocol = "minus",
  truncation = 5, # IMAGE_WIDTH
)
heli_transects <- generate.transects(heli_design)

heli_transects@line.length
plot(region, heli_transects, lwd = 0.5, col = 4)


grid_center <- make.coverage(region,
                             # spacing = 1000
                             n.grid.points = 1
)
plot(grid_center)
# plot(region, grid_center)
# remove coverage.scores column
grid_center@grid <- grid_center@grid %>% select(-coverage.scores)


# Create a buffer around the polygon boundary
buffer_1000m <- st_buffer(region@region, dist = -125)

# Identify points within 1000 meters of the boundary
selection_index <- st_disjoint(grid_center@grid, buffer_1000m, sparse = FALSE)
points_within_1000m <- grid_center@grid[selection_index,]

# Find the nearest points on the polygon for each point in the grid
nearest_lines <- st_nearest_points(points_within_1000m, buffer_1000m)
nearest_points <- st_cast(nearest_lines, "POINT")
new_coords <- st_coordinates(nearest_points[seq(2, length(nearest_points), 2)])

# Replace the coordinates in points1 with the new coordinates
st_geometry(grid_center@grid)[selection_index] <- st_sfc(lapply(1:nrow(points_within_1000m), function(i) st_point(new_coords[i,])))

# Create sub plot polygons
polygons <- create_sf_polygons(grid_center@grid, 100, 100)
plot(area)
plot(polygons, add = TRUE)

poly_region <- make.region(region.name = "sub plots", shape = polygons)
fixW_design <- make.design(
  region = poly_region,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = rep(heli_transects@line.length / length(poly_region@strata.name), length(poly_region@strata.name)), # OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 5, # IMAGE_WIDTH
)
fixW_transects <- generate.transects(fixW_design)
plot(region, fixW_transects, lwd = 0.5, col = 4)


grid_center <- make.coverage(region,
                             # spacing = 1000
                             n.grid.points = 16
)
plot(grid_center)
# plot(region, grid_center)
# remove coverage.scores column
grid_center@grid <- grid_center@grid %>% select(-coverage.scores)


# Create a buffer around the polygon boundary
buffer_1000m <- st_buffer(region@region, dist = -20)

# Identify points within 1000 meters of the boundary
selection_index <- st_disjoint(grid_center@grid, buffer_1000m, sparse = FALSE)
points_within_1000m <- grid_center@grid[selection_index,]

# Find the nearest points on the polygon for each point in the grid
nearest_lines <- st_nearest_points(points_within_1000m, buffer_1000m)
nearest_points <- st_cast(nearest_lines, "POINT")
new_coords <- st_coordinates(nearest_points[seq(2, length(nearest_points), 2)])

# Replace the coordinates in points1 with the new coordinates
st_geometry(grid_center@grid)[selection_index] <- st_sfc(lapply(1:nrow(points_within_1000m), function(i) st_point(new_coords[i,])))

# Create sub plot polygons
polygons <- create_sf_polygons(grid_center@grid, 50, 40)
plot(area)
plot(polygons, add = TRUE)

poly_region <- make.region(region.name = "sub plots", shape = polygons)
quad_design <- make.design(
  region = poly_region,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = rep(heli_transects@line.length / length(poly_region@strata.name), length(poly_region@strata.name)), # OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 5, # IMAGE_WIDTH
)
quad_transects <- generate.transects(quad_design)
plot(region, quad_transects, lwd = 0.5, col = 4)


# Helicopter transects
# Set the output file path using the here package
output_path <- here("Output", "Plots", "exHeliplot.tiff")

# Open a TIFF device
tiff(filename = output_path, width = 2250, height = 2250, units = "px", res = 300)

# Set graphical parameters for A4 size
par(mar = c(5, 5, 4, 2) + 0.1, cex = 0.9, cex.lab = 0.9, cex.axis = 0.9, lwd = 1, bty = "n", family = "Arial")

# Plot the data
plot(region, heli_transects, lwd = 0.5, col = 4)

# Close the TIFF device
dev.off()

# Fixedwing transects
# Set the output file path using the here package
output_path <- here("Output", "Plots", "exFWplot.tiff")

# Open a TIFF device
tiff(filename = output_path, width = 2250, height = 2250, units = "px", res = 300)

# Set graphical parameters for A4 size
par(mar = c(5, 5, 4, 2) + 0.1, cex = 0.9, cex.lab = 0.9, cex.axis = 0.9, lwd = 1, bty = "n", family = "Arial")

# Plot the data
plot(region, fixW_transects, lwd = 0.5, col = 4)

# Close the TIFF device
dev.off()

# Quadcopter transects
# Set the output file path using the here package
output_path <- here("Output", "Plots", "exQCplot.tiff")

# Open a TIFF device
tiff(filename = output_path, width = 2250, height = 2250, units = "px", res = 300)

# Set graphical parameters for A4 size
par(mar = c(5, 5, 4, 2) + 0.1, cex = 0.9, cex.lab = 0.9, cex.axis = 0.9, lwd = 1, bty = "n", family = "Arial")

# Plot the data
plot(region, quad_transects, lwd = 0.5, col = 4)

# Close the TIFF device
dev.off()



