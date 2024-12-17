### Example 1: Creating density and population for study_region and run simulation

# Load the required library
library(dsims)
library(sf)
library(dplyr)

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


## Systematic design
sys_design <- make.design(
  region = region,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = heli_transects@line.length, # OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 5, # IMAGE_WIDTH
)
sys_transects <- generate.transects(sys_design)



## Random design
rnd_design <- make.design(
  region = region,
  transect.type = "line",
  design = "random",
  samplers = numeric(0), # OR
  line.length = heli_transects@line.length, # OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 5, # IMAGE_WIDTH
)
rnd_transects <- generate.transects(rnd_design)



## Zigzag design
zigzag_design <- make.design(
  region = region,
  transect.type = "line",
  design = "eszigzag",
  samplers = numeric(0), # OR
  line.length = heli_transects@line.length, # OR
  spacing = numeric(0),
  design.angle = 90, # The design angle for the zigzag designs refers to the angle of a line which would run through the middle of each zigzag transect if the zigzags were to be generated within a rectangle. The design angle for zigzags should usually run along the longest dimension of the study region.
  edge.protocol = "minus",
  bounding.shape = "convex.hull", # rectangle or convex.hull. convex hull is generally more efficient.
  truncation = 5, # IMAGE_WIDTH
)
zigzag_transects <- generate.transects(zigzag_design)




## Zigzag with complementary line
zigzagcom_design <- make.design(
  region = region,
  transect.type = "line",
  design = "eszigzagcom", # eszigzag or eszigzagcom
  samplers = numeric(0), # OR
  line.length = heli_transects@line.length, # OR
  spacing = numeric(0),
  design.angle = 90,
  edge.protocol = "minus",
  bounding.shape = "convex.hull",
  truncation = 5, # IMAGE_WIDTH
)
zigzagcom_transects <- generate.transects(zigzagcom_design)

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



# Plot desings
par(mfrow = c(3, 3))
plot(region, heli_transects, lwd = 1, col = 4)
plot(region, drone_transects, lwd = 1, col = 4)
plot(region, sys_transects, lwd = 1, col = 4)
plot(region, rnd_transects, lwd = 1, col = 4)
plot(region, zigzag_transects, lwd = 1, col = 4)
plot(region, zigzagcom_transects, lwd = 1, col = 2)
par(mfrow = c(1, 1))
par(mfrow = c(3, 3))
plot(heli_design)
plot(drone_design)
plot(sys_design)
plot(rnd_design)
plot(zigzag_design)
plot(zigzagcom_design)
par(mfrow = c(1, 1))


# Plot zigzagcom_transects with alternating colors
plot(area, asp = 0.8, col = rgb(230, 230, 250, maxColorValue = 255))
colors <- c(4, 3) # Define the colors to alternate between
for (i in seq(2, length(zigzagcom_transects@samplers$geometry), by = 2)) {
  color <- colors[((i - 1) %/% 2) %% 2 + 1]
  plot(zigzagcom_transects@samplers$geometry[i], add = TRUE, col = color, lwd = 1)
  plot(zigzagcom_transects@samplers$geometry[i + 1], add = TRUE, col = color, lwd = 1)
}