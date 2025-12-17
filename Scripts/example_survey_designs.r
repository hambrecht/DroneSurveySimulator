
### Example 1: Creating density and population for study_region and run simulation

# Load the required library
library(dsims)
library(sf)
library(dplyr)
library(here)
# httpgd::hgd()

#' Function to find optimal polygon dimensions
#'
#' @param total_length Total length of the survey area.
#' @param number_blocks Number of blocks to divide the area into.
#' @param spacing Spacing between the blocks.
#' @return A tibble with the best block dimensions.
find_best_block_dim <- function(total_length, number_blocks, spacing) {
  block_length <- total_length / number_blocks
  possible_lines <- seq(1, block_length / spacing, 1)
  
  dimensions <- tibble(
    nr_lines = possible_lines,
    y_length = block_length / possible_lines,
    x_length = spacing * possible_lines,
    difference = abs((block_length / possible_lines) - (spacing * possible_lines))
  )
  
  best_block_dim <- dimensions %>%
    filter(total_length >= (y_length * nr_lines * number_blocks)) %>%
    arrange(difference) %>%
    slice(1)
  
  if (nrow(best_block_dim) == 0) {
    stop("No suitable dimensions found")
  } else {
    return(best_block_dim)
  }
}

find_best_block_dim(1170000, 1, 600)

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
area_m <- matrix(c(0,0,50000,0,50000,50000,0,50000,0,0), ncol=2, byrow=TRUE)
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
  spacing = 1200, # segments seperated by 1.2km
  seg.length = 10000, # segements of 10km
  design.angle = 0, # align transect with north south
  seg.threshold = 10, # any segments less than 10% of the segment length (i.e. 1km) will be discarded.
  edge.protocol = "minus",
  truncation = 500, # IMAGE_WIDTH
)
heli_transects <- generate.transects(heli_design)
output_path <- here("Output", "GPS", "example_Heli_transects.shp")
st_write(heli_transects@samplers, output_path, driver = "ESRI Shapefile")

total_length <- heli_transects@line.length
number_blocks <- 8
spacing <- 310
poly_dim <- NA
poly_dim$x_length <- 10680
poly_dim$y_length <- 10850
# Check if the polygons can fit into the area


# Checking that blocks fit within region
# if (region@area > (poly_dim$x_length * poly_dim$y_length * number_blocks)) {
  # Create a buffer around the polygon boundary
  buffer <- st_buffer(region@region, dist = -2000)
  
  buffer_region <- make.region(region.name = 'buffer', shape = buffer)
  
  # create coverage grid
  grid_center <- make.coverage(buffer_region,
                               # spacing = 1000
                               n.grid.points = number_blocks
  )
  
  # Check if the number of grid points is less than number_blocks
  if (nrow(grid_center@grid) < number_blocks) {
    grid_center <- make.coverage(buffer_region,
                                 # spacing = 1000
                                 n.grid.points = number_blocks + 1
    )
  }
  
  # plot(region, grid_center)
  # remove coverage.scores column
  grid_center@grid <- grid_center@grid %>% select(-coverage.scores)
  
  
  
  #
  # # Identify points within 1000 meters of the boundary
  # selection_index <- st_disjoint(grid_center@grid, buffer_1000m, sparse = FALSE)
  # points_within_1000m <- grid_center@grid[selection_index,]
  #
  # # Find the nearest points on the polygon for each point in the grid
  # nearest_lines <- st_nearest_points(points_within_1000m, buffer_1000m)
  # nearest_points <- st_cast(nearest_lines, "POINT")
  # new_coords <- st_coordinates(nearest_points[seq(2, length(nearest_points), 2)])
  #
  # # Replace the coordinates in points1 with the new coordinates
  # st_geometry(grid_center@grid)[selection_index] <- st_sfc(lapply(1:nrow(points_within_1000m), function(i) st_point(new_coords[i,])))
  
  # Create sub plot polygons
  FW_poly <- create_sf_polygons(grid_center@grid, poly_dim$x_length, poly_dim$y_length)
  FW_poly <- st_intersection(FW_poly, area)
  # Cast the geometries to POLYGON or MULTIPOLYGON
  FW_poly <- st_cast(FW_poly, "MULTIPOLYGON")
  
  # if (number_blocks < nrow(FW_poly)) {
  #   areas <- st_area(FW_poly)
  #   
  #   # Identify the index of the smallest polygon
  #   smallest_index <- which.min(areas)
  #   
  #   # Remove the smallest polygon
  #   FW_poly <- FW_poly[-smallest_index,]
  }
# }
plot(st_geometry(area))
plot(FW_poly[1], add = TRUE, col = "red")
print(paste0("Number of  required plots: ", number_blocks, " Number calculated plots: ",nrow(grid_center@grid)))



# create subplot region
FW_plots <- make.region(
  region.name = "study area",
  shape = FW_poly,
  strata.name = FW_poly$ID
)

# plot(FW_plots)
# create systematic flight lines withing fixed wing
FW_Sys_G_design <- make.design(
  region = FW_plots,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = numeric(0), # rep(total_length / length(FW_plots@strata.name), length(FW_plots@strata.name)) OR
  spacing = spacing,
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 155 # IMAGE_WIDTH
)
fixW_transects <- generate.transects(FW_Sys_G_design)
plot(region, fixW_transects, lwd = 0.5, col = 4)
output_path <- here("Output", "GPS", "example_FW_transects.shp")
st_write(fixW_transects@samplers, output_path, driver = "ESRI Shapefile")

number_blocks <- round(total_length / 26000)

poly_dim$x_length <- 2890
poly_dim$y_length <- 2790

# Checking that blocks fit within region
if (region@area > (poly_dim$x_length * poly_dim$y_length * number_blocks)) {
  # create coverage grid
  grid_center <- make.coverage(region,
                               # spacing = 1000
                               n.grid.points = number_blocks
  )
  # plot(region, grid_center)
  # remove coverage.scores column
  grid_center@grid <- grid_center@grid %>% select(-coverage.scores)
  
  
  # Create a buffer around the polygon boundary
  buffer_1000m <- st_buffer(region@region, dist = -1001)
  
  # Identify points within 1000 meters of the boundary
  selection_index <- st_disjoint(grid_center@grid, buffer_1000m, sparse = FALSE)
  points_within_1000m <- grid_center@grid[selection_index,]
  
  # Find the nearest points on the polygon for each point in the grid
  nearest_lines <- st_nearest_points(points_within_1000m, buffer_1000m)
  nearest_points <- st_cast(nearest_lines, "POINT")
  # new_coords <- st_coordinates(nearest_points[seq(2, length(nearest_points), 2)])
  # 
  # # Replace the coordinates in points1 with the new coordinates
  # st_geometry(grid_center@grid)[selection_index] <- st_sfc(lapply(seq_len(nrow(points_within_1000m)), function(i) st_point(new_coords[i,])))
  
  # Create sub plot polygons
  polygons <- create_sf_polygons(grid_center@grid, poly_dim$x_length, poly_dim$y_length)
  polygons <- st_intersection(polygons, area)
  
  # Convert MULTIPOLYGON to POLYGON
  polygons <- st_cast(polygons, "POLYGON")
}
plot(st_geometry(area))
plot(polygons[1], add = TRUE, col = "red")

number_blocks
nrow(polygons)



# create subplot region
QC_plots <- make.region(
  region.name = "study area",
  shape = polygons,
  strata.name = polygons$ID
)

# plot(QC_plots)
# create flight lines withing quadcopterplots
QC_Sys_design <- make.design(
  region = QC_plots,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = numeric(0), # rep(total_length / length(QC_plots@strata.name), length(QC_plots@strata.name)) OR
  spacing = spacing, # numeric(0)
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 40 # IMAGE_WIDTH
)
quad_transects <- generate.transects(QC_Sys_design)
plot(region, quad_transects, lwd = 0.5, col = 4)

output_path <- here("Output", "GPS", "example_QC_transects.shp")
st_write(quad_transects@samplers, output_path, driver = "ESRI Shapefile", overwrite=TRUE)


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





