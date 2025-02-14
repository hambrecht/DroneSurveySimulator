# Load necessary libraries
library(here)
library(dsims)
library(knitr)
library(RColorBrewer)
library(dplyr)
library(sf)
library(units)
library(geosphere)

# # Check if pbapply is installed
# if (!requireNamespace("pbapply", quietly = TRUE)) {
#   message("The 'pbapply' package is not installed. Installing it now...")
#   install.packages("pbapply")
# } else {
#   message("The 'pbapply' package is already installed.")
# }

# Define functions

#' Calculate Image Width Based on Altitude
#'
#' This function calculates the width of an image captured from a given altitude.
#'
#' @param ALTITUDE Numeric value representing the altitude (in meters). Must be positive.
#' @param CAMERA_FOV Numeric value representing the field of view of the camera in degrees. Default is 60 degrees.
#' @param CAMERA_ANGLE Numeric value representing the adjustment in degrees. Must be non-negative and not larger than FOV.
#'
#' @return Numeric value representing the rounded image width in meters.
#' @throws Error if ALTITUDE, FOV, or adjustment are not positive or if adjustment is larger than FOV.
#'
#' @examples
#' calculate_image_width(ALTITUDE = 100)
#' calculate_image_width(ALTITUDE = 150, FOV = 75, adjustment = 10)
#'
#' @export
calculate_image_width <- function(ALTITUDE, CAMERA_FOV = 25, CAMERA_ANGLE = 0) {
  if (ALTITUDE <= 0 ||
    CAMERA_FOV <= 0 ||
    CAMERA_ANGLE < 0) {
    stop("Altitude and camera FOV must be positive numbers, and camera angle must be non-negative")
  }

  # Adjust the FOV by doubling the adjustment and adding it to the FOV
  adjusted_FOV <- CAMERA_FOV + 2 * CAMERA_ANGLE

  # Calculate and round the image width
  round(2 * ALTITUDE * tan((adjusted_FOV * pi / 180) / 2), -1)
}

## design stats
## For details see: https://examples.distancesampling.org/dssd-getting-started/GettingStarted-distill.html#appendix-trackline-and-cyclic-trackline-lengths
## Line length = on effort line length
## The trackline length is the sum of the lengths of the transects plus the off-effort transit time required to complete the survey from the beginning of the first transect to the end of the last transect. The off-effort transit distance is calculated as the crow flies and may be longer in reality if transit is required around lakes, islands or coastlines etc.
## The cyclic trackline length is the trackline length plus the off-effort transit distance required to return from the end of the last transect to the beginning of the first transect.

#' Extract key metrics from each simulation summary
#'
#' @param design A design object containing the simulation summary.
#' @return A list of key metrics extracted from the simulation summary.
extract_design_metrics <- function(design) {
  list(
    design_type = design@design,
    mean_sampler_count = design@design.statistics$sampler.count[2, 'Total'],
    mean_cover_area = design@design.statistics$cov.area[2, 'Total'],
    mean_cover_percentage = design@design.statistics$p.cov.area[2, 'Total'],
    mean_line_length = design@design.statistics$line.length[2, 'Total'],
    mean_trackline = design@design.statistics$trackline[2, 'Total'],
    mean_cyclic_trackline = design@design.statistics$cyclictrackline[2, 'Total'],
    mean_on_effort = design@design.statistics$line.length[2, 'Total'],
    mean_off_effort = design@design.statistics$trackline[2, 'Total'] - design@design.statistics$line.length[2, 'Total'],
    mean_return2home = design@design.statistics$cyclictrackline[2, 'Total'] - design@design.statistics$trackline[2, 'Total'],
    mean_off_effort_return = design@design.statistics$cyclictrackline[2, 'Total'] - design@design.statistics$line.length[2, 'Total'],
    on_effort_percentage = round((design@design.statistics$line.length[2, 'Total'] / design@design.statistics$cyclictrackline[2, 'Total']) * 100, 2),
    off_effort_percentage = round(((design@design.statistics$trackline[2, 'Total'] - design@design.statistics$line.length[2, 'Total']) / design@design.statistics$cyclictrackline[2, 'Total']) * 100, 2),
    return2home_percentage = round(((design@design.statistics$cyclictrackline[2, 'Total'] - design@design.statistics$trackline[2, 'Total']) / design@design.statistics$cyclictrackline[2, 'Total']) * 100, 2),
    off_effort_return_percentage = round(((design@design.statistics$cyclictrackline[2, 'Total'] - design@design.statistics$line.length[2, 'Total']) / design@design.statistics$cyclictrackline[2, 'Total']) * 100, 2)
  )
}

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

#' Function to generate a random square polygon with random orientation
#'
#' @param x Width of the polygon.
#' @param y Height of the polygon.
#' @param crs Coordinate reference system.
#' @param area Area within which to generate the polygon.
#' @return A list containing the generated polygon.
generate_random_polygon <- function(x, y, crs, area) {
  center_x <- runif(1, min(st_bbox(area)$xmin), max(st_bbox(area)$xmax))
  center_y <- runif(1, min(st_bbox(area)$ymin), max(st_bbox(area)$ymax))

  coords <- matrix(c(-x / 2, -y / 2, x / 2, -y / 2, x / 2, y / 2, -x / 2, y / 2, -x / 2, -y / 2), ncol = 2, byrow = TRUE)

  coords[, 1] <- coords[, 1] + center_x
  coords[, 2] <- coords[, 2] + center_y

  polygon <- st_polygon(list(coords))
  list(polygon = st_sfc(polygon, crs = crs))
}

#' Function to check if the polygon is within acceptable bounds with a buffer
#'
#' @param polygon The polygon to check.
#' @param area The area within which the polygon should be.
#' @param buffer_distance The buffer distance to apply.
#' @return TRUE if the polygon is within bounds, FALSE otherwise.
check_polygon <- function(polygon, area, buffer_distance) {
  st_crs(polygon) <- st_crs(area)
  buffered_area <- st_buffer(area, -buffer_distance)
  intersection <- st_intersection(polygon, buffered_area)

  if (length(intersection) == 0) {
    return(FALSE)
  }

  if (st_area(intersection) == st_area(polygon)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Function to check if the polygon overlaps with any existing polygons
#'
#' @param polygon The polygon to check.
#' @param existing_polygons A list of existing polygons.
#' @return TRUE if the polygon overlaps with any existing polygons, FALSE otherwise.
check_overlap <- function(polygon, existing_polygons) {
  if (length(existing_polygons) == 0) {
    return(FALSE)
  }

  for (existing_polygon in existing_polygons) {
    if (length(st_intersection(polygon, existing_polygon)) > 0) {
      return(TRUE)
    }
  }

  return(FALSE)
}

#' Function to place N polygons within the area
#'
#' @param N Number of polygons to place.
#' @param x Width of each polygon.
#' @param y Height of each polygon.
#' @param area Area within which to place the polygons.
#' @param buffer_distance Buffer distance to apply.
#' @param max_attempts Maximum number of attempts to place a polygon.
#' @param max_retries Maximum number of retries to place all polygons.
#' @return An sf object containing the placed polygons.
# place_polygons <- function(N, x, y, area, buffer_distance, max_attempts = 1000, max_retries = 1000) {
#   crs <- st_crs(area)
#
#   for (retry in 1:max_retries) {
#     print(paste0("Retry ", retry))
#     polygons <- list()
#     success <- TRUE
#
#     # Place the first polygon in the center of the area
#     center_x <- (st_bbox(area)$xmin + st_bbox(area)$xmax) / 2
#     center_y <- (st_bbox(area)$ymin + st_bbox(area)$ymax) / 2
#     coords <- matrix(c(
#       center_x - x / 2, center_y - y / 2,
#       center_x + x / 2, center_y - y / 2,
#       center_x + x / 2, center_y + y / 2,
#       center_x - x / 2, center_y + y / 2,
#       center_x - x / 2, center_y - y / 2
#     ), ncol = 2, byrow = TRUE)
#     polygon <- st_polygon(list(coords))
#     polygons[[1]] <- st_sfc(polygon, crs = crs)
#
#     for (i in 2:N) {
#       attempts <- 0
#       repeat {
#         attempts <- attempts + 1
#         if (attempts > max_attempts) {
#           success <- FALSE
#           break
#         }
#         result <- generate_random_polygon(x, y, crs, area)
#         polygon <- result$polygon
#         if (check_polygon(polygon, area, buffer_distance) && !check_overlap(polygon, polygons)) {
#           polygons[[i]] <- polygon
#           break
#         }
#       }
#       if (!success) break
#     }
#
#     if (success) {
#       sf_polygons <- st_sf(geometry = do.call(c, polygons), crs = crs)
#       sf_polygons$ID <- LETTERS[1:N]
#       return(sf_polygons)
#     }
#   }
#
#   stop("Failed to place polygons after 100 retries")
# }


# Function to create a grid with cell size half the size of the polygon
create_grid <- function(area, x, y) {
  cell_size_x <- x / 2
  cell_size_y <- y / 2
  grid <- st_make_grid(area, cellsize = c(cell_size_x, cell_size_y), what = "centers")
  grid <- st_as_sf(grid)
  grid <- grid[st_within(grid, area, sparse = FALSE), ]
  return(grid)
}

# Function to place polygons using the grid cells
place_polygons_grid <- function(N, x, y, area, buffer_distance, max_attempts = 1000) {
  crs <- st_crs(area)
  grid <- create_grid(area, x, y)
  polygons <- list()
  success <- TRUE

  for (i in 1:N) {
    attempts <- 0
    repeat {
      attempts <- attempts + 1
      if (attempts > max_attempts) {
        success <- FALSE
        break
      }

      # Select a random grid cell
      center <- grid[sample(nrow(grid), 1), ]
      center_x <- st_coordinates(center)[1]
      center_y <- st_coordinates(center)[2]
      coords <- matrix(c(
        center_x - x / 2, center_y - y / 2,
        center_x + x / 2, center_y - y / 2,
        center_x + x / 2, center_y + y / 2,
        center_x - x / 2, center_y + y / 2,
        center_x - x / 2, center_y - y / 2
      ), ncol = 2, byrow = TRUE)
      polygon <- st_polygon(list(coords))
      polygon_sfc <- st_sfc(polygon, crs = crs)

      if (check_polygon(polygon_sfc, area, buffer_distance) && !check_overlap(polygon_sfc, polygons)) {
        polygons[[i]] <- polygon_sfc
        break
      }
    }
    if (!success) break
  }

  if (success) {
    sf_polygons <- st_sf(geometry = do.call(c, polygons), crs = crs)
    sf_polygons$ID <- LETTERS[1:N]
    return(sf_polygons)
  }

  stop("Failed to place polygons after max attempts")
}


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

#' Function to calculate distance between two points
#'
#' @param p1 First point.
#' @param p2 Second point.
#' @return The distance between the two points.
distance <- function(p1, p2) {
  sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
}

#' Find the longest dimension and its angle
#'
#' @param coords Coordinates of the vertices.
#' @return A list containing the angle and the two points defining the longest dimension.
find_longest_dimension_angle <- function(coords) {
  max_distance <- 0
  point1 <- point2 <- NULL

  # Calculate distances between all pairs of vertices
  for (i in 1:(nrow(coords) - 1)) {
    for (j in (i + 1):nrow(coords)) {
      dist <- distance(coords[i,], coords[j,])
      if (dist > max_distance) {
        max_distance <- dist
        point1 <- coords[i,]
        point2 <- coords[j,]
      }
    }
  }

  # Calculate the angle of the longest dimension
  dx <- point2[1] - point1[1]
  dy <- point2[2] - point1[2]
  angle <- atan2(dy, dx) * (180 / pi)

  # Adjust angle to be in the range [0, 180)
  if (angle < 0) {
    angle <- angle + 180
  }

  return(list(angle = angle, point1 = point1, point2 = point2))
}

#' Extract key metrics from each simulation summary
#'
#' @param sim A simulation object containing the summary.
#' @return A list of key metrics extracted from the simulation summary.
extract_metrics <- function(sim) {
  summary_data <- summary(sim, description.summary = FALSE)

  list(
    mean_estimate = summary_data@individuals$N$mean.Estimate,
    percent_bias = summary_data@individuals$N$percent.bias,
    rmse = summary_data@individuals$N$RMSE,
    ci_coverage_prob = summary_data@individuals$N$CI.coverage.prob,
    mean_se = summary_data@individuals$N$mean.se,
    sd_of_means = summary_data@individuals$N$sd.of.means,
    mean_cover_area = summary_data@individuals$summary$mean.Cover.Area,
    mean_effort = summary_data@individuals$summary$mean.Effort,
    mean_n = summary_data@individuals$summary$mean.n,
    mean_k = summary_data@individuals$summary$mean.k,
    mean_ER = summary_data@individuals$summary$mean.ER,
    mean_se_ER = summary_data@individuals$summary$mean.se.ER
  )
}

# Load density data
wmu_number_list <- c("501", "503", "512", "517", "528")
wmu_number <- wmu_number_list[2]
input_path <- here("Output", "Density", paste0("density", wmu_number, ".RData"))
load(file = input_path)

# load cover
output_path <- here("Output", "Simulation", paste0("cover-WMU", wmu_number, ".RData"))
load(file = output_path)




# Define constants
ALTITUDE <- 120 # Height in meters
CAMERA_HFOV <- 25 # Horizontal FOV in degrees. Max adjustment of 25 degrees. If more than 25 degrees then a third camera or gimbal would be needed to cover 0. Proposed intervals for adjustments are 0, 10, 20, 25
CAMERA_ANGLE <- 35 # Adjustment in degrees; max 25 with two fix cameras or 35 with gimbal in forested
IMAGE_WIDTH <- calculate_image_width(ALTITUDE, CAMERA_HFOV, CAMERA_ANGLE)
print(paste0("Half swath width is: ", IMAGE_WIDTH, " m"))

# Define cover grid spacing and repetition
COV_SPACE <- 500
COV_REPS <- 100

# Calculate the angle of longest Dimension of WMU, only needed for Zigzag
# COORDS <- st_coordinates(region@region)[, 1:2]
# LONGEST_DIMENSION <- find_longest_dimension_angle(COORDS)
# TRANSECT_ANGLE <- LONGEST_DIMENSION$angle

# create coverage grid
cover <- make.coverage(region,
                       spacing = COV_SPACE # OR
                       # n.grid.points = 1000
)
# plot(region, cover)


# Define survey design
## Helicopter design
H_SG_design <- make.design(
  region = region,
  transect.type = "line",
  design = "segmentedgrid",
  spacing = 1200, # segments seperated by 1.2km
  seg.length = 10000, # segements of 10km
  design.angle = 0, # align transect with north south
  seg.threshold = 10, # any segments less than 10% of the segment length (i.e. 1km) will be discarded.
  edge.protocol = "minus",
  truncation = 600, # IMAGE_WIDTH
  coverage.grid = cover
)
H_SG_transects <- generate.transects(H_SG_design)


### Coverage
#### You can re-run the coverage simulation using the following code. Note, your
#### results should vary slightly from mine, make sure you haven't set a seed!
H_SG_design <- run.coverage(H_SG_design, reps = COV_REPS)
total_length <- H_SG_design@design.statistics$line.length[2]
# total_length <- H-SG_transects@line.length


# Drone survey designs
## Fix-wing
# Compute polygon dimensions
number_blocks <- round(total_length / 367200) # 367.2km is the total distance superwake can fly, assuming a speed of 17m/s and a flight time of 6h.
spacing <- 800
poly_dim <- find_best_block_dim(total_length, number_blocks, spacing)
# Check if the polygons can fit into the area


# Checking that blocks fit within region
if (region@area > (poly_dim$x_length * poly_dim$y_length * number_blocks)) {
  # Create a buffer around the polygon boundary
  buffer <- st_buffer(region@region, dist = -1000)

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
  FW_poly <- st_intersection(FW_poly, wmu)
  # Cast the geometries to POLYGON or MULTIPOLYGON
  FW_poly <- st_cast(FW_poly, "MULTIPOLYGON")

  if (number_blocks < nrow(FW_poly)) {
    areas <- st_area(FW_poly)

    # Identify the index of the smallest polygon
    smallest_index <- which.min(areas)

    # Remove the smallest polygon
    FW_poly <- FW_poly[-smallest_index,]
  }
}
plot(st_geometry(wmu))
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
FW_Sys_design <- make.design(
  region = FW_plots,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = numeric(0), # rep(total_length / length(FW_plots@strata.name), length(FW_plots@strata.name)) OR
  spacing = spacing,
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 260, # IMAGE_WIDTH
  coverage.grid = cover
)
FW_Sys_design@truncation <- 260
FW_Sys_transects <- generate.transects(FW_Sys_design)
### Coverage
FW_Sys_design <- run.coverage(FW_Sys_design, reps = COV_REPS)


# Fixed wing zigzag flights
FW_ZZ_design <- make.design(
  region = FW_plots,
  transect.type = "line",
  design = "eszigzag",
  samplers = numeric(0), # OR
  line.length = FW_Sys_design@design.statistics$line.length[2,], # rep(total_length / length(FW_plots@strata.name), length(FW_plots@strata.name)) OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 260, # IMAGE_WIDTH
  coverage.grid = cover
)
FW_ZZ_design@truncation <- 260
FW_ZZ_transects <- generate.transects(FW_ZZ_design)
### Coverage
FW_ZZ_design <- run.coverage(FW_ZZ_design, reps = COV_REPS)


## Quadcopter
## Systematic subplot design
# number_blocks <- round(total_length / 26000)
# spacing <- 200

# # Checking that blocks fit within region
# if (region@area > (2000 * 2500 * number_blocks)) {
#   # create coverage grid
#   grid_center <- make.coverage(region,
#                                # spacing = 1000
#                                n.grid.points = number_blocks
#   )
#   # plot(region, grid_center)
#   # remove coverage.scores column
#   grid_center@grid <- grid_center@grid %>% select(-coverage.scores)


#   # Create a buffer around the polygon boundary
#   buffer_1000m <- st_buffer(region@region, dist = -1001)

#   # Identify points within 1000 meters of the boundary
#   selection_index <- st_disjoint(grid_center@grid, buffer_1000m, sparse = FALSE)
#   points_within_1000m <- grid_center@grid[selection_index,]

#   # Find the nearest points on the polygon for each point in the grid
#   nearest_lines <- st_nearest_points(points_within_1000m, buffer_1000m)
#   nearest_points <- st_cast(nearest_lines, "POINT")
#   new_coords <- st_coordinates(nearest_points[seq(2, length(nearest_points), 2)])

#   # Replace the coordinates in points1 with the new coordinates
#   st_geometry(grid_center@grid)[selection_index] <- st_sfc(lapply(seq_len(nrow(points_within_1000m)), function(i) st_point(new_coords[i,])))

#   # Create sub plot polygons
#   polygons <- create_sf_polygons(grid_center@grid, 2500, 2000)
# }
# plot(st_geometry(wmu))
# plot(polygons[1], add = TRUE, col = "red")

# number_blocks
# nrow(polygons)

# # create subplot region
# QC_plots <- make.region(
#   region.name = "study area",
#   shape = polygons,
#   strata.name = polygons$ID
# )

# # plot(QC_plots)
# # create flight lines withing quadcopterplots
# QC_Sys_nadir_design <- make.design(
#   region = QC_plots,
#   transect.type = "line",
#   design = "systematic",
#   samplers = numeric(0), # OR
#   line.length = numeric(0), # rep(total_length / length(QC_plots@strata.name), length(QC_plots@strata.name)) OR
#   spacing = spacing, # numeric(0)
#   design.angle = 0,
#   edge.protocol = "minus",
#   truncation = 50, # IMAGE_WIDTH
#   coverage.grid = cover
# )
# QC_Sys_nadir_transects <- generate.transects(QC_Sys_nadir_design)
# ### Coverage
# system.time(QC_Sys_nadir_design <- run.coverage(QC_Sys_nadir_design, reps = COV_REPS))[3] / 60


## Systematic subplot design
number_blocks <- round(total_length / 26000)
spacing <- 400

# Checking that blocks fit within region
if (region@area > (3200 * 3200 * number_blocks)) {
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
  new_coords <- st_coordinates(nearest_points[seq(2, length(nearest_points), 2)])

  # Replace the coordinates in points1 with the new coordinates
  st_geometry(grid_center@grid)[selection_index] <- st_sfc(lapply(seq_len(nrow(points_within_1000m)), function(i) st_point(new_coords[i,])))

  # Create sub plot polygons
  polygons <- create_sf_polygons(grid_center@grid, 3200, 3200)
  polygons <- st_intersection(polygons, wmu)
 
  # Convert MULTIPOLYGON to POLYGON
  polygons <- st_cast(polygons, "POLYGON")
}
plot(st_geometry(wmu))
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
  truncation = 260, # IMAGE_WIDTH
  coverage.grid = cover
)
QC_Sys_design@truncation <- 260

QC_Sys_transects <- generate.transects(QC_Sys_design)
plot(QC_Sys_transects)
### Coverage
system.time(QC_Sys_design <- run.coverage(QC_Sys_design, reps = COV_REPS))[3] / 60

# hist(get.coverage(H_SG_design))
# Plot desings
par(mfrow = c(2, 4))
plot(region, H_SG_transects, lwd = 0.5, col = 4)
plot(region, FW_Sys_transects, lwd = 0.5, col = 4)
plot(region, FW_ZZ_transects, lwd = 0.5, col = 4)
# plot(region, QC_Sys_nadir_transects, lwd = 0.5, col = 4)
plot(region, QC_Sys_transects, lwd = 0.5, col = 4)
par(mfrow = c(1, 1))
par(mfrow = c(2, 4))
plot(H_SG_design)
plot(FW_Sys_design)
plot(FW_ZZ_design)
# plot(QC_Sys_nadir_design)
plot(QC_Sys_design)
par(mfrow = c(1, 1))

## design stats
## For details see: https://examples.distancesampling.org/dssd-getting-started/GettingStarted-distill.html#appendix-trackline-and-cyclic-trackline-lengths
## Line lenght = on effort line length
## The trackline length is the sum of the lengths of the transects plus the off-effort transit time required to complete the survey from the beginning of the first transect to the end of the last transect. The off-effort transit distance is calculated as the crow flies and may be longer in reality if transit is required around lakes, islands or coastlines etc.
## The cyclic trackline length is the trackline length plus the off-effort transit distance required to return from the end of the last transect to the beginning of the first transect.
# Extract key metrics from each simulation summary
H_SG_design_metric <- extract_design_metrics(H_SG_design)
FW_Sys_design_metric <- extract_design_metrics(FW_Sys_design)
FW_ZZ_design_metric <- extract_design_metrics(FW_ZZ_design)
# QC_Sys_nadir_design_metric <- extract_design_metrics(QC_Sys_nadir_design)
QC_Sys_design_metric <- extract_design_metrics(QC_Sys_design)

# Combine metrics into a single dataframe
design_comparison_df <- data.frame(
  Simulation = c("H-SG", "FW-Sys", "FW-ZZ",  "QC_Sys"),
  Design = c(
    H_SG_design_metric$design_type,
    FW_Sys_design_metric$design_type[1],
    FW_ZZ_design_metric$design_type[1],
    # QC_Sys_nadir_design_metric$design_type[1],
    QC_Sys_design_metric$design_type[1]
  ),
  Mean_Sampler_Count = c(
    H_SG_design_metric$mean_sampler_count,
    FW_Sys_design_metric$mean_sampler_count,
    FW_ZZ_design_metric$mean_sampler_count,
    # QC_Sys_nadir_design_metric$mean_sampler_count,
    QC_Sys_design_metric$mean_sampler_count
  ),
  Mean_Cover_Area = c(
    H_SG_design_metric$mean_cover_area,
    FW_Sys_design_metric$mean_cover_area,
    FW_ZZ_design_metric$mean_cover_area,
    # QC_Sys_nadir_design_metric$mean_cover_area,
    QC_Sys_design_metric$mean_cover_area
  ),
  Mean_Cover_Percentage = c(
    H_SG_design_metric$mean_cover_percentage,
    FW_Sys_design_metric$mean_cover_percentage,
    FW_ZZ_design_metric$mean_cover_percentage,
    # QC_Sys_nadir_design_metric$mean_cover_percentage,
    QC_Sys_design_metric$mean_cover_percentage
  ),
  Mean_Line_Length = c(
    H_SG_design_metric$mean_line_length,
    FW_Sys_design_metric$mean_line_length,
    FW_ZZ_design_metric$mean_line_length,
    # QC_Sys_nadir_design_metric$mean_line_length,
    QC_Sys_design_metric$mean_line_length
  ),
  Mean_Trackline_Length = c(
    H_SG_design_metric$mean_trackline,
    FW_Sys_design_metric$mean_trackline,
    FW_ZZ_design_metric$mean_trackline,
    # QC_Sys_nadir_design_metric$mean_trackline,
    QC_Sys_design_metric$mean_trackline
  ),
  Mean_Cyclic_Trackline_Length = c(
    H_SG_design_metric$mean_cyclic_trackline,
    FW_Sys_design_metric$mean_cyclic_trackline,
    FW_ZZ_design_metric$mean_cyclic_trackline,
    # QC_Sys_nadir_design_metric$mean_cyclic_trackline,
    QC_Sys_design_metric$mean_cyclic_trackline
  ),
  Mean_On_Effort = c(
    H_SG_design_metric$mean_on_effort,
    FW_Sys_design_metric$mean_on_effort,
    FW_ZZ_design_metric$mean_on_effort,
    # QC_Sys_nadir_design_metric$mean_on_effort,
    QC_Sys_design_metric$mean_on_effort
  ),
  Mean_Off_Effort = c(
    H_SG_design_metric$mean_off_effort,
    FW_Sys_design_metric$mean_off_effort,
    FW_ZZ_design_metric$mean_off_effort,
    # QC_Sys_nadir_design_metric$mean_off_effort,
    QC_Sys_design_metric$mean_off_effort
  ),
  Mean_Return_to_Home = c(
    H_SG_design_metric$mean_return2home,
    FW_Sys_design_metric$mean_return2home,
    FW_ZZ_design_metric$mean_return2home,
    # QC_Sys_nadir_design_metric$mean_return2home,
    QC_Sys_design_metric$mean_return2home
  ),
  Mean_Off_Effort_Return = c(
    H_SG_design_metric$mean_off_effort_return,
    FW_Sys_design_metric$mean_off_effort_return,
    FW_ZZ_design_metric$mean_off_effort_return,
    # QC_Sys_nadir_design_metric$mean_off_effort_return,
    QC_Sys_design_metric$mean_off_effort_return
  ),
  On_Effort_Percentage = c(
    H_SG_design_metric$on_effort_percentage,
    FW_Sys_design_metric$on_effort_percentage,
    FW_ZZ_design_metric$on_effort_percentage,
    # QC_Sys_nadir_design_metric$on_effort_percentage,
    QC_Sys_design_metric$on_effort_percentage
  ),
  Off_Effort_Percentage = c(
    H_SG_design_metric$off_effort_percentage,
    FW_Sys_design_metric$off_effort_percentage,
    FW_ZZ_design_metric$off_effort_percentage,
    # QC_Sys_nadir_design_metric$off_effort_percentage,
    QC_Sys_design_metric$off_effort_percentage
  ),
  Return_to_Home_Percentage = c(
    H_SG_design_metric$return2home_percentage,
    FW_Sys_design_metric$return2home_percentage,
    FW_ZZ_design_metric$return2home_percentage,
    # QC_Sys_nadir_design_metric$return2home_percentage,
    QC_Sys_design_metric$return2home_percentage
  ),
  Off_Effort_Return_Percentage = c(
    H_SG_design_metric$off_effort_return_percentage,
    FW_Sys_design_metric$off_effort_return_percentage,
    FW_ZZ_design_metric$off_effort_return_percentage,
    # QC_Sys_nadir_design_metric$off_effort_return_percentage,
    QC_Sys_design_metric$off_effort_return_percentage
  ),
  Number_of_Plots = c(
    length(H_SG_design_metric$design_type),
    length(FW_Sys_design_metric$design_type),
    length(FW_ZZ_design_metric$design_type),
    # length(QC_Sys_nadir_design_metric$design_type),
    length(QC_Sys_design_metric$design_type)
  )
)
# Print the comparison dataframe
# print(comparison_df)
kable(design_comparison_df)
# drop all but the first 8 rows


# Save simulation data
output_path <- here("Output", "Simulation", paste0("cover-WMU", wmu_number, ".RData"))
# output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number,"-T",IMAGE_WIDTH,"H-SG-DF", detectF@key.function, ".RData"))
save(cover, H_SG_design, FW_Sys_design, FW_ZZ_design, QC_Sys_design, H_SG_transects, FW_Sys_transects, FW_ZZ_transects,  QC_Sys_transects, design_comparison_df, file = output_path)

