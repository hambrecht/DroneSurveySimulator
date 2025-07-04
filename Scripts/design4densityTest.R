#' Drone Survey Simulator: Design and Coverage Analysis
#'
#' This script is part of the Drone Survey Simulator project and is used to
#' design and analyze survey coverage for wildlife management units (WMUs).
#' It includes functions for calculating image dimensions, extracting design
#' metrics, generating polygons, and creating systematic survey designs.
#' The script also performs coverage simulations and compares different survey
#' designs based on key metrics.
#'
#' ## Key Features:
#' - **Image Width Calculation**: Calculates the width of an image captured
#'   from a given altitude and camera field of view.
#' - **Design Metrics Extraction**: Extracts key metrics such as mean sampler
#'   count, coverage area, and effort percentages from survey designs.
#' - **Polygon Generation**: Generates random or grid-based polygons within
#'   a specified area, ensuring no overlap and adherence to buffer constraints.
#' - **Systematic Survey Design**: Creates systematic line transect designs
#'   for drone surveys, including baseline and quadcopter designs.
#' - **Coverage Simulation**: Simulates survey coverage and evaluates
#'   performance metrics for different designs.
#' - **Design Comparison**: Compares multiple survey designs based on metrics
#'   such as on-effort and off-effort percentages, trackline lengths, and
#'   coverage area.
#'
#' ## Workflow:
#' 1. **Load Libraries and Data**: Loads necessary R libraries and input data
#'    for density and coverage.
#' 2. **Define Constants**: Sets parameters such as altitude, camera field of
#'    view, and grid spacing.
#' 3. **Calculate Image Width**: Computes the image width based on altitude
#'    and camera specifications.
#' 4. **Generate Coverage Grid**: Creates a grid for coverage analysis.
#' 5. **Design Survey**: Defines systematic survey designs with varying
#'    numbers of samplers and transects.
#' 6. **Run Coverage Simulation**: Simulates coverage for each design and
#'    extracts key metrics.
#' 7. **Compare Designs**: Combines metrics into a dataframe for comparison
#'    and visualization.
#' 8. **Save Results**: Saves simulation data and comparison results to files.
#'
#' ## Outputs:
#' - **Coverage Simulation Data**: Saved as RData files for further analysis.
#' - **Comparison Table**: A CSV file summarizing design metrics for different
#'   survey designs.
#'
#' ## Notes:
#' - The script uses the `dsims` package for survey design and coverage
#'   simulation.
#' - Ensure that the input data files for density and coverage are available
#'   in the specified paths.
#' - The script includes functions for both random and grid-based polygon
#'   placement, with checks for overlap and boundary constraints.
#' - Coverage simulations are run with a specified number of repetitions
#'   (`COV_REPS`) to ensure robust results.
#'
#' ## References:
#' - Distance Sampling: https://examples.distancesampling.org/dssd-getting-started/
#' - `dsims` Package Documentation: https://cran.r-project.org/web/packages/dsims/
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
  if (CAMERA_ANGLE > CAMERA_FOV) {
    stop("Camera angle must not exceed the camera field of view")
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
    mean_sampler_count = design@design.statistics$sampler.count[2, "Total"],
    mean_cover_area = design@design.statistics$cov.area[2, "Total"],
    # mean_cover_percentage = design@design.statistics$p.cov.area[2, 'Total'],
    mean_cover_percentage = round(design@design.statistics$cov.area[2, "Total"] / region@area * 100, 2),
    mean_line_length = design@design.statistics$line.length[2, "Total"],
    mean_trackline = design@design.statistics$trackline[2, "Total"],
    mean_cyclic_trackline = design@design.statistics$cyclictrackline[2, "Total"],
    mean_on_effort = design@design.statistics$line.length[2, "Total"],
    mean_off_effort = design@design.statistics$trackline[2, "Total"] - design@design.statistics$line.length[2, "Total"],
    mean_return2home = design@design.statistics$cyclictrackline[2, "Total"] - design@design.statistics$trackline[2, "Total"],
    mean_off_effort_return = design@design.statistics$cyclictrackline[2, "Total"] - design@design.statistics$line.length[2, "Total"],
    on_effort_percentage = round((design@design.statistics$line.length[2, "Total"] / design@design.statistics$cyclictrackline[2, "Total"]) * 100, 2),
    off_effort_percentage = round(((design@design.statistics$trackline[2, "Total"] - design@design.statistics$line.length[2, "Total"]) / design@design.statistics$cyclictrackline[2, "Total"]) * 100, 2),
    return2home_percentage = round(((design@design.statistics$cyclictrackline[2, "Total"] - design@design.statistics$trackline[2, "Total"]) / design@design.statistics$cyclictrackline[2, "Total"]) * 100, 2),
    off_effort_return_percentage = round(((design@design.statistics$cyclictrackline[2, "Total"] - design@design.statistics$line.length[2, "Total"]) / design@design.statistics$cyclictrackline[2, "Total"]) * 100, 2)
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
      if (nrow(grid) > 0) {
        center <- grid[sample(nrow(grid), 1), ]
      } else {
        stop("The grid is empty. Cannot sample a random grid cell.")
      }
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
  # Generate unique IDs based on coordinates and row index
  generate_id <- function(center, index) {
    coords <- st_coordinates(center)
    paste0("ID_", index, "_", format(coords[1], scientific = FALSE), "_", format(coords[2], scientific = FALSE))
  }

  ids <- sapply(seq_len(nrow(center_points)), function(i) {
    generate_id(center_points[i, ], i)
  })

  ids <- sapply(seq_len(nrow(center_points)), function(i) {
    generate_id(center_points[i, ])
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
      dist <- distance(coords[i, ], coords[j, ])
      if (dist > max_distance) {
        max_distance <- dist
        point1 <- coords[i, ]
        point2 <- coords[j, ]
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
wmu_number <- wmu_number_list[1]
input_path <- here("Output", "Density", paste0("density", wmu_number, ".RData"))
load(file = input_path)

# load cover
output_path <- here("Output", "Simulation", paste0("cover-WMU", wmu_number, ".RData"))
load(file = output_path)

output_path <- here("Output", "Simulation", paste0("cover-WMU", wmu_number, "-varEffort.RData"))
# output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number,"-T",IMAGE_WIDTH,"H-SG-DF", detectF@key.function, ".RData"))
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
# cover <- make.coverage(region,
#                        spacing = COV_SPACE # OR
#                        # n.grid.points = 1000
# )
# plot(region, cover)


# Define survey design
# Baseline design
l10_design <- make.design(
  region = region,
  transect.type = "line",
  design = "systematic",
  samplers = 10,
  design.angle = 0, # align transect with north south
  edge.protocol = "minus",
  truncation = 260, # IMAGE_WIDTH
  coverage.grid = cover
)
l10_transects <- generate.transects(l10_design)
l10_design <- run.coverage(l10_design, reps = COV_REPS)

l20_design <- make.design(
  region = region,
  transect.type = "line",
  design = "systematic",
  samplers = 20,
  design.angle = 0, # align transect with north south
  edge.protocol = "minus",
  truncation = 260, # IMAGE_WIDTH
  coverage.grid = cover
)
l20_transects <- generate.transects(l20_design)

### Coverage
#### You can re-run the coverage simulation using the following code. Note, your
#### results should vary slightly from mine, make sure you haven't set a seed!
l10_design <- run.coverage(l10_design, reps = COV_REPS)
l20_design <- run.coverage(l20_design, reps = COV_REPS)
total_length <- H_SG_design@design.statistics$line.length[2]
# total_length <- H-SG_transects@line.length


## Systematic subplot design
total_length <- H_SG_design@design.statistics$line.length[2]
number_blocks <- round(total_length / 26000)
print(number_blocks)
# spacing <- 400
spacing <- 260 * 2

poly_dim <- find_best_block_dim(total_length, number_blocks, spacing) # poly_dim, 7 lines, 3715x3640m

# less effort
# Loop to reduce number_blocks by a quarter each time
# Create a list of block counts by reducing the original number of blocks by a factor of 2^(i-1)
block_counts <- lapply(1:4, function(i) round(number_blocks * (1 / (2^(i - 1)))))
block_counts[[5]] <- 46

# Iterate over the list of block counts
for (current_number_blocks in block_counts) {
  print(current_number_blocks)

  # Check if the current number of blocks fits within the region's area
  if (region@area > (poly_dim$x_length *
    poly_dim$y_length *
    current_number_blocks)) {
    # Create a coverage grid with the current number of blocks
    grid_center <- make.coverage(region, n.grid.points = current_number_blocks)
    grid_center@grid <- grid_center@grid %>% select(-coverage.scores)

    # Create a buffer around the polygon boundary to exclude points near the edges
    buffer_1000m <- st_buffer(region@region, dist = -1001)

    # Identify points within 1000 meters of the boundary
    selection_index <- st_disjoint(grid_center@grid, buffer_1000m, sparse = FALSE)
    points_within_1000m <- grid_center@grid[selection_index, ]

    # Check if there are points within the boundary
    if (exists("grid_center") && !is.null(grid_center@grid) && exists("selection_index") && !is.null(selection_index) && nrow(points_within_1000m) > 0) {
      # Find the nearest points on the polygon for each point in the grid
      nearest_lines <- st_nearest_points(points_within_1000m, buffer_1000m)
      nearest_points <- st_cast(nearest_lines, "POINT")
      new_coords <- st_coordinates(nearest_points[seq(2, length(nearest_points), 2)])

      # Replace the coordinates in points1 with the new coordinates
      st_geometry(grid_center@grid)[selection_index] <- st_sfc(lapply(seq_len(nrow(points_within_1000m)), function(i) st_point(new_coords[i, ])))
    }

    # Create sub plot polygons using the adjusted grid
    polygons <- create_sf_polygons(grid_center@grid, poly_dim$x_length, poly_dim$y_length)
    polygons <- st_intersection(polygons, wmu)
    polygons <- st_cast(polygons, "POLYGON")
  }

  # Create a subplot region using the generated polygons
  QC_plots <- make.region(
    region.name = "study area",
    shape = polygons,
    strata.name = polygons$ID
  )

  # Create flight lines within the quadcopter plots
  assign(paste0("QC_Sys_design_", current_number_blocks), make.design(
    region = QC_plots,
    transect.type = "line",
    design = "systematic",
    samplers = numeric(0),
    line.length = numeric(0),
    spacing = spacing,
    design.angle = 0,
    edge.protocol = "minus",
    truncation = 260,
    coverage.grid = cover
  ))

  # Generate transects for the current design
  assign(paste0("QC_Sys_transects_", current_number_blocks), generate.transects(get(paste0("QC_Sys_design_", current_number_blocks))))

  # Run coverage simulation for the current design
  assign(paste0("QC_Sys_design_", current_number_blocks), run.coverage(get(paste0("QC_Sys_design_", current_number_blocks)), reps = COV_REPS))
}


# hist(get.coverage(H_SG_design))
# Plot desings
par(mfrow = c(2, 3))
plot(region, QC_Sys_transects_8, lwd = 0.5, col = 4)
plot(region, QC_Sys_transects_15, lwd = 0.5, col = 4)
plot(region, QC_Sys_transects_30, lwd = 0.5, col = 4)
plot(region, QC_Sys_transects_46, lwd = 0.5, col = 4)
plot(region, QC_Sys_transects_61, lwd = 0.5, col = 4)
par(mfrow = c(1, 1))
par(mfrow = c(2, 3))
plot(QC_Sys_design_8)
plot(QC_Sys_design_15)
plot(QC_Sys_design_30)
plot(QC_Sys_design_46)
plot(QC_Sys_design_61)
par(mfrow = c(1, 1))

## design stats
## For details see: https://examples.distancesampling.org/dssd-getting-started/GettingStarted-distill.html#appendix-trackline-and-cyclic-trackline-lengths
## Line lenght = on effort line length
## The trackline length is the sum of the lengths of the transects plus the off-effort transit time required to complete the survey from the beginning of the first transect to the end of the last transect. The off-effort transit distance is calculated as the crow flies and may be longer in reality if transit is required around lakes, islands or coastlines etc.
## The cyclic trackline length is the trackline length plus the off-effort transit distance required to return from the end of the last transect to the beginning of the first transect.
# Extract key metrics from each simulation summary
l10_design_metric <- extract_design_metrics(l10_design)
l20_design_metric <- extract_design_metrics(l20_design)
QC_Sys_design_8_metric <- extract_design_metrics(QC_Sys_design_8)
QC_Sys_design_15_metric <- extract_design_metrics(QC_Sys_design_15)
QC_Sys_design_30_metric <- extract_design_metrics(QC_Sys_design_30)
QC_Sys_design_46_metric <- extract_design_metrics(QC_Sys_design_46)
QC_Sys_design_61_metric <- extract_design_metrics(QC_Sys_design_61)

# Combine metrics into a single dataframe
design_comparison_df <- data.frame(
  Simulation = c("10L", "20L", "QC-8", "QC-15", "QC-30", "QC-46", "QC-61"),
  Design = c(
    l10_design_metric$design_type[1],
    l20_design_metric$design_type[1],
    QC_Sys_design_8_metric$design_type[1],
    QC_Sys_design_15_metric$design_type[1],
    QC_Sys_design_30_metric$design_type[1],
    QC_Sys_design_46_metric$design_type[1],
    QC_Sys_design_61_metric$design_type[1]
  ),
  Mean_Sampler_Count = c(
    l10_design_metric$mean_sampler_count,
    l20_design_metric$mean_sampler_count,
    QC_Sys_design_8_metric$mean_sampler_count,
    QC_Sys_design_15_metric$mean_sampler_count,
    QC_Sys_design_30_metric$mean_sampler_count,
    QC_Sys_design_46_metric$mean_sampler_count,
    QC_Sys_design_61_metric$mean_sampler_count
  ),
  Mean_Cover_Area = c(
    l10_design_metric$mean_cover_area,
    l20_design_metric$mean_cover_area,
    QC_Sys_design_8_metric$mean_cover_area,
    QC_Sys_design_15_metric$mean_cover_area,
    QC_Sys_design_30_metric$mean_cover_area,
    QC_Sys_design_46_metric$mean_cover_area,
    QC_Sys_design_61_metric$mean_cover_area
  ),
  Mean_Cover_Percentage = c(
    l10_design_metric$mean_cover_percentage,
    l20_design_metric$mean_cover_percentage,
    QC_Sys_design_8_metric$mean_cover_percentage,
    QC_Sys_design_15_metric$mean_cover_percentage,
    QC_Sys_design_30_metric$mean_cover_percentage,
    QC_Sys_design_46_metric$mean_cover_percentage,
    QC_Sys_design_61_metric$mean_cover_percentage
  ),
  Mean_Line_Length = c(
    l10_design_metric$mean_line_length,
    l20_design_metric$mean_line_length,
    QC_Sys_design_8_metric$mean_line_length,
    QC_Sys_design_15_metric$mean_line_length,
    QC_Sys_design_30_metric$mean_line_length,
    QC_Sys_design_46_metric$mean_line_length,
    QC_Sys_design_61_metric$mean_line_length
  ),
  Mean_Trackline_Length = c(
    l10_design_metric$mean_trackline,
    l20_design_metric$mean_trackline,
    QC_Sys_design_8_metric$mean_trackline,
    QC_Sys_design_15_metric$mean_trackline,
    QC_Sys_design_30_metric$mean_trackline,
    QC_Sys_design_46_metric$mean_trackline,
    QC_Sys_design_61_metric$mean_trackline
  ),
  Mean_Cyclic_Trackline_Length = c(
    l10_design_metric$mean_cyclic_trackline,
    l20_design_metric$mean_cyclic_trackline,
    QC_Sys_design_8_metric$mean_cyclic_trackline,
    QC_Sys_design_15_metric$mean_cyclic_trackline,
    QC_Sys_design_30_metric$mean_cyclic_trackline,
    QC_Sys_design_46_metric$mean_cyclic_trackline,
    QC_Sys_design_61_metric$mean_cyclic_trackline
  ),
  Mean_On_Effort = c(
    l10_design_metric$mean_on_effort,
    l20_design_metric$mean_on_effort,
    QC_Sys_design_8_metric$mean_on_effort,
    QC_Sys_design_15_metric$mean_on_effort,
    QC_Sys_design_30_metric$mean_on_effort,
    QC_Sys_design_46_metric$mean_on_effort,
    QC_Sys_design_61_metric$mean_on_effort
  ),
  Mean_Off_Effort = c(
    l10_design_metric$mean_off_effort,
    l20_design_metric$mean_off_effort,
    QC_Sys_design_8_metric$mean_off_effort,
    QC_Sys_design_15_metric$mean_off_effort,
    QC_Sys_design_30_metric$mean_off_effort,
    QC_Sys_design_46_metric$mean_off_effort,
    QC_Sys_design_61_metric$mean_off_effort
  ),
  Mean_Return_to_Home = c(
    l10_design_metric$mean_return2home,
    l20_design_metric$mean_return2home,
    QC_Sys_design_8_metric$mean_return2home,
    QC_Sys_design_15_metric$mean_return2home,
    QC_Sys_design_30_metric$mean_return2home,
    QC_Sys_design_46_metric$mean_return2home,
    QC_Sys_design_61_metric$mean_return2home
  ),
  Mean_Off_Effort_Return = c(
    l10_design_metric$mean_off_effort_return,
    l20_design_metric$mean_off_effort_return,
    QC_Sys_design_8_metric$mean_off_effort_return,
    QC_Sys_design_15_metric$mean_off_effort_return,
    QC_Sys_design_30_metric$mean_off_effort_return,
    QC_Sys_design_46_metric$mean_off_effort_return,
    QC_Sys_design_61_metric$mean_off_effort_return
  ),
  On_Effort_Percentage = c(
    l10_design_metric$on_effort_percentage,
    l20_design_metric$on_effort_percentage,
    QC_Sys_design_8_metric$on_effort_percentage,
    QC_Sys_design_15_metric$on_effort_percentage,
    QC_Sys_design_30_metric$on_effort_percentage,
    QC_Sys_design_46_metric$on_effort_percentage,
    QC_Sys_design_61_metric$on_effort_percentage
  ),
  Off_Effort_Percentage = c(
    l10_design_metric$off_effort_percentage,
    l20_design_metric$off_effort_percentage,
    QC_Sys_design_8_metric$off_effort_percentage,
    QC_Sys_design_15_metric$off_effort_percentage,
    QC_Sys_design_30_metric$off_effort_percentage,
    QC_Sys_design_46_metric$off_effort_percentage,
    QC_Sys_design_61_metric$off_effort_percentage
  ),
  Return_to_Home_Percentage = c(
    l10_design_metric$return2home_percentage,
    l20_design_metric$return2home_percentage,
    QC_Sys_design_8_metric$return2home_percentage,
    QC_Sys_design_15_metric$return2home_percentage,
    QC_Sys_design_30_metric$return2home_percentage,
    QC_Sys_design_46_metric$return2home_percentage,
    QC_Sys_design_61_metric$return2home_percentage
  ),
  Off_Effort_Return_Percentage = c(
    l10_design_metric$off_effort_return_percentage,
    l20_design_metric$off_effort_return_percentage,
    QC_Sys_design_8_metric$off_effort_return_percentage,
    QC_Sys_design_15_metric$off_effort_return_percentage,
    QC_Sys_design_30_metric$off_effort_return_percentage,
    QC_Sys_design_46_metric$off_effort_return_percentage,
    QC_Sys_design_61_metric$off_effort_return_percentage
  ),
  Number_of_Plots = c(
    length(l10_design_metric$design_type),
    length(l20_design_metric$design_type),
    length(QC_Sys_design_8_metric$design_type),
    length(QC_Sys_design_15_metric$design_type),
    length(QC_Sys_design_30_metric$design_type),
    length(QC_Sys_design_46_metric$design_type),
    length(QC_Sys_design_61_metric$design_type)
  )
)
# Print the comparison dataframe
# print(comparison_df)
# Display a comparison table of design metrics for different survey designs
# This table summarizes key metrics such as mean sampler count, coverage area,
# line length, and effort percentages for each design.
kable(design_comparison_df)
# drop all but the first 8 rows


# Save simulation data
output_path <- here("Output", "Simulation", paste0("cover-WMU", wmu_number, "-varEffort.RData"))
save(l10_design, l20_design, QC_Sys_design_8, QC_Sys_design_15, QC_Sys_design_30, QC_Sys_design_46, QC_Sys_design_61, design_comparison_df, file = output_path)

# save comparison_df
output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number, "-varEffort-comparsiondf.csv"))
write.csv(design_comparison_df, file = output_path, row.names = FALSE)
