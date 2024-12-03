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
  if (ALTITUDE <= 0 || CAMERA_FOV <= 0 || CAMERA_ANGLE < 0) {
    stop("Altitude and camera FOV must be positive numbers, and camera angle must be non-negative")
  }

  # Adjust the FOV by doubling the adjustment and adding it to the FOV
  adjusted_FOV <- CAMERA_FOV + 2 * CAMERA_ANGLE

  # Calculate and round the image width
  round(2 * ALTITUDE * tan((adjusted_FOV * pi / 180) / 2), -1)
}


## design stats
## For details see: https://examples.distancesampling.org/dssd-getting-started/GettingStarted-distill.html#appendix-trackline-and-cyclic-trackline-lengths
## Line lenght = on effort line length
## The trackline length is the sum of the lengths of the transects plus the off-effort transit time required to complete the survey from the beginning of the first transect to the end of the last transect. The off-effort transit distance is calculated as the crow flies and may be longer in reality if transit is required around lakes, islands or coastlines etc.
## The cyclic trackline length is the trackline length plus the off-effort transit distance required to return from the end of the last transect to the beginning of the first transect.
# Extract key metrics from each simulation summary
extract_design_metrics <- function(design) {
  list(
    design_type = design@design,
    mean_sampler_count = design@design.statistics$sampler.count[2],
    mean_cover_area = design@design.statistics$cov.area[2],
    mean_cover_percentage = design@design.statistics$p.cov.area[2],
    mean_line_length = design@design.statistics$line.length[2],
    mean_trackline = design@design.statistics$trackline[2],
    mean_cyclic_trackline = design@design.statistics$cyclictrackline[2],
    mean_on_effort = design@design.statistics$line.length[2],
    mean_off_effort = design@design.statistics$trackline[2] - design@design.statistics$line.length[2],
    mean_return2home = design@design.statistics$cyclictrackline[2] - design@design.statistics$trackline[2],
    mean_off_effort_return_percentage = design@design.statistics$cyclictrackline[2] - design@design.statistics$line.length[2],
    on_effort_percentage = (design@design.statistics$line.length[2] / design@design.statistics$cyclictrackline[2]) * 100,
    off_effort_percentage = ((design@design.statistics$trackline[2] - design@design.statistics$line.length[2]) / design@design.statistics$cyclictrackline[2]) * 100,
    return2home_percentage = ((design@design.statistics$cyclictrackline[2] - design@design.statistics$trackline[2]) / design@design.statistics$cyclictrackline[2]) * 100,
    off_effort_return_percentage = ((design@design.statistics$cyclictrackline[2] - design@design.statistics$line.length[2]) / design@design.statistics$cyclictrackline[2]) * 100
  )
}

# Function to find optimal polygon dimensions
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


# Function to generate a random square polygon with random orientation
generate_random_polygon <- function(x, y, crs, area) {
  center_x <- runif(1, min(st_bbox(area)$xmin), max(st_bbox(area)$xmax))
  center_y <- runif(1, min(st_bbox(area)$ymin), max(st_bbox(area)$ymax))

  coords <- matrix(c(-x / 2, -y / 2, x / 2, -y / 2, x / 2, y / 2, -x / 2, y / 2, -x / 2, -y / 2), ncol = 2, byrow = TRUE)

  coords[, 1] <- coords[, 1] + center_x
  coords[, 2] <- coords[, 2] + center_y

  polygon <- st_polygon(list(coords))
  list(polygon = st_sfc(polygon, crs = crs))
}

# Function to check if the polygon is within acceptable bounds with a buffer
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

# Function to check if the polygon overlaps with any existing polygons
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

# Function to place N polygons within the area
place_polygons <- function(N, x, y, area, buffer_distance, max_attempts = 10000, max_retries = 100) {
  crs <- st_crs(area)

  for (retry in 1:max_retries) {
    print(paste0("Retry ", retry))
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
        result <- generate_random_polygon(x, y, crs, area)
        polygon <- result$polygon
        if (check_polygon(polygon, area, buffer_distance) && !check_overlap(polygon, polygons)) {
          polygons[[i]] <- polygon
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
  }

  stop("Failed to place polygons after 100 retries")
}

# create polygons when grid of centroids is provided
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
  polygons <- lapply(1:nrow(center_points), function(i) {
    create_polygon(center_points[i, ], x_dim, y_dim)
  })
  
   # Generate chess-like IDs
  generate_id <- function(index) {
    row <- ceiling(index / 8)
    col <- index %% 8
    if (col == 0) col <- 8
    paste0(LETTERS[row], col)
  }
  
  ids <- sapply(1:nrow(center_points), generate_id)
  
  # Combine polygons into an sf object and add IDs
  sf_polygons <- st_sf(
    geometry = st_sfc(polygons),
    crs = crs,
    ID = ids
  )
  
  sf_polygons
}



# Function to calculate distance between two points
distance <- function(p1, p2) {
  sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
}

# Find the longest dimension and its angle
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

# Extract key metrics from each simulation summary
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
wmu_number_list <- c("501", "503", "512", "528") #' 517'
wmu_number <- wmu_number_list[2]
input_path <- here("Output", "Density", paste0("density", wmu_number, ".RData"))
load(file = input_path)

# Define constants
ALTITUDE <- 120 # Height in meters
CAMERA_HFOV <- 25 # Horizontal FOV in degrees. Max adjustment of 25 degrees. If more than 25 degrees then a third camera or gimbal would be needed to cover 0. Proposed intervals for adjustments are 0, 10, 20, 25
CAMERA_ANGLE <- 0 # Adjustment in degrees; max 25 with two fix cameras or 35 with gimbal in forested
IMAGE_WIDTH <- calculate_image_width(ALTITUDE, CAMERA_HFOV, CAMERA_ANGLE)
print(paste0("Half swath width is: ", IMAGE_WIDTH, " m"))

# Define cover grid spacing and repetition
COV_SPACE = 1000
COV_REPS = 10

# Calculate the angle of longest Dimension of WMU
COORDS <- st_coordinates(region@region)[, 1:2]
LONGEST_DIMENSION <- find_longest_dimension_angle(COORDS)
TRANSECT_ANGLE <- LONGEST_DIMENSION$angle


# Create population description
pop_desc <- make.population.description(
  region = region,
  density = density,
  N = total_abundance,
  fixed.N = T
)

# Define and visualise detection function
detect_hr_overview <- make.detectability(
  key.function = "hr",
  scale.param = rep(200, 10),
  shape.param = seq(0.1, 4, 0.4),
  truncation = IMAGE_WIDTH
)
COLORS <- brewer.pal(10, "Paired")
plot(detect_hr_overview, pop_desc, col = COLORS)
legend(x = "topright", legend = seq(0.1, 4, 0.4), col = COLORS, lty = 1, cex = 0.8)

detect_hr <- make.detectability(
  key.function = "hr",
  scale.param = 100, # heli:50
  shape.param = 3, # heli:2
  truncation = IMAGE_WIDTH # heli:20
)
plot(detect_hr, pop_desc, legend= FALSE)

# Define and visualise uniform detection function
detect_uf <- make.detectability(
  key.function = "uf",
  scale.param = 0.8, # accounting for canopy cover
  truncation = IMAGE_WIDTH
)
plot(detect_uf, pop_desc)

# Set detection function for simulation
detectF <- detect_uf

# create coverage grid
cover <- make.coverage(region,
  spacing = COV_SPACE # OR
  # n.grid.points = 1000
)
# plot(region, cover)


# Define survey design
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
  truncation = IMAGE_WIDTH, # IMAGE_WIDTH
  coverage.grid = cover
)
heli_transects <- generate.transects(heli_design)


### Coverage
#### You can re-run the coverage simulation using the following code. Note, your
#### results should vary slightly from mine, make sure you haven't set a seed!
heli_design <- run.coverage(heli_design, reps = COV_REPS)
total_length <- heli_design@design.statistics$line.length[2]


## Systematic design
sys_design <- make.design(
  region = region,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = total_length, # OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = IMAGE_WIDTH, # IMAGE_WIDTH
  coverage.grid = cover
)
sys_transects <- generate.transects(sys_design)

sys_design <- run.coverage(sys_design, reps = COV_REPS)


## Random design
rnd_design <- make.design(
  region = region,
  transect.type = "line",
  design = "random",
  samplers = numeric(0), # OR
  line.length = total_length, # OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = IMAGE_WIDTH, # IMAGE_WIDTH
  coverage.grid = cover
)
rnd_transects <- generate.transects(rnd_design)
rnd_design <- run.coverage(rnd_design, reps = COV_REPS)



## Zigzag design
zigzag_design <- make.design(
  region = region,
  transect.type = "line",
  design = "eszigzag",
  samplers = numeric(0), # OR
  line.length = total_length, # OR
  spacing = numeric(0),
  design.angle = TRANSECT_ANGLE, # The design angle for the zigzag designs refers to the angle of a line which would run through the middle of each zigzag transect if the zigzags were to be generated within a rectangle. The design angle for zigzags should usually run along the longest dimension of the study region.
  edge.protocol = "minus",
  bounding.shape = "convex.hull", # rectangle or convex.hull. convex hull is generally more efficient.
  truncation = IMAGE_WIDTH, # IMAGE_WIDTH
  coverage.grid = cover
)
zigzag_transects <- generate.transects(zigzag_design)

### Coverage
zigzag_design <- run.coverage(zigzag_design, reps = COV_REPS)


## Zigzag with complementary line
zigzagcom_design <- make.design(
  region = region,
  transect.type = "line",
  design = "eszigzagcom", # eszigzag or eszigzagcom
  samplers = numeric(0), # OR
  line.length = total_length, # OR
  spacing = numeric(0),
  design.angle = TRANSECT_ANGLE,
  edge.protocol = "minus",
  bounding.shape = "convex.hull",
  truncation = IMAGE_WIDTH, # IMAGE_WIDTH
  coverage.grid = cover
)
zigzagcom_transects <- generate.transects(zigzagcom_design)
### Coverage
zigzagcom_design <- run.coverage(zigzagcom_design, reps = COV_REPS)

# Drone survey designs
## Fix-wing
# Compute polygon dimensions
number_blocks <- round(total_length/367200) # 367.2km is the total distance superwake can fly, assuming a speed of 17m/s and a flight time of 6h.
spacing <- 500
poly_dim <- find_best_block_dim(total_length, number_blocks, spacing)

# Checking that blocks fit within region
if(region@area> (poly_dim[2]*poly_dim[3]*number_blocks)){
  # Create polygons
  fixW_poly <- place_polygons(number_blocks, poly_dim$x_length, poly_dim$y_length, wmu, buffer_distance = 100)

}

# create subplot region
fixW_plots <- make.region(
  region.name = "study area",
  shape = fixW_poly,
  strata.name = fixW_poly$ID
)

# plot(fixW_plots)
# create systematic flight lines withing fixed wing
fixW_sys_design <- make.design(
  region = fixW_plots,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = rep(total_length / length(fixW_plots@strata.name), length(fixW_plots@strata.name)), # OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = IMAGE_WIDTH, # IMAGE_WIDTH
  coverage.grid = cover
)
fixW_sys_transects <- generate.transects(fixW_sys_design)
### Coverage
fixW_sys_design <- run.coverage(fixW_sys_design, reps = COV_REPS)




# Fixed wing zigzag flights
fixW_zigzag_design <- make.design(
  region = fixW_plots,
  transect.type = "line",
  design = "eszigzag",
  samplers = numeric(0), # OR
  line.length = rep(total_length / length(fixW_plots@strata.name), length(fixW_plots@strata.name)), # OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = IMAGE_WIDTH, # IMAGE_WIDTH
  coverage.grid = cover
)
fixW_zigzag_transects <- generate.transects(fixW_zigzag_design)
### Coverage
fixW_zigzag_design <- run.coverage(fixW_zigzag_design, reps = COV_REPS)



## Quadcopter
## Systematic subplot design
number_blocks <- round(total_length/26000)
spacing <- 200

# Checking that blocks fit within region
if(region@area> (2000*2500*number_blocks)){
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
  new_coords <- st_coordinates(nearest_points[seq(2,length(nearest_points),2)])

  # Replace the coordinates in points1 with the new coordinates
  st_geometry(grid_center@grid)[selection_index] <- st_sfc(lapply(1:nrow(points_within_1000m), function(i) st_point(new_coords[i, ])))

  # Create sub plot polygons
  polygons <- create_sf_polygons(grid_center@grid, 2500, 2000)
}
# plot(st_geometry(wmu))
# plot(polygons[1], add = TRUE, col = "red")

# create subplot region
quadcopter_plots <- make.region(
  region.name = "study area",
  shape = polygons,
  strata.name = polygons$ID
)

# plot(quadcopter_plots)
# create flight lines withing quadcopterplots
quadcopter_design <- make.design(
  region = quadcopter_plots,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = rep(total_length / length(quadcopter_plots@strata.name), length(quadcopter_plots@strata.name)), # OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = IMAGE_WIDTH, # IMAGE_WIDTH
  coverage.grid = cover
)
quadcopter_transects <- generate.transects(quadcopter_design)
### Coverage
system.time(quadcopter_design <- run.coverage(quadcopter_design, reps = COV_REPS))



# Plot desings
par(mfrow = c(2, 4))
plot(region, heli_transects, lwd = 0.5, col = 4)
plot(region, sys_transects, lwd = 0.5, col = 4)
plot(region, rnd_transects, lwd = 0.5, col = 4)
plot(region, zigzag_transects, lwd = 0.5, col = 4)
plot(region, zigzagcom_transects, lwd = 0.5, col = 4)
plot(region, fixW_sys_transects, lwd = 0.5, col = 4)
plot(region, fixW_zigzag_transects, lwd = 0.5, col = 4)
plot(region, quadcopter_transects, lwd = 0.5, col = 4)
par(mfrow = c(1, 1))
par(mfrow = c(2, 4))
plot(heli_design)
plot(sys_design)
plot(rnd_design)
plot(zigzag_design)
plot(zigzagcom_design)
plot(fixW_sys_design)
plot(fixW_zigzag_design)
plot(quadcopter_design)
par(mfrow = c(1, 1))

## design stats
## For details see: https://examples.distancesampling.org/dssd-getting-started/GettingStarted-distill.html#appendix-trackline-and-cyclic-trackline-lengths
## Line lenght = on effort line length
## The trackline length is the sum of the lengths of the transects plus the off-effort transit time required to complete the survey from the beginning of the first transect to the end of the last transect. The off-effort transit distance is calculated as the crow flies and may be longer in reality if transit is required around lakes, islands or coastlines etc.
## The cyclic trackline length is the trackline length plus the off-effort transit distance required to return from the end of the last transect to the beginning of the first transect.
# Extract key metrics from each simulation summary
heli_design_metric <- extract_design_metrics(heli_design)
sys_design_metric <- extract_design_metrics(sys_design)
rnd_design_metric <- extract_design_metrics(rnd_design)
zigzag_design_metric <- extract_design_metrics(zigzag_design)
zigzagcom_design_metric <- extract_design_metrics(zigzagcom_design)
fixW_sys_design_metric <- extract_design_metrics(fixW_sys_design)
fixW_zigzag_design_metric <- extract_design_metrics(fixW_zigzag_design)
quadcopter_design_metric <- extract_design_metrics(quadcopter_design)

# Combine metrics into a single dataframe
rm(design_comparison_df)
design_comparison_df <- data.frame(
  Simulation = c("Heli", "Sys", "Rnd", "Zig", "Zagcom", "FixW-Sys", "FixW-Zig", "Quadcopter"),
  Design = c(
    heli_design_metric$design_type,
    sys_design_metric$design_type,
    rnd_design_metric$design_type,
    zigzag_design_metric$design_type,
    zigzagcom_design_metric$design_type,
    fixW_sys_design_metric$design_type,
    fixW_zigzag_design_metric$design_type,
    quadcopter_design_metric$design_type
  ),
  Mean_Sampler_Count = c(
    heli_design_metric$mean_sampler_count,
    sys_design_metric$mean_sampler_count,
    rnd_design_metric$mean_sampler_count,
    zigzag_design_metric$mean_sampler_count,
    zigzagcom_design_metric$mean_sampler_count,
    fixW_sys_design_metric$mean_sampler_count,
    fixW_zigzag_design_metric$mean_sampler_count,
    quadcopter_design_metric$mean_sampler_count
  ),
  Mean_Cover_Area = c(
    heli_design_metric$mean_cover_area,
    sys_design_metric$mean_cover_area,
    rnd_design_metric$mean_cover_area,
    zigzag_design_metric$mean_cover_area,
    zigzagcom_design_metric$mean_cover_area,
    fixW_sys_design_metric$mean_cover_area,
    fixW_zigzag_design_metric$mean_cover_area,
    quadcopter_design_metric$mean_cover_area
  ),
  Mean_Cover_Percentage = c(
    heli_design_metric$mean_cover_percentage,
    sys_design_metric$mean_cover_percentage,
    rnd_design_metric$mean_cover_percentage,
    zigzag_design_metric$mean_cover_percentage,
    zigzagcom_design_metric$mean_cover_percentage,
    fixW_sys_design_metric$mean_cover_percentage,
    fixW_zigzag_design_metric$mean_cover_percentage,
    quadcopter_design_metric$mean_cover_percentage
  ),
  Mean_Line_Length = c(
    heli_design_metric$mean_line_length,
    sys_design_metric$mean_line_length,
    rnd_design_metric$mean_line_length,
    zigzag_design_metric$mean_line_length,
    zigzagcom_design_metric$mean_line_length,
    fixW_sys_design_metric$mean_line_length,
    fixW_zigzag_design_metric$mean_line_length,
    quadcopter_design_metric$mean_line_length
  ),
  Mean_Trackline_Length = c(
    heli_design_metric$mean_trackline,
    sys_design_metric$mean_trackline,
    rnd_design_metric$mean_trackline,
    zigzag_design_metric$mean_trackline,
    zigzagcom_design_metric$mean_trackline,
    fixW_sys_design_metric$mean_trackline,
    fixW_zigzag_design_metric$mean_trackline,
    quadcopter_design_metric$mean_trackline
  ),
  Mean_Cyclic_Trackline_Length = c(
    heli_design_metric$mean_cyclic_trackline,
    sys_design_metric$mean_cyclic_trackline,
    rnd_design_metric$mean_cyclic_trackline,
    zigzag_design_metric$mean_cyclic_trackline,
    zigzagcom_design_metric$mean_cyclic_trackline,
    fixW_sys_design_metric$mean_cyclic_trackline,
    fixW_zigzag_design_metric$mean_cyclic_trackline,
    quadcopter_design_metric$mean_cyclic_trackline
  ),
  Mean_On_Effort = c(
    heli_design_metric$mean_on_effort,
    sys_design_metric$mean_on_effort,
    rnd_design_metric$mean_on_effort,
    zigzag_design_metric$mean_on_effort,
    zigzagcom_design_metric$mean_on_effort,
    fixW_sys_design_metric$mean_on_effort,
    fixW_zigzag_design_metric$mean_on_effort,
    quadcopter_design_metric$mean_on_effort
  ),
  Mean_Off_Effort = c(
    heli_design_metric$mean_off_effort,
    sys_design_metric$mean_off_effort,
    rnd_design_metric$mean_off_effort,
    zigzag_design_metric$mean_off_effort,
    zigzagcom_design_metric$mean_off_effort,
    fixW_sys_design_metric$mean_off_effort,
    fixW_zigzag_design_metric$mean_off_effort,
    quadcopter_design_metric$mean_off_effort
  ),
  Mean_Return_to_Home = c(
    heli_design_metric$mean_return2home,
    sys_design_metric$mean_return2home,
    rnd_design_metric$mean_return2home,
    zigzag_design_metric$mean_return2home,
    zigzagcom_design_metric$mean_return2home,
    fixW_sys_design_metric$mean_return2home,
    fixW_zigzag_design_metric$mean_return2home,
    quadcopter_design_metric$mean_return2home
  ),
  Mean_Off_Effort_Return = c(
    heli_design_metric$mean_off_effort_return,
    sys_design_metric$mean_off_effort_return,
    rnd_design_metric$mean_off_effort_return,
    zigzag_design_metric$mean_off_effort_return,
    zigzagcom_design_metric$mean_off_effort_return,
    fixW_sys_design_metric$mean_off_effort_return,
    fixW_zigzag_design_metric$mean_off_effort_return,
    quadcopter_design_metric$mean_off_effort_return
  ),
  On_Effort_Percentage = c(
    heli_design_metric$on_effort_percentage,
    sys_design_metric$on_effort_percentage,
    rnd_design_metric$on_effort_percentage,
    zigzag_design_metric$on_effort_percentage,
    zigzagcom_design_metric$on_effort_percentage,
    fixW_sys_design_metric$on_effort_percentage,
    fixW_zigzag_design_metric$on_effort_percentage,
    quadcopter_design_metric$on_effort_percentage
  ),
  Off_Effort_Percentage = c(
    heli_design_metric$off_effort_percentage,
    sys_design_metric$off_effort_percentage,
    rnd_design_metric$off_effort_percentage,
    zigzag_design_metric$off_effort_percentage,
    zigzagcom_design_metric$off_effort_percentage,
    fixW_sys_design_metric$off_effort_percentage,
    fixW_zigzag_design_metric$off_effort_percentage,
    quadcopter_design_metric$off_effort_percentage
  ),
  Return_to_Home_Percentage = c(
    heli_design_metric$return2home_percentage,
    sys_design_metric$return2home_percentage,
    rnd_design_metric$return2home_percentage,
    zigzag_design_metric$return2home_percentage,
    zigzagcom_design_metric$return2home_percentage,
    zigzagcom_design_metric$off_effort_percentage,
    fixW_zigzag_design_metric$return2home_percentage,
    quadcopter_design_metric$return2home_percentage
  ),
  Off_Effort_Return_Percentage = c(
    heli_design_metric$off_effort_return_percentage,
    sys_design_metric$off_effort_return_percentage,
    rnd_design_metric$off_effort_return_percentage,
    zigzag_design_metric$off_effort_return_percentage,
    zigzagcom_design_metric$off_effort_return_percentage,
    fixW_sys_design_metric$off_effort_return_percentage,
    fixW_zigzag_design_metric$off_effort_return_percentage,
    quadcopter_design_metric$off_effort_return_percentage
  )
)
# unknown error: doublicated rows
design_comparison_df <- design_comparison_df[1:8,]
# Print the comparison dataframe
# print(comparison_df)
kable(design_comparison_df)
# drop all but the first 8 rows


# Save simulation data
output_path <- here("Output", "Simulation", paste0("cover-WMU", wmu_number,".RData"))
# output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number,"-T",IMAGE_WIDTH,"heli-DF", detectF@key.function, ".RData"))
save(heli_design, sys_design, rnd_design, zigzag_design, zigzagcom_design, fixW_sys_design, fixW_zigzag_design, quadcopter_design, heli_transects,sys_transects, rnd_transects,zigzag_transects, zigzagcom_transects, fixW_sys_transects, fixW_zigzag_transects, quadcopter_transects, design_comparison_df, file = output_path)



# # Define analysis models
# ddf_analyses <- make.ds.analysis(
#   dfmodel = ~1,
#   key = "hn",
#   criteria = "AIC",
#   truncation = IMAGE_WIDTH
# )
ddf_analyses <- make.ds.analysis(
  dfmodel = list(~1, ~1),
  key = c("hn", "hr"),
  criteria = "AIC",
  truncation = IMAGE_WIDTH
)


# Create and run the simulation
sim_heli <- make.simulation(
  reps = 999,
  design = heli_design,
  population.description = pop_desc,
  detectability = detectF,
  ds.analysis = ddf_analyses
)

sim_rnd <- make.simulation(
  reps = 999,
  design = rnd_design,
  population.description = pop_desc,
  detectability = detectF,
  ds.analysis = ddf_analyses
)

sim_sys <- make.simulation(
  reps = 999,
  design = sys_design,
  population.description = pop_desc,
  detectability = detectF,
  ds.analysis = ddf_analyses
)
sim_zig <- make.simulation(
  reps = 999,
  design = zigzag_design,
  population.description = pop_desc,
  detectability = detectF,
  ds.analysis = ddf_analyses
)
sim_zagcom <- make.simulation(
  reps = 999,
  design = zigzagcom_design,
  population.description = pop_desc,
  detectability = detectF,
  ds.analysis = ddf_analyses
)

# Drone sims
example_population <- generate.population(object = pop_desc, detectability = detect_uf, region = region)
# plot(example_population, region)
# termine abundance in each subplot
# Convert points dataframe to sf object
points_sf <- st_as_sf(example_population@population, coords = c("x", "y"), crs = st_crs(wmu))

# Perform spatial join to count points within each polygon
points_within_fixW_polygons <- st_join(points_sf, fixW_sys_design@region@region, join = st_within)

# Count the number of points in each polygon
fixW_points_count <- points_within_fixW_polygons %>%
  group_by(ID) %>%  # Replace `id` with the actual column name identifying polygons
  summarise(count = n()) %>%
  filter(!is.na(ID))  # Remove rows with NA in the id column

# View the result
print(fixW_points_count)
fixW_points_count$count

fixW_density <- density
fixW_density@density.surface[[1]] <- st_intersection(density@density.surface[[1]], fixW_sys_design@region@region)
fixW_density@region.name <- fixW_sys_design@region@region.name
fixW_density@strata.name <- fixW_sys_design@region@strata.name
fixW_density@density.surface[[1]] <- fixW_density@density.surface[[1]] %>%
  mutate(strata = ID) %>%
  select(-ID)

pop_desc_fixW <- make.population.description(
  region = fixW_plots,
  density = fixW_density,
  N = fixW_points_count$count,
  fixed.N = TRUE
)

ddf_analyses_fixW <- make.ds.analysis(
  dfmodel = list(~1, ~1),
  key = c("hn", "hr"),
  criteria = "AIC",
  truncation = IMAGE_WIDTH,
  group.strata = data.frame(design.id = fixW_plots@strata.name, analysis.id = rep("A", length(fixW_plots@strata.name)))
)

sim_fixW_sys <- make.simulation(
  reps = 999,
  design = fixW_sys_design,
  population.description = pop_desc_fixW,
  detectability = detectF,
  ds.analysis = ddf_analyses_fixW
)

sim_fixW_zigzag <- make.simulation(
  reps = 999,
  design = fixW_zigzag_design,
  population.description = pop_desc_fixW,
  detectability = detectF,
  ds.analysis = ddf_analyses_fixW
)
# termine abundance in each subplot

# Perform spatial join to count points within each polygon
points_within_quad_polygons <- st_join(points_sf, quadcopter_design@region@region, join = st_within)

# Count the number of points in each polygon
quad_points_count <- points_within_quad_polygons %>%
  group_by(ID) %>%  # Replace `id` with the actual column name identifying polygons
  summarise(count = n()) %>%
  filter(!is.na(ID))  # Remove rows with NA in the id column

# View the result
print(quad_points_count)
quad_points_count$count

quadcopter_density <- density
quadcopter_density@density.surface[[1]] <- st_intersection(density@density.surface[[1]], quadcopter_design@region@region)
quadcopter_density@region.name <- quadcopter_design@region@region.name
quadcopter_density@strata.name <- quadcopter_design@region@strata.name
quadcopter_density@density.surface[[1]] <- quadcopter_density@density.surface[[1]] %>%
  mutate(strata = ID) %>%
  select(-ID)


pop_desc_quadcopter <- make.population.description(
  region = quadcopter_design@region,
  density = quadcopter_density,
  N = quad_points_count$count,
  fixed.N = TRUE
)

ddf_analyses_quadcopter <- make.ds.analysis(
  dfmodel = list(~1, ~1),
  key = c("hn", "hr"),
  criteria = "AIC",
  truncation =  IMAGE_WIDTH,
  group.strata = data.frame(design.id = quadcopter_design@region@strata.name, analysis.id = rep("A", length(quadcopter_design@region@strata.name)))
)

sim_quad <- make.simulation(
  reps = 999,
  design = quadcopter_design,
  population.description = pop_desc_quadcopter,
  detectability = detectF,
  ds.analysis = ddf_analyses_quadcopter
)


heli_survey <- run.survey(sim_heli)
sys_survey <- run.survey(sim_sys)
zig_survey <- run.survey(sim_zig)
zagcom_survey <- run.survey(sim_zagcom)
sim_fixW_sys <- run.survey(sim_fixW_sys)
sim_fixW_zigzag <- run.survey(sim_fixW_zigzag)
sim_quad <- run.survey(sim_quad)

# plot(heli_survey, region)
# plot(sys_survey, region)
# plot(zig_survey, region)
# plot(zagcom_survey, region)


# Run the full simulation
sim_heli <- run.simulation(simulation = sim_heli, run.parallel = T)
sim_sys <- run.simulation(simulation = sim_sys, run.parallel = T)
sim_zig <- run.simulation(simulation = sim_zig, run.parallel = T)
sim_zagcom <- run.simulation(simulation = sim_zagcom, run.parallel = T)

# Save simulation data
output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number,"-T",IMAGE_WIDTH,"-DF", detectF@key.function, ".RData"))
# output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number,"-T",IMAGE_WIDTH,"heli-DF", detectF@key.function, ".RData"))
save(sim_heli, sim_sys, sim_zig, sim_zagcom, file = output_path)

# output_path <- here("Output", "Simulation")
# # save.sim.results(sim_heli, output_path)

# Display results
# summary(sim_heli, description.summary = FALSE)
# summary(sim_sys, description.summary = FALSE)
# summary(sim_zig, description.summary = FALSE)
# summary(sim_zagcom, description.summary = FALSE)

# par(mfrow = c(2, 2))
# histogram.N.ests(sim_heli, xlim = c(7500, 11000))
# histogram.N.ests(sim_sys, xlim = c(7500, 11000))
# histogram.N.ests(sim_zig, xlim = c(7500, 11000))
# histogram.N.ests(sim_zagcom, xlim = c(7500, 11000))


# # Extract metrics for each simulation
# metrics_heli <- extract_metrics(sim_heli)
# metrics_sys <- extract_metrics(sim_sys)
# metrics_zig <- extract_metrics(sim_zig)
# metrics_zagcom <- extract_metrics(sim_zagcom)

# # Combine metrics into a single dataframe
# comparison_df <- data.frame(
#   Simulation = c("Heli", "Sys", "Zig", "Zagcom"),
#   Mean_Estimate = c(metrics_heli$mean_estimate, metrics_sys$mean_estimate, metrics_zig$mean_estimate, metrics_zagcom$mean_estimate),
#   Percent_Bias = c(metrics_heli$percent_bias, metrics_sys$percent_bias, metrics_zig$percent_bias, metrics_zagcom$percent_bias),
#   RMSE = c(metrics_heli$rmse, metrics_sys$rmse, metrics_zig$rmse, metrics_zagcom$rmse),
#   CI_Coverage_Prob = c(metrics_heli$ci_coverage_prob, metrics_sys$ci_coverage_prob, metrics_zig$ci_coverage_prob, metrics_zagcom$ci_coverage_prob),
#   Mean_SE = c(metrics_heli$mean_se, metrics_sys$mean_se, metrics_zig$mean_se, metrics_zagcom$mean_se),
#   SD_of_Means = c(metrics_heli$sd_of_means, metrics_sys$sd_of_means, metrics_zig$sd_of_means, metrics_zagcom$sd_of_means),
#   Mean_Cover_Area = c(metrics_heli$mean_cover_area, metrics_sys$mean_cover_area, metrics_zig$mean_cover_area, metrics_zagcom$mean_cover_area),
#   Mean_Effort = c(metrics_heli$mean_effort, metrics_sys$mean_effort, metrics_zig$mean_effort, metrics_zagcom$mean_effort),
#   Mean_n = c(metrics_heli$mean_n, metrics_sys$mean_n, metrics_zig$mean_n, metrics_zagcom$mean_n),
#   Mean_k = c(metrics_heli$mean_k, metrics_sys$mean_k, metrics_zig$mean_k, metrics_zagcom$mean_k),
#   Mean_ER = c(metrics_heli$mean_ER, metrics_sys$mean_ER, metrics_zig$mean_ER, metrics_zagcom$mean_ER),
#   Mean_se_ER = c(metrics_heli$mean_se_ER, metrics_sys$mean_se_ER, metrics_zig$mean_se_ER, metrics_zagcom$mean_se_ER)
# )

# # Print the comparison dataframe
# # print(comparison_df)
# kable(comparison_df)






# ## Testing truncation distances. Required with drone???
# # Investigate truncation distances
# truncation_distances <- c(
#   calculate_image_width(100), calculate_image_width(200),
#   calculate_image_width(300), calculate_image_width(400),
#   calculate_image_width(500)
# )

# results_list <- vector("list", length(truncation_distances))
# summary_list <- vector("list", length(truncation_distances))

# for (i in seq_along(truncation_distances)) {
#   cat(sprintf("\nRunning for truncation = %d", truncation_distances[i]))

#   new_ds_analyses <- make.ds.analysis(
#     dfmodel = list(~1, ~1),
#     key = c("hn", "hr"),
#     criteria = "AIC",
#     truncation = truncation_distances[i]
#   )

#   sim@ds.analysis <- new_ds_analyses
#   results_list[[i]] <- run.simulation(sim, run.parallel = F)
#   summary_list[[i]] <- summary(results_list[[i]], description.summary = FALSE)
# }

# names(results_list) <- paste0("t", truncation_distances)
# names(summary_list) <- paste0("t", truncation_distances)





# # Extracting results statistics

# N <- unlist(lapply(summary_list, function(x) {
#   x@individuals$N$mean.Estimate
# }))
# n <- unlist(lapply(summary_list, function(x) {
#   x@individuals$summary$mean.n
# }))
# se <- unlist(lapply(summary_list, function(x) {
#   x@individuals$N$mean.se
# }))
# sd_N <- unlist(lapply(summary_list, function(x) {
#   x@individuals$N$sd.of.means
# }))
# bias <- unlist(lapply(summary_list, function(x) {
#   x@individuals$N$percent.bias
# }))
# RMSE <- unlist(lapply(summary_list, function(x) {
#   x@individuals$N$RMSE
# }))
# cov <- unlist(lapply(summary_list, function(x) {
#   x@individuals$N$CI.coverage.prob
# }))

# sim_data <- data.frame(
#   trunc = truncation_distances,
#   n = round(n),
#   N = round(N),
#   se = round(se, 2),
#   sd.N = round(sd_N, 2),
#   bias = round(bias, 2),
#   RMSE = round(RMSE, 2),
#   cov = round(cov * 100, 1)
# )

# kable(sim_data,
#   col.names = c("$Truncation$", "$mean\\ n$", "$mean\\ \\hat{N}$", "$mean\\ se$", "$SD(\\hat{N})$", "$\\% Bias$", "$RMSE$", "$\\%\\ CI\\ Coverage$"),
#   row.names = FALSE,
#   align = c("c", "c", "c", "c", "c", "c", "c", "c"),
#   caption = "Simulation Results for the simple half normal detection probability: The truncation distance, mean number of detections, mean estimated population size (N), mean standard error of $\\hat{N}$, the standard deviation of $\\hat{N}$, percentage bias, root mean squared error, percentage of times the true value of N was captured in the confidence intervals.",
#   table.placement = "!h",
#   format = "simple"
# )
