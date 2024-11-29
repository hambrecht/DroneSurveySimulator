# Load necessary libraries
library(here)
library(dsims)
library(knitr)
library(RColorBrewer)
library(dplyr)
library(sf)
library(units)
library(geosphere)

#TODO: try replicate https://github.com/DistanceDevelopment/dsims/issues/76

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

correct_degrees <- function(angle) {
  corrected_angle <- (angle - 2) %% 180
  return(corrected_angle)
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
  angle <- runif(1, 0, 180)
  center_x <- runif(1, min(st_bbox(area)$xmin), max(st_bbox(area)$xmax))
  center_y <- runif(1, min(st_bbox(area)$ymin), max(st_bbox(area)$ymax))

  coords <- matrix(c(-x / 2, -y / 2, x / 2, -y / 2, x / 2, y / 2, -x / 2, y / 2, -x / 2, -y / 2), ncol = 2, byrow = TRUE)
  rotation_matrix <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), ncol = 2)
  rotated_coords <- coords %*% rotation_matrix

  rotated_coords[, 1] <- rotated_coords[, 1] + center_x
  rotated_coords[, 2] <- rotated_coords[, 2] + center_y

  polygon <- st_polygon(list(rotated_coords))
  list(polygon = st_sfc(polygon, crs = crs), angle = angle)
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
place_polygons <- function(N, x, y, area, buffer_distance) {
  polygons <- list()
  angles <- numeric(N)
  crs <- st_crs(area)
  for (i in 1:N) {
    repeat {
      result <- generate_random_polygon(x, y, crs, area)
      polygon <- result$polygon
      angle <- result$angle
      if (check_polygon(polygon, area, buffer_distance) && !check_overlap(polygon, polygons)) {
        polygons[[i]] <- polygon
        angles[i] <- angle
        break
      }
    }
  }
  sf_polygons <- st_sf(geometry = do.call(c, polygons), crs = crs)
  sf_polygons$angle <- angles
  sf_polygons$id <- LETTERS[1:N]
  sf_polygons
}

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
  scale.param = rep(150, 10),
  shape.param = seq(0.1, 2, 0.2),
  truncation = 600
)
COLORS <- brewer.pal(10, "Paired")
plot(detect_hr_overview, pop_desc, col = COLORS)
legend(x = "topright", legend = seq(0.1, 2, 0.2), col = COLORS, lty = 1, cex = 0.8)

detect_hr <- make.detectability(
  key.function = "hr",
  scale.param = 150,
  shape.param = 1.9,
  truncation = 600
)
plot(detect_hr, pop_desc)

# Define and visualise uniform detection function
detect_uf <- make.detectability(
  key.function = "uf",
  scale.param = 0.9, # accounting for canopy cover
  truncation = IMAGE_WIDTH
)
plot(detect_uf, pop_desc)

detectF <- detect_uf

# Create example population
example_population <- generate.population(object = pop_desc, detectability = detect_uf, region = region)
plot(example_population, region)

# create coverage grid
cover <- make.coverage(region,
  spacing = 1000
  # n.grid.points = 1000
)
plot(region, cover)

# subsample design

# example design
example_design <- make.design(
  region = region,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = numeric(0), # OR
  spacing = 1200,
  design.angle = correct_degrees(0),
  edge.protocol = "minus",
  truncation = 600 # IMAGE_WIDTH
)
example_transects <- generate.transects(example_design)
# retrieve individuals from helisurvey within suplots to det
# Define parameters
total_length <- example_transects@line.length
number_blocks <- round(total_length/26000)
spacing <- 200


region@area> (2000*2500*number_blocks)
# best_block_dim <- find_best_block_dim(total_length, number_blocks, spacing)
# print(best_block_dim)
# polygons <- place_polygons(number_blocks, best_block_dim$x_length, best_block_dim$y_length, wmu, buffer_distance = 100)
# polygons <- place_polygons(number_blocks, 2500, 2000, wmu, buffer_distance = 100)
# create coverage grid
grid_center <- make.coverage(region,
  # spacing = 1000
  n.grid.points = number_blocks
)
plot(region, grid_center)
grid_center@grid <- grid_center@grid %>% select(-coverage.scores)


# Assuming `polygon` is your sf polygon and `points` is your sf points

# Step 1: Create a buffer around the polygon boundary
buffer_1000m <- st_buffer(region@region, dist = -1001)

rm(points_within_1000m)
# Step 2: Identify points within 1000 meters of the boundary
selection_index <- st_disjoint(grid_center@grid, buffer_1000m, sparse = FALSE)
points_within_1000m <- grid_center@grid[selection_index,]

# Step 3: Find the nearest points on the polygon for each point in the grid
nearest_lines <- st_nearest_points(points_within_1000m, buffer_1000m)
nearest_points <- st_cast(nearest_lines, "POINT")
new_coords <- st_coordinates(nearest_points[seq(2,length(nearest_points),2)])

# Step 3: Replace the coordinates in points1 with the new coordinates
st_geometry(grid_center@grid)[selection_index] <- st_sfc(lapply(1:nrow(points_within_1000m), function(i) st_point(new_coords[i, ])))


polygons <- create_sf_polygons(grid_center@grid, 2500, 2000)

plot(st_geometry(wmu))
plot(polygons[1], add = TRUE, col = "red")

subplots <- make.region(
  region.name = "study area",
  shape = polygons,
  strata.name = polygons$ID
)

plot(subplots)


# termine abundance in each subplot
# Convert points dataframe to sf object
points_sf <- st_as_sf(example_population@population, coords = c("x", "y"), crs = st_crs(polygons))

# Perform spatial join to count points within each polygon
points_within_polygons <- st_join(points_sf, polygons, join = st_within)

# Count the number of points in each polygon
points_count <- points_within_polygons %>%
  group_by(ID) %>%  # Replace `id` with the actual column name identifying polygons
  summarise(count = n()) %>%
  filter(!is.na(ID))  # Remove rows with NA in the id column

# View the result
print(points_count)
points_count$count

empty_density <- make.density(region = subplots, x.space = 500, constant = 0.001)
plot(empty_density, region)

sub_density <- density
sub_density@density.surface[[1]] <- st_intersection(density@density.surface[[1]], subplots@region)
sub_density@region.name <- subplots@region.name
sub_density@strata.name <- subplots@strata.name
sub_density@density.surface[[1]] <- sub_density@density.surface[[1]] %>%
  mutate(strata = ID) %>%
  select(-ID)

plot(sub_density, region)

# Create subpopulation description
sub_pop_desc <- make.population.description(
  region = subplots,
  density = sub_density,
  N = points_count$count,
  fixed.N = TRUE
)

# # create coverage grid
# sub_cover <- make.coverage(subplots,
#   spacing = 1000
#   # n.grid.points = 1000
# )
# plot(subplots, sub_cover)

## Systematic subplot design
subplots_design <- make.design(
  region = subplots,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = numeric(0), # OR
  spacing = spacing,
  design.angle = 0,
  edge.protocol = "minus",
  truncation = IMAGE_WIDTH, # IMAGE_WIDTH
  # coverage.grid = cover
)

subplots_transects <- generate.transects(subplots_design)
plot(region, subplots_transects, lwd = 0.5, col = 4)
# subplots_design <- run.coverage(subplots_design, reps = 10)
# plot(subplots_design)

ddf_analyses_sub <- make.ds.analysis(
  dfmodel = ~1,
  key = c("hn", "hr"),
  criteria = "AIC",
  truncation =  IMAGE_WIDTH,
  group.strata = data.frame(design.id = subplots@strata.name, analysis.id = rep("A", length(subplots@strata.name)))
)

# subplots
sim_sub <- make.simulation(
  reps = 999,
  design = subplots_design,
  population.description = sub_pop_desc,
  detectability = detect_uf,
  ds.analysis = ddf_analyses_sub
)
 # survey
slotNames(sim_sub)
sim_sub@design
summary(sim_sub, use.max.reps = TRUE, description.summary = FALSE)
sub_survey <- run.survey(sim_sub)
sub_survey
plot(sub_survey, subplots)

# Run the full simulation
sim_sub <- run.simulation(simulation = sim_sub, run.parallel = TRUE)

# Results
summary(sim_sub, description.summary = FALSE)
histogram.N.ests(sim_sub, xlim = c(7500, 9500))
histogram.N.ests(sim_sub)
sum(points_count$count)
