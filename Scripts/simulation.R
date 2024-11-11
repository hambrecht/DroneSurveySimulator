# Load necessary libraries
library(here)
library(dsims)
library(knitr)
library(RColorBrewer)
library(dplyr)
library(sf)
library(units)

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
wmu_number <- wmu_number_list[1]
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

# create coverage grid
cover <- make.coverage(region,
  spacing = 1000
  # n.grid.points = 1000
)
plot(region, cover)

detectF <- detect_hr

# Define survey design
## Helicopter design
heli_design <- make.design(
  region = region,
  transect.type = "line",
  design = "segmentedgrid",
  spacing = 1200, # segments seperated by 1.2km
  seg.length = 10000, # segements of 10km
  design.angle = correct_degrees(0), # align transect with north south
  seg.threshold = 10, # any segments less than 10% of the segment length (i.e. 1km) will be discarded.
  edge.protocol = "minus",
  truncation = 600, # IMAGE_WIDTH
  coverage.grid = cover
)
heli_transects <- generate.transects(heli_design)
# plot(region, heli_transects, lwd = 0.5, col = 4)
### Coverage
#### You can re-run the coverage simulation using the following code. Note, your
#### results should vary slightly from mine, make sure you haven't set a seed!
heli_design <- run.coverage(heli_design, reps = 10)
# plot(heli_design)

## Systematic design
sys_design <- make.design(
  region = region,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = heli_transects@line.length, # OR
  spacing = numeric(0),
  design.angle = correct_degrees(0),
  edge.protocol = "minus",
  truncation = 600, # IMAGE_WIDTH
  coverage.grid = cover
)
sys_transects <- generate.transects(sys_design)
# plot(region, sys_transects, lwd = 0.5, col = 4)
sys_design <- run.coverage(sys_design, reps = 10)
# plot(sys_design)


## Zigzag design
zigzag_design <- make.design(
  region = region,
  transect.type = "line",
  design = "eszigzag",
  samplers = numeric(0), # OR
  line.length = numeric(0), # OR
  spacing = 1200,
  design.angle = TRANSECT_ANGLE, # The design angle for the zigzag designs refers to the angle of a line which would run through the middle of each zigzag transect if the zigzags were to be generated within a rectangle. The design angle for zigzags should usually run along the longest dimension of the study region.
  edge.protocol = "minus",
  bounding.shape = "convex.hull", # rectangle or convex.hull. convex hull is generally more efficient.
  truncation = 600, # IMAGE_WIDTH
  coverage.grid = cover
)
zigzag_transects <- generate.transects(zigzag_design)
# plot(region, zigzag_transects, lwd = 0.5, col = 4)
### Coverage
zigzag_design <- run.coverage(zigzag_design, reps = 10)
# plot(zigzag_design)

## Zigzag with complementary line
zigzagcom_design <- make.design(
  region = region,
  transect.type = "line",
  design = "eszigzagcom", # eszigzag or eszigzagcom
  samplers = numeric(0), # OR
  line.length = heli_transects@line.length, # OR
  spacing = numeric(0),
  design.angle = TRANSECT_ANGLE,
  edge.protocol = "minus",
  bounding.shape = "convex.hull",
  truncation = 600, # IMAGE_WIDTH
  coverage.grid = cover
)
zigzagcom_transects <- generate.transects(zigzagcom_design)
# plot(region, zigzagcom_transects, lwd = 0.5, col = 4)
### Coverage
zigzagcom_design <- run.coverage(zigzagcom_design, reps = 10)
# plot(zigzagcom_design)




# # Define analysis models
# ddf_analyses <- make.ds.analysis(
#   dfmodel = ~1,
#   key = "hn",
#   criteria = "AIC",
#   truncation = 600
# )
ddf_analyses <- make.ds.analysis(
  dfmodel = list(~1, ~1),
  key = c("hn", "hr"),
  criteria = "AIC",
  truncation = 600
)



# Create and run the simulation
sim_heli <- make.simulation(
  reps = 999,
  design = heli_design,
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


heli_survey <- run.survey(sim_heli)
sys_survey <- run.survey(sim_sys)
zig_survey <- run.survey(sim_zig)
zagcom_survey <- run.survey(sim_zagcom)

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
# output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number,"-T",IMAGE_WIDTH,"-DF", detectF@key.function, ".RData"))
output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number,"-T",IMAGE_WIDTH,"heli-DF", detectF@key.function, ".RData"))
save(sim_heli, sim_sys, sim_zig, sim_zagcom, file = output_path)

# output_path <- here("Output", "Simulation")
# # save.sim.results(sim_heli, output_path)

# # Display results
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
