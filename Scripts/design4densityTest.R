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



#' Function to calculate distance between two points
#'
#' @param p1 First point.
#' @param p2 Second point.
#' @return The distance between the two points.
distance <- function(p1, p2) {
  sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
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

# Define the 10x10km study area
area_m <- matrix(c(0,0,10000,0,10000,10000,0,10000,0,0), ncol=2, byrow=TRUE)
area <- sf::st_polygon(list(area_m))

# Plot the study area and the polygons
plot(area)

# Create regions for the study area and subplots
region <- make.region(region.name = "study area", shape = area)

# H30T Equivalent Focal Length: 52mm, DFOV: 45.2°, Photo Resolution: 1280×1024, assumed HFOV: 35
# H20T DFOV: 40°, Resolution: 640×512, 
# Mavic 3T DFOV: 61°, Equivalent Focal Length: 40mm, Resolution 640 × 512 
# Define constants
ALTITUDE <- 120 # Height in meters
CAMERA_HFOV <- 35 # Horizontal FOV in degrees. Max adjustment of 25 degrees. If more than 25 degrees then a third camera or gimbal would be needed to cover 0. Proposed intervals for adjustments are 0, 10, 20, 25
CAMERA_ANGLE <- 0 # Adjustment in degrees; max 25 with two fix cameras or 35 with gimbal in forested
IMAGE_WIDTH <- calculate_image_width(ALTITUDE, CAMERA_HFOV, CAMERA_ANGLE)
print(paste0("Half swath width is: ", IMAGE_WIDTH, " m"))

# Define cover grid spacing and repetition
COV_SPACE <- 500
COV_REPS <- 100

# create coverage grid
cover <- make.coverage(region,
                       spacing = COV_SPACE # OR
                       # n.grid.points = 1000
)
plot(region, cover)


# Define survey design
poly_dim <- c()
poly_dim$x_length <- 2000 # 2km x 2km plot
poly_dim$y_length <- 2000 # 2km x 2km plot


# less effort
# Loop to reduce number_blocks by a quarter each time
# Create a list of block counts by reducing the original number of blocks by a factor of 2^(i-1)
block_counts <- (2:5)^2

# Iterate over the list of block counts
for (current_number_blocks in block_counts) {
  print(current_number_blocks)
  
  # Compute grid dimensions
  n_side <- sqrt(current_number_blocks)
  if (n_side != floor(n_side)) stop("Number of plots must be a perfect square.")

  # Calculate spacing between plots
  total_plot_span <- n_side * poly_dim$x_length
  total_gap_span <- 10000 - total_plot_span
  gap <- total_gap_span / (n_side + 1)

  # Create empty list to store polygons
  plots <- list()

  # Generate subplots
  for (i in 0:(n_side - 1)) {
    for (j in 0:(n_side - 1)) {
      xmin <- gap + i * (poly_dim$x_length + gap)
      ymin <- gap + j * (poly_dim$x_length + gap)
      xmax <- xmin + poly_dim$x_length
      ymax <- ymin + poly_dim$x_length
      plots[[length(plots) + 1]] <- st_polygon(list(rbind(
        c(xmin, ymin),
        c(xmin, ymax),
        c(xmax, ymax),
        c(xmax, ymin),
        c(xmin, ymin)
      )))
    }
  }

  # Convert to sf object
  subplots_sf <- st_sf(geometry = st_sfc(plots))


  # Create a subplot region using the generated polygons
  QC_plots <- make.region(
    region.name = "study area",
    shape = subplots_sf
  )

  # Create flight lines within the quadcopter plots
  assign(paste0("QC_gimbal_design_", current_number_blocks), make.design(
    region = QC_plots,
    transect.type = "line",
    design = "systematic",
    samplers = numeric(0),
    line.length = numeric(0),
    spacing = 260*2,
    design.angle = 0,
    edge.protocol = "minus",
    truncation = 260,
    coverage.grid = cover
  ))

    assign(paste0("QC_200_design_", current_number_blocks), make.design(
    region = QC_plots,
    transect.type = "line",
    design = "systematic",
    samplers = numeric(0),
    line.length = numeric(0),
    spacing = 200,
    design.angle = 0,
    edge.protocol = "minus",
    truncation = IMAGE_WIDTH/2,
    coverage.grid = cover
  ))

    assign(paste0("QC_0_design_", current_number_blocks), make.design(
    region = QC_plots,
    transect.type = "line",
    design = "systematic",
    samplers = numeric(0),
    line.length = numeric(0),
    spacing = IMAGE_WIDTH,
    design.angle = 0,
    edge.protocol = "minus",
    truncation = IMAGE_WIDTH/2,
    coverage.grid = cover
  ))
    assign(paste0("QC_10_design_", current_number_blocks), make.design(
    region = QC_plots,
    transect.type = "line",
    design = "systematic",
    samplers = numeric(0),
    line.length = numeric(0),
    spacing = IMAGE_WIDTH*(1-0.1),
    design.angle = 0,
    edge.protocol = "minus",
    truncation = IMAGE_WIDTH/2,
    coverage.grid = cover
  ))

    assign(paste0("QC_65_design_", current_number_blocks), make.design(
    region = QC_plots,
    transect.type = "line",
    design = "systematic",
    samplers = numeric(0),
    line.length = numeric(0),
    spacing = IMAGE_WIDTH*(1-0.65),
    design.angle = 0,
    edge.protocol = "minus",
    truncation = IMAGE_WIDTH/2,
    coverage.grid = cover
  ))

    # Generate transects for the current design
    assign(paste0("QC_gimbal_transects_", current_number_blocks), generate.transects(get(paste0("QC_gimbal_design_", current_number_blocks))))
    assign(paste0("QC_200_transects_", current_number_blocks), generate.transects(get(paste0("QC_200_design_", current_number_blocks))))
    assign(paste0("QC_0_transects_", current_number_blocks), generate.transects(get(paste0("QC_0_design_", current_number_blocks))))
    assign(paste0("QC_10_transects_", current_number_blocks), generate.transects(get(paste0("QC_10_design_", current_number_blocks))))
    assign(paste0("QC_65_transects_", current_number_blocks), generate.transects(get(paste0("QC_65_design_", current_number_blocks))))

    # Run coverage simulation for the current design
    assign(paste0("QC_gimbal_design_", current_number_blocks), run.coverage(get(paste0("QC_gimbal_design_", current_number_blocks)), reps = COV_REPS))
    assign(paste0("QC_200_design_", current_number_blocks), run.coverage(get(paste0("QC_200_design_", current_number_blocks)), reps = COV_REPS))
    assign(paste0("QC_0_design_", current_number_blocks), run.coverage(get(paste0("QC_0_design_", current_number_blocks)), reps = COV_REPS))
    assign(paste0("QC_10_design_", current_number_blocks), run.coverage(get(paste0("QC_10_design_", current_number_blocks)), reps = COV_REPS))
    assign(paste0("QC_65_design_", current_number_blocks), run.coverage(get(paste0("QC_65_design_", current_number_blocks)), reps = COV_REPS))

}


# hist(get.coverage(H_SG_design))
# Plot desings
par(mfrow = c(2, 3))
plot(region, QC_0_transects_4, lwd = 0.5, col = 4)
plot(region, QC_0_transects_9, lwd = 0.5, col = 4)
plot(region, QC_0_transects_16, lwd = 0.5, col = 4)
plot(region, QC_0_transects_25, lwd = 0.5, col = 4)
par(mfrow = c(1, 1))
par(mfrow = c(2, 3))
plot(QC_0_design_4)
plot(QC_0_design_9)
plot(QC_0_design_16)
plot(QC_0_design_25)
par(mfrow = c(1, 1))

## design stats
## For details see: https://examples.distancesampling.org/dssd-getting-started/GettingStarted-distill.html#appendix-trackline-and-cyclic-trackline-lengths
## Line lenght = on effort line length
## The trackline length is the sum of the lengths of the transects plus the off-effort transit time required to complete the survey from the beginning of the first transect to the end of the last transect. The off-effort transit distance is calculated as the crow flies and may be longer in reality if transit is required around lakes, islands or coastlines etc.
## The cyclic trackline length is the trackline length plus the off-effort transit distance required to return from the end of the last transect to the beginning of the first transect.
# Extract key metrics from each simulation summary
QC_gimbal_design_4_metric <- extract_design_metrics(QC_gimbal_design_4)
QC_gimbal_design_9_metric <- extract_design_metrics(QC_gimbal_design_9)
QC_gimbal_design_16_metric <- extract_design_metrics(QC_gimbal_design_16)
QC_gimbal_design_25_metric <- extract_design_metrics(QC_gimbal_design_25)
QC_200_design_4_metric <- extract_design_metrics(QC_200_design_4)
QC_200_design_9_metric <- extract_design_metrics(QC_200_design_9)
QC_200_design_16_metric <- extract_design_metrics(QC_200_design_16)
QC_200_design_25_metric <- extract_design_metrics(QC_200_design_25)
QC_0_design_4_metric <- extract_design_metrics(QC_0_design_4)
QC_0_design_9_metric <- extract_design_metrics(QC_0_design_9)
QC_0_design_16_metric <- extract_design_metrics(QC_0_design_16)
QC_0_design_25_metric <- extract_design_metrics(QC_0_design_25)
QC_10_design_4_metric <- extract_design_metrics(QC_10_design_4)
QC_10_design_9_metric <- extract_design_metrics(QC_10_design_9)
QC_10_design_16_metric <- extract_design_metrics(QC_10_design_16)
QC_10_design_25_metric <- extract_design_metrics(QC_10_design_25)
QC_65_design_4_metric <- extract_design_metrics(QC_65_design_4)
QC_65_design_9_metric <- extract_design_metrics(QC_65_design_9)
QC_65_design_16_metric <- extract_design_metrics(QC_65_design_16)
QC_65_design_25_metric <- extract_design_metrics(QC_65_design_25)

# Combine metrics into a single dataframe
design_comparison_df <- data.frame(
  Simulation = c("gimbal4","gimbal9","gimbal16","gimba25","2004","2009","20016","20025","04","09","016","025","104","109","1016","1025","654","659","6516","6525"),
  Design = c(
    QC_gimbal_design_4_metric$design_type[1],
    QC_gimbal_design_9_metric$design_type[1],
    QC_gimbal_design_16_metric$design_type[1],
    QC_gimbal_design_25_metric$design_type[1],
    QC_200_design_4_metric$design_type[1],
    QC_200_design_9_metric$design_type[1],
    QC_200_design_16_metric$design_type[1],
    QC_200_design_25_metric$design_type[1],
    QC_0_design_4_metric$design_type[1],
    QC_0_design_9_metric$design_type[1],
    QC_0_design_16_metric$design_type[1],
    QC_0_design_25_metric$design_type[1],
    QC_10_design_4_metric$design_type[1],
    QC_10_design_9_metric$design_type[1],
    QC_10_design_16_metric$design_type[1],
    QC_10_design_25_metric$design_type[1],
    QC_65_design_4_metric$design_type[1],
    QC_65_design_9_metric$design_type[1],
    QC_65_design_16_metric$design_type[1],
    QC_65_design_25_metric$design_type[1]
  ),
  Mean_Sampler_Count = c(
    QC_gimbal_design_4_metric$mean_sampler_count,
    QC_gimbal_design_9_metric$mean_sampler_count,
    QC_gimbal_design_16_metric$mean_sampler_count,
    QC_gimbal_design_25_metric$mean_sampler_count,
    QC_200_design_4_metric$mean_sampler_count,
    QC_200_design_9_metric$mean_sampler_count,
    QC_200_design_16_metric$mean_sampler_count,
    QC_200_design_25_metric$mean_sampler_count,
    QC_0_design_4_metric$mean_sampler_count,
    QC_0_design_9_metric$mean_sampler_count,
    QC_0_design_16_metric$mean_sampler_count,
    QC_0_design_25_metric$mean_sampler_count,
    QC_10_design_4_metric$mean_sampler_count,
    QC_10_design_9_metric$mean_sampler_count,
    QC_10_design_16_metric$mean_sampler_count,
    QC_10_design_25_metric$mean_sampler_count,
    QC_65_design_4_metric$mean_sampler_count,
    QC_65_design_9_metric$mean_sampler_count,
    QC_65_design_16_metric$mean_sampler_count,
    QC_65_design_25_metric$mean_sampler_count
  ),
  Mean_Cover_Area = c(
    QC_gimbal_design_4_metric$mean_cover_area,
    QC_gimbal_design_9_metric$mean_cover_area,
    QC_gimbal_design_16_metric$mean_cover_area,
    QC_gimbal_design_25_metric$mean_cover_area,
    QC_200_design_4_metric$mean_cover_area,
    QC_200_design_9_metric$mean_cover_area,
    QC_200_design_16_metric$mean_cover_area,
    QC_200_design_25_metric$mean_cover_area,
    QC_0_design_4_metric$mean_cover_area,
    QC_0_design_9_metric$mean_cover_area,
    QC_0_design_16_metric$mean_cover_area,
    QC_0_design_25_metric$mean_cover_area,
    QC_10_design_4_metric$mean_cover_area,
    QC_10_design_9_metric$mean_cover_area,
    QC_10_design_16_metric$mean_cover_area,
    QC_10_design_25_metric$mean_cover_area,
    QC_65_design_4_metric$mean_cover_area,
    QC_65_design_9_metric$mean_cover_area,
    QC_65_design_16_metric$mean_cover_area,
    QC_65_design_25_metric$mean_cover_area
  ),
  Mean_Cover_Percentage = c(
    QC_gimbal_design_4_metric$mean_cover_percentage,
    QC_gimbal_design_9_metric$mean_cover_percentage,
    QC_gimbal_design_16_metric$mean_cover_percentage,
    QC_gimbal_design_25_metric$mean_cover_percentage,
    QC_200_design_4_metric$mean_cover_percentage,
    QC_200_design_9_metric$mean_cover_percentage,
    QC_200_design_16_metric$mean_cover_percentage,
    QC_200_design_25_metric$mean_cover_percentage,
    QC_0_design_4_metric$mean_cover_percentage,
    QC_0_design_9_metric$mean_cover_percentage,
    QC_0_design_16_metric$mean_cover_percentage,
    QC_0_design_25_metric$mean_cover_percentage,
    QC_10_design_4_metric$mean_cover_percentage,
    QC_10_design_9_metric$mean_cover_percentage,
    QC_10_design_16_metric$mean_cover_percentage,
    QC_10_design_25_metric$mean_cover_percentage,
    QC_65_design_4_metric$mean_cover_percentage,
    QC_65_design_9_metric$mean_cover_percentage,
    QC_65_design_16_metric$mean_cover_percentage,
    QC_65_design_25_metric$mean_cover_percentage
  ),
  Mean_Line_Length = c(
    QC_gimbal_design_4_metric$mean_line_length,
    QC_gimbal_design_9_metric$mean_line_length,
    QC_gimbal_design_16_metric$mean_line_length,
    QC_gimbal_design_25_metric$mean_line_length,
    QC_200_design_4_metric$mean_line_length,
    QC_200_design_9_metric$mean_line_length,
    QC_200_design_16_metric$mean_line_length,
    QC_200_design_25_metric$mean_line_length,
    QC_0_design_4_metric$mean_line_length,
    QC_0_design_9_metric$mean_line_length,
    QC_0_design_16_metric$mean_line_length,
    QC_0_design_25_metric$mean_line_length,
    QC_10_design_4_metric$mean_line_length,
    QC_10_design_9_metric$mean_line_length,
    QC_10_design_16_metric$mean_line_length,
    QC_10_design_25_metric$mean_line_length,
    QC_65_design_4_metric$mean_line_length,
    QC_65_design_9_metric$mean_line_length,
    QC_65_design_16_metric$mean_line_length,
    QC_65_design_25_metric$mean_line_length
  ),
  Mean_Trackline_Length = c(
    QC_gimbal_design_4_metric$mean_trackline,
    QC_gimbal_design_9_metric$mean_trackline,
    QC_gimbal_design_16_metric$mean_trackline,
    QC_gimbal_design_25_metric$mean_trackline,
    QC_200_design_4_metric$mean_trackline,
    QC_200_design_9_metric$mean_trackline,
    QC_200_design_16_metric$mean_trackline,
    QC_200_design_25_metric$mean_trackline,
    QC_0_design_4_metric$mean_trackline,
    QC_0_design_9_metric$mean_trackline,
    QC_0_design_16_metric$mean_trackline,
    QC_0_design_25_metric$mean_trackline,
    QC_10_design_4_metric$mean_trackline,
    QC_10_design_9_metric$mean_trackline,
    QC_10_design_16_metric$mean_trackline,
    QC_10_design_25_metric$mean_trackline,
    QC_65_design_4_metric$mean_trackline,
    QC_65_design_9_metric$mean_trackline,
    QC_65_design_16_metric$mean_trackline,
    QC_65_design_25_metric$mean_trackline
  ),
  Mean_Cyclic_Trackline_Length = c(
    QC_gimbal_design_4_metric$mean_cyclic_trackline,
    QC_gimbal_design_9_metric$mean_cyclic_trackline,
    QC_gimbal_design_16_metric$mean_cyclic_trackline,
    QC_gimbal_design_25_metric$mean_cyclic_trackline,
    QC_200_design_4_metric$mean_cyclic_trackline,
    QC_200_design_9_metric$mean_cyclic_trackline,
    QC_200_design_16_metric$mean_cyclic_trackline,
    QC_200_design_25_metric$mean_cyclic_trackline,
    QC_0_design_4_metric$mean_cyclic_trackline,
    QC_0_design_9_metric$mean_cyclic_trackline,
    QC_0_design_16_metric$mean_cyclic_trackline,
    QC_0_design_25_metric$mean_cyclic_trackline,
    QC_10_design_4_metric$mean_cyclic_trackline,
    QC_10_design_9_metric$mean_cyclic_trackline,
    QC_10_design_16_metric$mean_cyclic_trackline,
    QC_10_design_25_metric$mean_cyclic_trackline,
    QC_65_design_4_metric$mean_cyclic_trackline,
    QC_65_design_9_metric$mean_cyclic_trackline,
    QC_65_design_16_metric$mean_cyclic_trackline,
    QC_65_design_25_metric$mean_cyclic_trackline
  ),
  Mean_On_Effort = c(
    QC_gimbal_design_4_metric$mean_on_effort,
    QC_gimbal_design_9_metric$mean_on_effort,
    QC_gimbal_design_16_metric$mean_on_effort,
    QC_gimbal_design_25_metric$mean_on_effort,
    QC_200_design_4_metric$mean_on_effort,
    QC_200_design_9_metric$mean_on_effort,
    QC_200_design_16_metric$mean_on_effort,
    QC_200_design_25_metric$mean_on_effort,
    QC_0_design_4_metric$mean_on_effort,
    QC_0_design_9_metric$mean_on_effort,
    QC_0_design_16_metric$mean_on_effort,
    QC_0_design_25_metric$mean_on_effort,
    QC_10_design_4_metric$mean_on_effort,
    QC_10_design_9_metric$mean_on_effort,
    QC_10_design_16_metric$mean_on_effort,
    QC_10_design_25_metric$mean_on_effort,
    QC_65_design_4_metric$mean_on_effort,
    QC_65_design_9_metric$mean_on_effort,
    QC_65_design_16_metric$mean_on_effort,
    QC_65_design_25_metric$mean_on_effort
  ),
  Mean_Off_Effort = c(
    QC_gimbal_design_4_metric$mean_off_effort,
    QC_gimbal_design_9_metric$mean_off_effort,
    QC_gimbal_design_16_metric$mean_off_effort,
    QC_gimbal_design_25_metric$mean_off_effort,
    QC_200_design_4_metric$mean_off_effort,
    QC_200_design_9_metric$mean_off_effort,
    QC_200_design_16_metric$mean_off_effort,
    QC_200_design_25_metric$mean_off_effort,
    QC_0_design_4_metric$mean_off_effort,
    QC_0_design_9_metric$mean_off_effort,
    QC_0_design_16_metric$mean_off_effort,
    QC_0_design_25_metric$mean_off_effort,
    QC_10_design_4_metric$mean_off_effort,
    QC_10_design_9_metric$mean_off_effort,
    QC_10_design_16_metric$mean_off_effort,
    QC_10_design_25_metric$mean_off_effort,
    QC_65_design_4_metric$mean_off_effort,
    QC_65_design_9_metric$mean_off_effort,
    QC_65_design_16_metric$mean_off_effort,
    QC_65_design_25_metric$mean_off_effort
  ),
  Mean_Return_to_Home = c(
    QC_gimbal_design_4_metric$mean_return2home,
    QC_gimbal_design_9_metric$mean_return2home,
    QC_gimbal_design_16_metric$mean_return2home,
    QC_gimbal_design_25_metric$mean_return2home,
    QC_200_design_4_metric$mean_return2home,
    QC_200_design_9_metric$mean_return2home,
    QC_200_design_16_metric$mean_return2home,
    QC_200_design_25_metric$mean_return2home,
    QC_0_design_4_metric$mean_return2home,
    QC_0_design_9_metric$mean_return2home,
    QC_0_design_16_metric$mean_return2home,
    QC_0_design_25_metric$mean_return2home,
    QC_10_design_4_metric$mean_return2home,
    QC_10_design_9_metric$mean_return2home,
    QC_10_design_16_metric$mean_return2home,
    QC_10_design_25_metric$mean_return2home,
    QC_65_design_4_metric$mean_return2home,
    QC_65_design_9_metric$mean_return2home,
    QC_65_design_16_metric$mean_return2home,
    QC_65_design_25_metric$mean_return2home
  ),
  Mean_Off_Effort_Return = c(
    QC_gimbal_design_4_metric$mean_off_effort_return,
    QC_gimbal_design_9_metric$mean_off_effort_return,
    QC_gimbal_design_16_metric$mean_off_effort_return,
    QC_gimbal_design_25_metric$mean_off_effort_return,
    QC_200_design_4_metric$mean_off_effort_return,
    QC_200_design_9_metric$mean_off_effort_return,
    QC_200_design_16_metric$mean_off_effort_return,
    QC_200_design_25_metric$mean_off_effort_return,
    QC_0_design_4_metric$mean_off_effort_return,
    QC_0_design_9_metric$mean_off_effort_return,
    QC_0_design_16_metric$mean_off_effort_return,
    QC_0_design_25_metric$mean_off_effort_return,
    QC_10_design_4_metric$mean_off_effort_return,
    QC_10_design_9_metric$mean_off_effort_return,
    QC_10_design_16_metric$mean_off_effort_return,
    QC_10_design_25_metric$mean_off_effort_return,
    QC_65_design_4_metric$mean_off_effort_return,
    QC_65_design_9_metric$mean_off_effort_return,
    QC_65_design_16_metric$mean_off_effort_return,
    QC_65_design_25_metric$mean_off_effort_return
  ),
  On_Effort_Percentage = c(
    QC_gimbal_design_4_metric$on_effort_percentage,
    QC_gimbal_design_9_metric$on_effort_percentage,
    QC_gimbal_design_16_metric$on_effort_percentage,
    QC_gimbal_design_25_metric$on_effort_percentage,
    QC_200_design_4_metric$on_effort_percentage,
    QC_200_design_9_metric$on_effort_percentage,
    QC_200_design_16_metric$on_effort_percentage,
    QC_200_design_25_metric$on_effort_percentage,
    QC_0_design_4_metric$on_effort_percentage,
    QC_0_design_9_metric$on_effort_percentage,
    QC_0_design_16_metric$on_effort_percentage,
    QC_0_design_25_metric$on_effort_percentage,
    QC_10_design_4_metric$on_effort_percentage,
    QC_10_design_9_metric$on_effort_percentage,
    QC_10_design_16_metric$on_effort_percentage,
    QC_10_design_25_metric$on_effort_percentage,
    QC_65_design_4_metric$on_effort_percentage,
    QC_65_design_9_metric$on_effort_percentage,
    QC_65_design_16_metric$on_effort_percentage,
    QC_65_design_25_metric$on_effort_percentage
  ),
  Off_Effort_Percentage = c(
    QC_gimbal_design_4_metric$off_effort_percentage,
    QC_gimbal_design_9_metric$off_effort_percentage,
    QC_gimbal_design_16_metric$off_effort_percentage,
    QC_gimbal_design_25_metric$off_effort_percentage,
    QC_200_design_4_metric$off_effort_percentage,
    QC_200_design_9_metric$off_effort_percentage,
    QC_200_design_16_metric$off_effort_percentage,
    QC_200_design_25_metric$off_effort_percentage,
    QC_0_design_4_metric$off_effort_percentage,
    QC_0_design_9_metric$off_effort_percentage,
    QC_0_design_16_metric$off_effort_percentage,
    QC_0_design_25_metric$off_effort_percentage,
    QC_10_design_4_metric$off_effort_percentage,
    QC_10_design_9_metric$off_effort_percentage,
    QC_10_design_16_metric$off_effort_percentage,
    QC_10_design_25_metric$off_effort_percentage,
    QC_65_design_4_metric$off_effort_percentage,
    QC_65_design_9_metric$off_effort_percentage,
    QC_65_design_16_metric$off_effort_percentage,
    QC_65_design_25_metric$off_effort_percentage
  ),
  Return_to_Home_Percentage = c(
    QC_gimbal_design_4_metric$return2home_percentage,
    QC_gimbal_design_9_metric$return2home_percentage,
    QC_gimbal_design_16_metric$return2home_percentage,
    QC_gimbal_design_25_metric$return2home_percentage,
    QC_200_design_4_metric$return2home_percentage,
    QC_200_design_9_metric$return2home_percentage,
    QC_200_design_16_metric$return2home_percentage,
    QC_200_design_25_metric$return2home_percentage,
    QC_0_design_4_metric$return2home_percentage,
    QC_0_design_9_metric$return2home_percentage,
    QC_0_design_16_metric$return2home_percentage,
    QC_0_design_25_metric$return2home_percentage,
    QC_10_design_4_metric$return2home_percentage,
    QC_10_design_9_metric$return2home_percentage,
    QC_10_design_16_metric$return2home_percentage,
    QC_10_design_25_metric$return2home_percentage,
    QC_65_design_4_metric$return2home_percentage,
    QC_65_design_9_metric$return2home_percentage,
    QC_65_design_16_metric$return2home_percentage,
    QC_65_design_25_metric$return2home_percentage
  ),
  Off_Effort_Return_Percentage = c(
    QC_gimbal_design_4_metric$off_effort_return_percentage,
    QC_gimbal_design_9_metric$off_effort_return_percentage,
    QC_gimbal_design_16_metric$off_effort_return_percentage,
    QC_gimbal_design_25_metric$off_effort_return_percentage,
    QC_200_design_4_metric$off_effort_return_percentage,
    QC_200_design_9_metric$off_effort_return_percentage,
    QC_200_design_16_metric$off_effort_return_percentage,
    QC_200_design_25_metric$off_effort_return_percentage,
    QC_0_design_4_metric$off_effort_return_percentage,
    QC_0_design_9_metric$off_effort_return_percentage,
    QC_0_design_16_metric$off_effort_return_percentage,
    QC_0_design_25_metric$off_effort_return_percentage,
    QC_10_design_4_metric$off_effort_return_percentage,
    QC_10_design_9_metric$off_effort_return_percentage,
    QC_10_design_16_metric$off_effort_return_percentage,
    QC_10_design_25_metric$off_effort_return_percentage,
    QC_65_design_4_metric$off_effort_return_percentage,
    QC_65_design_9_metric$off_effort_return_percentage,
    QC_65_design_16_metric$off_effort_return_percentage,
    QC_65_design_25_metric$off_effort_return_percentage
  ),
  Number_of_Plots = c(
    length(QC_gimbal_design_4_metric$design_type),
    length(QC_gimbal_design_9_metric$design_type),
    length(QC_gimbal_design_16_metric$design_type),
    length(QC_gimbal_design_25_metric$design_type),
    length(QC_200_design_4_metric$design_type),
    length(QC_200_design_9_metric$design_type),
    length(QC_200_design_16_metric$design_type),
    length(QC_200_design_25_metric$design_type),
    length(QC_0_design_4_metric$design_type),
    length(QC_0_design_9_metric$design_type),
    length(QC_0_design_16_metric$design_type),
    length(QC_0_design_25_metric$design_type),
    length(QC_10_design_4_metric$design_type),
    length(QC_10_design_9_metric$design_type),
    length(QC_10_design_16_metric$design_type),
    length(QC_10_design_25_metric$design_type),
    length(QC_65_design_4_metric$design_type),
    length(QC_65_design_9_metric$design_type),
    length(QC_65_design_16_metric$design_type),
    length(QC_65_design_25_metric$design_type)
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
output_path <- here("Output", "Simulation", paste0("designsQC.RData"))
save(QC_gimbal_design_4,
    QC_gimbal_design_9,
    QC_gimbal_design_16,
    QC_gimbal_design_25,
    QC_200_design_4,
    QC_200_design_9,
    QC_200_design_16,
    QC_200_design_25,
    QC_0_design_4,
    QC_0_design_9,
    QC_0_design_16,
    QC_0_design_25,
    QC_10_design_4,
    QC_10_design_9,
    QC_10_design_16,
    QC_10_design_25,
    QC_65_design_4,
    QC_65_design_9,
    QC_65_design_16,
    QC_65_design_25, file = output_path)

# save comparison_df
output_path <- here("Output", "Simulation", paste0("desingsQC-comparsiondf.csv"))
write.csv(design_comparison_df, file = output_path, row.names = FALSE)
