# Load necessary libraries
library(here)
library(dsims)
library(knitr)
library(RColorBrewer)
library(dplyr)
library(sf)
library(units)

# Check if pbapply is installed
if (!requireNamespace("pbapply", quietly = TRUE, dependencies = TRUE)) {
 message("The 'pbapply' package is not installed. Installing it now...")
 install.packages("pbapply")
} else {
 message("The 'pbapply' package is already installed.")
}

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


# Define constants
ALTITUDE <- 120 # Height in meters
CAMERA_HFOV <- 35 # Horizontal FOV in degrees. Max adjustment of 25 degrees. If more than 25 degrees then a third camera or gimbal would be needed to cover 0. Proposed intervals for adjustments are 0, 10, 20, 25
CAMERA_ANGLE <- 0 # Adjustment in degrees; max 25 with two fix cameras or 35 with gimbal in forested
IMAGE_WIDTH <- calculate_image_width(ALTITUDE, CAMERA_HFOV, CAMERA_ANGLE)
print(paste0("Half swath width is: ", IMAGE_WIDTH, " m"))
SIM_REPS <- 999


# Set detection function for simulation


# List all objects containing 'QC_Sys_design_'
input_path <- here("Output", "Simulation", paste0("designsQC.RData"))
load(input_path)
st_crs(region@region) <- 32633 # Set the CRS for the region
region@units <- "m"

# Define a curved path (e.g., a semi-ellipse or arc)
centres <- list(
  c(5000, 3000),
  c(5789, 3108),
  c(6551, 3421),
  c(7236, 3906),
  c(7765, 4509),
  c(8039, 5165),
  c(7968, 5803),
  c(7506, 6354),
  c(6678, 6758),
  c(5591, 6972),
  c(4408, 6972),
  c(3321, 6758),
  c(2493, 6354),
  c(2031, 5803),
  c(1960, 5165),
  c(2234, 4509),
  c(2763, 3906),
  c(3448, 3421),
  c(4210, 3108),
  c(4999, 3000)
)

density <- make.density(region = region,
                        x.space = 500,
                        constant = 1)

# Add each hotspot along the curve
for (centre in centres) {
  density <- add.hotspot(object = density,
                         centre = centre,
                         sigma = 1500,
                         amplitude = 2)
}
density <- add.hotspot(object = density,
                       centre = c(5000, 8000),
                       sigma = 1500,
                       amplitude = -4)
plot(density, region, scale = 2)

pop_desc <- make.population.description(
      region = region,
      density = density,
      N = 200, # Total population size
      fixed.N = T
    )



# Define and visualise detection function
detect_G <- make.detectability(
  key.function = "hn",
  scale.param = 170,
  # shape.param = 3,
  truncation = 260
)
plot(detect_G, pop_desc, legend = FALSE)

# Define and visualise uniform detection function
detect_NADIR <- make.detectability(
  key.function = "uf",
  scale.param = 1, # accounting for canopy cover
  truncation = IMAGE_WIDTH
)
plot(detect_NADIR, pop_desc)

ddf_analyses_G <- make.ds.analysis(
  dfmodel = ~1,
  key = "hr",
  criteria = "AIC",
  truncation = 260,
  group.strata = data.frame(design.id = QC_Sys_design@region@strata.name, analysis.id = rep("A", length(QC_Sys_design@region@strata.name)))
)
ddf_analyses_nadir <- make.ds.analysis(
  dfmodel = ~1,
  key = "hr",
  criteria = "AIC",
  truncation = IMAGE_WIDTH,
  group.strata = data.frame(design.id = QC_Sys_design@region@strata.name, analysis.id = rep("A", length(QC_Sys_design@region@strata.name)))
)

ABUNDANCE_LIST <- c(5,10,20,30,40)

loaded_objects <- ls(pattern = "^QC_")
dev.off() # clear plots from memory
for (design_name in loaded_objects) {
  print(design_name)
  design <- get(design_name)
  st_crs(design@region@region) <- 32633 # Set the CRS for the region
  design@region@units <- "m"

  # create design density
  design_density <- make.density(region = design@region,
                        x.space = 500,
                        constant = 1)

  # Add each hotspot along the curve
  for (centre in centres) {
    design_density <- add.hotspot(object = design_density,
                          centre = centre,
                          sigma = 1500,
                          amplitude = 2)
  }
  design_density <- add.hotspot(object = design_density,
                        centre = c(5000, 8000),
                        sigma = 1500,
                        amplitude = -4)
  # plot(design_density, design@region, scale = 2)
  design_density@density.surface[[1]]$density <- design_density@density.surface[[1]]$density * 0.01 # fixes memory issues

  ddf_analyses <- make.ds.analysis(
  dfmodel = ~1,
  key = "hr",
  criteria = "AIC",
  truncation = IMAGE_WIDTH,
  group.strata = data.frame(design.id = design@region@strata.name, analysis.id = rep("A", length(design@region@strata.name)))
  )


  # If design_name does contain "gimbal" then use detect_G, otherwise use detect_NADIR

  if (grepl("gimbal", design_name, ignore.case = TRUE)) {
    detect_fun <- detect_G
    ddf_analyses@truncation <- detect_G@truncation # set truncation distance suitable for gimbal
  } else {
    detect_fun <- detect_NADIR
  }

  for (ABUNDANCE in ABUNDANCE_LIST) {
    print(ABUNDANCE)

    # Create population description
    ex_pop_desc <- make.population.description(
      region = region,
      density = density,
      N = 200, # Total population size
      fixed.N = T
    )

    example_population <- generate.population(object = ex_pop_desc, detectability = detect_NADIR, region = region)
    # termine abundance in each subplot
    # Convert points dataframe to sf object
    points_sf <- st_as_sf(example_population@population, coords = c("x", "y"), crs = st_crs(region@region))
    plot(design@region@region)
    plot(points_sf$geometry, add = T, col = "red", pch = 20)


    # Perform spatial join to count points within each polygon
    points_within_design <- st_join(points_sf, design@region@region, join = st_within, crs = st_crs(region@region))
    # Filter out points that do not fall within any polygon (i.e., remove NAs)
    points_within_design <- points_within_design[!is.na(points_within_design[[ncol(points_within_design)]]), ]

    within_idx <- lengths(st_within(points_sf, design@region@region)) > 0
    points_within_design <- points_sf[within_idx, ]

    # Find which polygon each point falls within
    within_list <- st_within(points_sf, design@region@region)

    # Assign polygon index (or NA) to each point
    strata_names <- design@region@strata.name
    points_sf$strata.name <- sapply(within_list, function(x) if(length(x) > 0) strata_names[x[1]] else NA)


    # Count points per polygon
    points_count <- points_sf %>%
      filter(!is.na(strata.name)) %>%
      group_by(strata.name) %>%
      summarise(count = n()) %>%
      st_drop_geometry()

    # print(points_count)

    # Create population description
    pop_desc <- make.population.description(
      region = design@region,
      density = design_density,
      N = points_count$count, # Total population size
      fixed.N = F
    )

    QC_Sys_sim_density <- make.simulation(
      reps = SIM_REPS,
      design = design,
      population.description = pop_desc,
      detectability = detect_fun,
      ds.analysis = ddf_analyses
    )
    stop()
    # survey <- run.survey(QC_Sys_sim_density)
    QC_Sys_sim_density <- run.simulation(simulation = QC_Sys_sim_density, run.parallel = TRUE, max.cores = 20)

    output_path <- here("Output", "Simulation", paste0(design_name, "-density_sim-A", ABUNDANCE, ".RData"))
    save(QC_Sys_sim_density, file = output_path)
    # clear memory
    rm(QC_Sys_sim_density, points_sf, points_within_design, within_list, points_count, pop_desc, ex_pop_desc, example_population)
    gc()
    
  }
}

print('Done')