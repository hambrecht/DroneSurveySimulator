# Load necessary libraries
library(here)
library(dsims)
library(knitr)
library(RColorBrewer)
library(dplyr)
library(sf)
library(units)
library(geosphere)
library(parallelly)

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

# Load density data
wmu_number_list <- c("501", "503", "512", "517", "528") 
wmu_number <- wmu_number_list[4]
input_path <- here("Output", "Density", paste0("density", wmu_number, ".RData"))
load(file = input_path)
input_path <- here("Output", "Simulation", paste0("cover-WMU", wmu_number, ".RData"))
load(file = input_path)

# Define constants
ALTITUDE <- 120 # Height in meters
CAMERA_HFOV <- 25 # Horizontal FOV in degrees. Max adjustment of 25 degrees. If more than 25 degrees then a third camera or gimbal would be needed to cover 0. Proposed intervals for adjustments are 0, 10, 20, 25
CAMERA_ANGLE <- 25 # Adjustment in degrees; max 25 with two fix cameras or 35 with gimbal in forested
IMAGE_WIDTH <- calculate_image_width(ALTITUDE, CAMERA_HFOV, CAMERA_ANGLE)
print(paste0("Half swath width is: ", IMAGE_WIDTH, " m"))

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
  scale.param = rep(50, 10),
  shape.param = seq(0.1, 4, 0.4),
  truncation = 550
)
COLORS <- brewer.pal(10, "Paired")
plot(detect_hr_overview, pop_desc, col = COLORS)
legend(x = "topright", legend = seq(0.1, 4, 0.4), col = COLORS, lty = 1, cex = 0.8)

detect_heli <- make.detectability(
  key.function = "hr",
  scale.param = 50, # heli:50
  shape.param = 2, # heli:2
  truncation = 500 # heli:20
)
plot(detect_heli, pop_desc, legend = FALSE)

detect_fixW2 <- make.detectability(
  key.function = "hn",
  scale.param = 170, # heli:50
  # shape.param = 1.3, # heli:2
  truncation = 180 # heli:20
)
plot(detect_fixW2, pop_desc, legend = FALSE)

detect_fixWG <- make.detectability(
  key.function = "hn",
  scale.param = 170, # heli:50
  # shape.param = 3, # heli:2
  truncation = 260 # heli:20
)
plot(detect_fixWG, pop_desc, legend = FALSE)

# Define and visualise uniform detection function
detect_uf <- make.detectability(
  key.function = "uf",
  scale.param = 0.8, # accounting for canopy cover
  truncation = 50
)
plot(detect_uf, pop_desc)

# Set detection function for simulation



## Simulation

SIM_REPS <- 999
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
  truncation = 500
)


# Create and run the simulation
sim_heli <- make.simulation(
  reps = SIM_REPS,
  design = heli_design,
  population.description = pop_desc,
  detectability = detect_heli,
  ds.analysis = ddf_analyses
)

sim_rnd <- make.simulation(
  reps = SIM_REPS,
  design = rnd_design,
  population.description = pop_desc,
  detectability = detect_heli,
  ds.analysis = ddf_analyses
)

sim_sys <- make.simulation(
  reps = SIM_REPS,
  design = sys_design,
  population.description = pop_desc,
  detectability = detect_heli,
  ds.analysis = ddf_analyses
)
sim_zig <- make.simulation(
  reps = SIM_REPS,
  design = zigzag_design,
  population.description = pop_desc,
  detectability = detect_heli,
  ds.analysis = ddf_analyses
)
sim_zagcom <- make.simulation(
  reps = SIM_REPS,
  design = zigzagcom_design,
  population.description = pop_desc,
  detectability = detect_heli,
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
  region = fixW_sys_design@region,
  density = fixW_density,
  N = fixW_points_count$count,
  fixed.N = TRUE
)

ddf_analyses_fixW_2C <- make.ds.analysis(
  dfmodel = list(~1, ~1),
  key = c("hn", "hr"),
  criteria = "AIC",
  truncation = 180,
  group.strata = data.frame(design.id = fixW_sys_design@region@strata.name, analysis.id = rep("A", length(fixW_sys_design@region@strata.name)))
)

ddf_analyses_fixW_G <- make.ds.analysis(
  dfmodel = list(~1, ~1),
  key = c("hn", "hr"),
  criteria = "AIC",
  truncation = 260,
  group.strata = data.frame(design.id = fixW_sys_design@region@strata.name, analysis.id = rep("A", length(fixW_sys_design@region@strata.name)))
)
fixW_sys_design@truncation <- 180
fixW_zigzag_design@truncation <- 180
sim_fixW_sys_2C <- make.simulation(
  reps = SIM_REPS,
  design = fixW_sys_design,
  population.description = pop_desc_fixW,
  detectability = detect_fixW2,
  ds.analysis = ddf_analyses_fixW_2C
)

sim_fixW_zigzag_2C <- make.simulation(
  reps = SIM_REPS,
  design = fixW_zigzag_design,
  population.description = pop_desc_fixW,
  detectability = detect_fixW2,
  ds.analysis = ddf_analyses_fixW_2C
)
fixW_sys_design@truncation <- 260
fixW_zigzag_design@truncation <- 260
sim_fixW_sys_G <- make.simulation(
  reps = SIM_REPS,
  design = fixW_sys_design,
  population.description = pop_desc_fixW,
  detectability = detect_fixWG,
  ds.analysis = ddf_analyses_fixW_G
)

sim_fixW_zigzag_G <- make.simulation(
  reps = SIM_REPS,
  design = fixW_zigzag_design,
  population.description = pop_desc_fixW,
  detectability = detect_fixWG,
  ds.analysis = ddf_analyses_fixW_G
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
  truncation = 50,
  group.strata = data.frame(design.id = quadcopter_design@region@strata.name, analysis.id = rep("A", length(quadcopter_design@region@strata.name)))
)

quadcopter_design@truncation <-50
sim_quad <- make.simulation(
  reps = SIM_REPS,
  design = quadcopter_design,
  population.description = pop_desc_quadcopter,
  detectability = detect_uf,
  ds.analysis = ddf_analyses_quadcopter
)


heli_survey <- run.survey(sim_heli)
rnd_survey <- run.survey(sim_rnd)
sys_survey <- run.survey(sim_sys)
zig_survey <- run.survey(sim_zig)
zagcom_survey <- run.survey(sim_zagcom)
fixW_sys_survey_2C <- run.survey(sim_fixW_sys_2C)
fixW_zigzag_survey_2C <- run.survey(sim_fixW_zigzag_2C)
fixW_sys_survey_G <- run.survey(sim_fixW_sys_G)
fixW_zigzag_survey_G <- run.survey(sim_fixW_zigzag_G)
quad_survey <- run.survey(sim_quad)

plot(heli_survey, region)
plot(rnd_survey, region)
plot(sys_survey, region)
plot(zig_survey, region)
plot(zagcom_survey, region)
plot(fixW_sys_survey_2C, region)
plot(fixW_zigzag_survey_2C, region)
plot(fixW_sys_survey_G, region)
plot(fixW_zigzag_survey_G, region)
plot(quad_survey, region)


# Run the full simulation
sim_heli <- run.simulation(simulation = sim_heli, run.parallel = T, max.cores=20)
sim_rnd <- run.simulation(simulation = sim_rnd, run.parallel = T, max.cores=20)
sim_sys <- run.simulation(simulation = sim_sys, run.parallel = T, max.cores=20)
sim_zig <- run.simulation(simulation = sim_zig, run.parallel = T, max.cores=20)
sim_zagcom <- run.simulation(simulation = sim_zagcom, run.parallel = T, max.cores=20)
sim_fixW_sys_2C <- run.simulation(simulation = sim_fixW_sys_2C, run.parallel = T, max.cores=20)
sim_fixW_zigzag_2C <- run.simulation(simulation = sim_fixW_zigzag_2C, run.parallel = T, max.cores=20)
sim_fixW_sys_G <- run.simulation(simulation = sim_fixW_sys_G, run.parallel = T, max.cores=20)
sim_fixW_zigzag_G <- run.simulation(simulation = sim_fixW_zigzag_G, run.parallel = T, max.cores=20)
sim_quad <- run.simulation(simulation = sim_quad, run.parallel = T, max.cores=20)

# Save simulation data
output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number, ".RData"))
# output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number,"-T",IMAGE_WIDTH,"heli-DF", detectF@key.function, ".RData"))
save(sim_heli, sim_rnd, sim_sys, sim_zig, sim_zagcom, sim_fixW_sys_2C, sim_fixW_zigzag_2C, sim_fixW_sys_G, sim_fixW_zigzag_G, sim_quad, file = output_path)



# Display results
summary(sim_heli, description.summary = FALSE)
summary(sim_rnd, description.summary = FALSE)
summary(sim_sys, description.summary = FALSE)
summary(sim_zig, description.summary = FALSE)
summary(sim_zagcom, description.summary = FALSE)
summary(sim_fixW_sys_2C, description.summary = FALSE)
summary(sim_fixW_zigzag_2C, description.summary = FALSE)
summary(sim_fixW_sys_G, description.summary = FALSE)
summary(sim_fixW_zigzag_G, description.summary = FALSE)
summary(sim_quad, description.summary = FALSE)
total_abundance
histogram.N.ests(sim_fixW_sys, use.max.reps = TRUE)
par(mfrow = c(2, 4))
histogram.N.ests(sim_heli, xlim = c(7100, 9600))
histogram.N.ests(sim_sys, xlim = c(7100, 9600))
histogram.N.ests(sim_rnd, xlim = c(7100, 9600))
histogram.N.ests(sim_zig, xlim = c(7100, 9600))
histogram.N.ests(sim_zagcom, xlim = c(7100, 9600))
histogram.N.ests(sim_fixW_sys_2C, xlim = c(2800, 3700))
histogram.N.ests(sim_fixW_zigzag_2C, xlim = c(2800, 3700))
histogram.N.ests(sim_fixW_sys_G, xlim = c(2800, 3700))
histogram.N.ests(sim_fixW_zigzag_G, xlim = c(2800, 3700))
histogram.N.ests(sim_quad, xlim = c(1000, 1410))
par(mfrow = c(1, 1))

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
