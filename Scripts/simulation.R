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

# Load density data
wmu_number_list <- c("501", "503", "512", "517", "528")
wmu_number <- wmu_number_list[2]
input_path <- here("Output", "Density", paste0("density", wmu_number, ".RData"))
load(file = input_path)
input_path <- here("Output", "Simulation", paste0("cover-WMU", wmu_number, ".RData"))
load(file = input_path)

# Save simulation data
input_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number, ".RData"))
load(file = input_path)




# Define constants
ALTITUDE <- 120 # Height in meters
CAMERA_HFOV <- 25 # Horizontal FOV in degrees. Max adjustment of 25 degrees. If more than 25 degrees then a third camera or gimbal would be needed to cover 0. Proposed intervals for adjustments are 0, 10, 20, 25
CAMERA_ANGLE <- 35 # Adjustment in degrees; max 25 with two fix cameras or 35 with gimbal in forested
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

detect_H <- make.detectability(
  key.function = "hn",
  scale.param = 140, # heli:50
  # shape.param = 2, # heli:2
  truncation = 600 # heli:20
)
plot(detect_H, pop_desc, legend = FALSE)

detect_FW2 <- make.detectability(
  key.function = "hn",
  scale.param = 170,
  # shape.param = 1.3,
  truncation = 180 
)
plot(detect_FW2, pop_desc, legend = FALSE)

detect_FWG <- make.detectability(
  key.function = "hn",
  scale.param = 170,
  # shape.param = 3,
  truncation = 260
)
plot(detect_FWG, pop_desc, legend = FALSE)

# Define and visualise uniform detection function
detect_uf <- make.detectability(
  key.function = "uf",
  scale.param = 1, # accounting for canopy cover
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
  dfmodel = ~1,
  key = "hr",
  criteria = "AIC",
  truncation = 500
)


# Create and run the simulation
H_SG_sim <- make.simulation(
  reps = SIM_REPS,
  design = H_SG_design,
  population.description = pop_desc,
  detectability = detect_H,
  ds.analysis = ddf_analyses
)

# Drone sims
example_population <- generate.population(object = pop_desc, detectability = detect_uf, region = region)
# plot(example_population, region)
# termine abundance in each subplot
# Convert points dataframe to sf object
points_sf <- st_as_sf(example_population@population, coords = c("x", "y"), crs = st_crs(wmu))

# Perform spatial join to count points within each polygon
points_within_FW_polygons <- st_join(points_sf, FW_Sys_design@region@region, join = st_within)

# Count the number of points in each polygon
FW_points_count <- points_within_FW_polygons %>%
  group_by(ID) %>%  # Replace `id` with the actual column name identifying polygons
  summarise(count = n()) %>%
  filter(!is.na(ID))  # Remove rows with NA in the id column

# View the result
print(FW_points_count)
FW_points_count$count

FW_density <- density
FW_density@density.surface[[1]] <- st_intersection(density@density.surface[[1]], FW_Sys_design@region@region)
FW_density@region.name <- FW_Sys_design@region@region.name
FW_density@strata.name <- FW_Sys_design@region@strata.name
FW_density@density.surface[[1]] <- FW_density@density.surface[[1]] %>%
  mutate(strata = ID) %>%
  select(-ID)

pop_desc_FW <- make.population.description(
  region = FW_Sys_design@region,
  density = FW_density,
  N = FW_points_count$count,
  fixed.N = TRUE
)

ddf_analyses_FW_2C <- make.ds.analysis(
  dfmodel = ~1,
  key = "hr",
  criteria = "AIC",
  truncation = 180,
  group.strata = data.frame(design.id = FW_Sys_design@region@strata.name, analysis.id = rep("A", length(FW_Sys_design@region@strata.name)))
)

ddf_analyses_FW_G <- make.ds.analysis(
  dfmodel = ~1,
  key = "hr",
  criteria = "AIC",
  truncation = 260,
  group.strata = data.frame(design.id = FW_Sys_design@region@strata.name, analysis.id = rep("A", length(FW_Sys_design@region@strata.name)))
)
FW_Sys_design@truncation <- 180
FW_ZZ_design@truncation <- 180
FW_Sys_2C_sim <- make.simulation(
  reps = SIM_REPS,
  design = FW_Sys_design,
  population.description = pop_desc_FW,
  detectability = detect_FW2,
  ds.analysis = ddf_analyses_FW_2C
)

FW_ZZ_2C_sim <- make.simulation(
  reps = SIM_REPS,
  design = FW_ZZ_design,
  population.description = pop_desc_FW,
  detectability = detect_FW2,
  ds.analysis = ddf_analyses_FW_2C
)
FW_Sys_design@truncation <- 260
FW_ZZ_design@truncation <- 260
FW_Sys_G_sim <- make.simulation(
  reps = SIM_REPS,
  design = FW_Sys_design,
  population.description = pop_desc_FW,
  detectability = detect_FWG,
  ds.analysis = ddf_analyses_FW_G
)

FW_ZZ_G_sim <- make.simulation(
  reps = SIM_REPS,
  design = FW_ZZ_design,
  population.description = pop_desc_FW,
  detectability = detect_FWG,
  ds.analysis = ddf_analyses_FW_G
)
# termine abundance in each subplot

# # Perform spatial join to count points within each polygon
# points_within_QC_polygons <- st_join(points_sf, QC_Sys_nadir_design@region@region, join = st_within)

# QC_points_count <- QC_Sys_nadir_design@region@strata.name
# # Count the number of points in each polygon
# QC_points_count <- points_within_QC_polygons %>%
#   group_by(ID) %>%  # Replace `id` with the actual column name identifying polygons
#   summarise(count = n()) %>%
#   filter(!is.na(ID))  # Remove rows with NA in the id column

# # View the result
# print(QC_points_count)
# QC_points_count$ID

# id_df <- data.frame(
#   ID = QC_Sys_nadir_design@region@strata.name,
#   count = 0
# )
# merged_df <- full_join(id_df, QC_points_count, by = "ID")
# # Replace NA values in the count column with 0
# merged_df$count <- ifelse(is.na(merged_df$count.y), merged_df$count.x, merged_df$count.y)
# QC_pop_count <- merged_df %>% select(ID, count)

# QC_Sys_density <- density
# QC_Sys_density@density.surface[[1]] <- st_intersection(density@density.surface[[1]], QC_Sys_nadir_design@region@region)
# QC_Sys_density@region.name <- QC_Sys_nadir_design@region@region.name
# QC_Sys_density@strata.name <- QC_Sys_nadir_design@region@strata.name
# QC_Sys_density@density.surface[[1]] <- QC_Sys_density@density.surface[[1]] %>%
#   mutate(strata = ID) %>%
#   select(-ID)




# Perform spatial join to count points within each polygon
points_within_QC_polygons <- st_join(points_sf, QC_Sys_design@region@region, join = st_within)

QC_points_count <- QC_Sys_design@region@strata.name
# Count the number of points in each polygon
QC_points_count <- points_within_QC_polygons %>%
  group_by(ID) %>%  # Replace `id` with the actual column name identifying polygons
  summarise(count = n()) %>%
  filter(!is.na(ID))  # Remove rows with NA in the id column

# View the result
print(QC_points_count)
QC_points_count$ID

id_df <- data.frame(
  ID = QC_Sys_design@region@strata.name,
  count = 0
)
merged_df <- full_join(id_df, QC_points_count, by = "ID")
# Replace NA values in the count column with 0
merged_df$count <- ifelse(is.na(merged_df$count.y), merged_df$count.x, merged_df$count.y)
QC_pop_count <- merged_df %>% select(ID, count)

QC_Sys_density <- density
QC_Sys_density@density.surface[[1]] <- st_intersection(density@density.surface[[1]], QC_Sys_design@region@region)
QC_Sys_density@region.name <- QC_Sys_design@region@region.name
QC_Sys_density@strata.name <- QC_Sys_design@region@strata.name
QC_Sys_density@density.surface[[1]] <- QC_Sys_density@density.surface[[1]] %>%
  mutate(strata = ID) %>%
  select(-ID)


pop_desc_QC_nadir <- make.population.description(
  region = QC_Sys_design@region,
  density = QC_Sys_density,
  N = QC_pop_count$count,
  fixed.N = TRUE
)

ddf_analyses_QC_nadir <- make.ds.analysis(
  dfmodel = ~1,
  key = "hr",
  criteria = "AIC",
  truncation = 50,
  group.strata = data.frame(design.id = QC_Sys_design@region@strata.name, analysis.id = rep("A", length(QC_Sys_design@region@strata.name)))
)

QC_Sys_design@truncation <- 50
QC_Sys_nadir_sim <- make.simulation(
  reps = SIM_REPS,
  design = QC_Sys_design,
  population.description = pop_desc_QC_nadir,
  detectability = detect_uf,
  ds.analysis = ddf_analyses_QC_nadir
)


pop_desc_QC <- make.population.description(
  region = QC_Sys_design@region,
  density = QC_Sys_density,
  N = QC_pop_count$count,
  fixed.N = TRUE
)

ddf_analyses_QC <- make.ds.analysis(
  dfmodel = ~1,
  key = "hr",
  criteria = "AIC",
  truncation = 260,
  group.strata = data.frame(design.id = QC_Sys_design@region@strata.name, analysis.id = rep("A", length(QC_Sys_design@region@strata.name)))
)

QC_Sys_design@truncation <- 260
QC_Sys_sim <- make.simulation(
  reps = SIM_REPS,
  design = QC_Sys_design,
  population.description = pop_desc_QC,
  detectability = detect_FWG,
  ds.analysis = ddf_analyses_QC
)


H_SG_survey <- run.survey(H_SG_sim)
FW_Sys_survey_2C <- run.survey(FW_Sys_2C_sim)
FW_ZZ_survey_2C <- run.survey(FW_ZZ_2C_sim)
FW_Sys_survey_G <- run.survey(FW_Sys_G_sim)
FW_ZZ_survey_G <- run.survey(FW_ZZ_G_sim)
QC_Sys_nadir_survey <- run.survey(QC_Sys_nadir_sim)
QC_Sys_survey <- run.survey(QC_Sys_sim)

plot(H_SG_survey, region)
plot(FW_Sys_survey_2C, region)
plot(FW_ZZ_survey_2C, region)
plot(FW_Sys_survey_G, region)
plot(FW_ZZ_survey_G, region)
plot(QC_Sys_nadir_survey, region)
plot(QC_Sys_survey, region)


# Run the full simulation
H_SG_sim <- run.simulation(simulation = H_SG_sim, run.parallel = T, max.cores = 20)
FW_Sys_2C_sim <- run.simulation(simulation = FW_Sys_2C_sim, run.parallel = T, max.cores = 20)
FW_ZZ_2C_sim <- run.simulation(simulation = FW_ZZ_2C_sim, run.parallel = T, max.cores = 20)
FW_Sys_G_sim <- run.simulation(simulation = FW_Sys_G_sim, run.parallel = T, max.cores = 20)
FW_ZZ_G_sim <- run.simulation(simulation = FW_ZZ_G_sim, run.parallel = T, max.cores = 20)
QC_Sys_nadir_sim <- run.simulation(simulation = QC_Sys_nadir_sim, run.parallel = T, max.cores = 20)
QC_Sys_sim <- run.simulation(simulation = QC_Sys_sim, run.parallel = T, max.cores = 20)


# Save simulation data
output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number, ".RData"))
# output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number,"-T",IMAGE_WIDTH,"H_SG-DF", detectF@key.function, ".RData"))
save(H_SG_sim, FW_Sys_2C_sim, FW_Sys_2C_sim, FW_ZZ_2C_sim, FW_Sys_G_sim, FW_ZZ_G_sim, QC_Sys_nadir_sim, QC_Sys_sim, file = output_path)



# Display results
summary(H_SG_sim, description.summary = FALSE)
summary(FW_Sys_2C_sim, description.summary = FALSE)
summary(FW_ZZ_2C_sim, description.summary = FALSE)
summary(FW_Sys_G_sim, description.summary = FALSE)
summary(FW_ZZ_G_sim, description.summary = FALSE)
summary(QC_Sys_nadir_sim, description.summary = FALSE)
summary(QC_Sys_sim, description.summary = FALSE)
total_abundance
histogram.N.ests(Rnd_sim, use.max.reps = TRUE)

# Define labels
labels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
xlims <- c(49500, 80000)
par(mfrow = c(2, 5))

# Plot each histogram and add the corresponding label
histogram.N.ests(H_SG_sim, xlim = xlims)
usr <- par("usr")
text(x = usr[2], y = usr[4], labels = labels[1], adj = c(1.2, 1.2), col = "black", cex = 1.2, font = 2)

histogram.N.ests(FW_Sys_2C_sim, xlim = c(30000, 33500))
usr <- par("usr")
text(x = usr[2], y = usr[4], labels = labels[6], adj = c(1.2, 1.2), col = "black", cex = 1.2, font = 2)

histogram.N.ests(FW_ZZ_2C_sim, xlim = c(30000, 33500))
usr <- par("usr")
text(x = usr[2], y = usr[4], labels = labels[7], adj = c(1.2, 1.2), col = "black", cex = 1.2, font = 2)

histogram.N.ests(FW_Sys_G_sim, xlim = c(30000, 33500))
usr <- par("usr")
text(x = usr[2], y = usr[4], labels = labels[8], adj = c(1.2, 1.2), col = "black", cex = 1.2, font = 2)

histogram.N.ests(FW_ZZ_G_sim, xlim = c(30000, 33500))
usr <- par("usr")
text(x = usr[2], y = usr[4], labels = labels[9], adj = c(1.2, 1.2), col = "black", cex = 1.2, font = 2)

histogram.N.ests(QC_Sys_sim, xlim = c(6500, 8900))
usr <- par("usr")
text(x = usr[2], y = usr[4], labels = labels[10], adj = c(1.2, 1.2), col = "black", cex = 1.2, font = 2)

histogram.N.ests(QC_Sys_nadir_sim, xlim = c(6500, 8900))
usr <- par("usr")
text(x = usr[2], y = usr[4], labels = labels[10], adj = c(1.2, 1.2), col = "black", cex = 1.2, font = 2)


par(mfrow = c(1, 1))



# Density
N_factor <- c(1.0, 0.75, 0.5, 0.25)
FACTOR <- 1
H_SG_sim_density <- H_SG_sim
FW_Sys_G_sim_density <- FW_Sys_G_sim
QC_Sys_sim_density <- QC_Sys_sim
# Save simulation data
output_path <- here("Output", "Simulation", paste0("density_sim-WMU", wmu_number, "-D", N_factor[FACTOR], ".RData"))
# output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number,"-T",IMAGE_WIDTH,"H_SG-DF", detectF@key.function, ".RData"))
save(H_SG_sim_density, FW_Sys_G_sim_density, QC_Sys_sim_density, file = output_path)

for (FACTOR in 2:length(N_factor)) {
  output_path <- here("Output", "Simulation", paste0("density_sim-WMU", wmu_number, "-D", N_factor[FACTOR], ".RData"))
  load(output_path)
  #pop_desc_density <- pop_desc
  #pop_desc_density@N <- pop_desc@N * N_factor[FACTOR]
  #H_SG_sim_density <- make.simulation(
  #  reps = SIM_REPS,
  #  design = H_SG_design,
  #  population.description = pop_desc_density,
  #  detectability = detect_H,
  #  ds.analysis = ddf_analyses
  #)
  
  #pop_desc_FW_density <- pop_desc_FW
  #pop_desc_FW_density@N <-pop_desc_FW@N * N_factor[FACTOR]
  #FW_Sys_G_sim_density <- make.simulation(
  #  reps = SIM_REPS,
  #  design = FW_Sys_design,
  #  population.description = pop_desc_FW_density,
  #  detectability = detect_FWG,
  #  ds.analysis = ddf_analyses_FW_G
  #)
  
  pop_desc_QC_density <- pop_desc_QC
  pop_desc_QC_density@N <- pop_desc_QC@N* N_factor[FACTOR]
  QC_Sys_sim_density <- make.simulation(
    reps = SIM_REPS,
    design = QC_Sys_design,
    population.description = pop_desc_QC_density,
    detectability = detect_FWG,
    ds.analysis = ddf_analyses_QC
  )
  
  #H_SG_sim_density <- run.simulation(simulation = H_SG_sim_density, run.parallel = T, max.cores = 20)
  #FW_Sys_G_sim_density <- run.simulation(simulation = FW_Sys_G_sim_density, run.parallel = T, max.cores = 20)
  QC_Sys_sim_density <- run.simulation(simulation = QC_Sys_sim_density, run.parallel = T, max.cores = 20)
  
  # Save simulation data
  # output_path <- here("Output", "Simulation", paste0("density_sim-WMU", wmu_number, "-D", N_factor[FACTOR], ".RData"))
  # output_path <- here("Output", "Simulation", paste0("simulation-WMU", wmu_number,"-T",IMAGE_WIDTH,"H_SG-DF", detectF@key.function, ".RData"))
  save(H_SG_sim_density, FW_Sys_G_sim_density, QC_Sys_sim_density, file = output_path)
}


# # Extract metrics for each simulation
# metrics_H_SG <- extract_metrics(H_SG_sim)
# metrics_Sys <- extract_metrics(Sys_sim)
# metrics_ZZ <- extract_metrics(ZZ_sim)
# metrics_ZZC <- extract_metrics(ZZC_sim)

# # Combine metrics into a single dataframe
# comparison_df <- data.frame(
#   Simulation = c("H_SG", "Sys", "Zig", "Zagcom"),
#   Mean_Estimate = c(metrics_H_SG$mean_estimate, metrics_Sys$mean_estimate, metrics_ZZ$mean_estimate, metrics_ZZC$mean_estimate),
#   Percent_Bias = c(metrics_H_SG$percent_bias, metrics_Sys$percent_bias, metrics_ZZ$percent_bias, metrics_ZZC$percent_bias),
#   RMSE = c(metrics_H_SG$rmse, metrics_Sys$rmse, metrics_ZZ$rmse, metrics_ZZC$rmse),
#   CI_Coverage_Prob = c(metrics_H_SG$ci_coverage_prob, metrics_Sys$ci_coverage_prob, metrics_ZZ$ci_coverage_prob, metrics_ZZC$ci_coverage_prob),
#   Mean_SE = c(metrics_H_SG$mean_se, metrics_Sys$mean_se, metrics_ZZ$mean_se, metrics_ZZC$mean_se),
#   SD_of_Means = c(metrics_H_SG$sd_of_means, metrics_Sys$sd_of_means, metrics_ZZ$sd_of_means, metrics_ZZC$sd_of_means),
#   Mean_Cover_Area = c(metrics_H_SG$mean_cover_area, metrics_Sys$mean_cover_area, metrics_ZZ$mean_cover_area, metrics_ZZC$mean_cover_area),
#   Mean_Effort = c(metrics_H_SG$mean_effort, metrics_Sys$mean_effort, metrics_ZZ$mean_effort, metrics_ZZC$mean_effort),
#   Mean_n = c(metrics_H_SG$mean_n, metrics_Sys$mean_n, metrics_ZZ$mean_n, metrics_ZZC$mean_n),
#   Mean_k = c(metrics_H_SG$mean_k, metrics_Sys$mean_k, metrics_ZZ$mean_k, metrics_ZZC$mean_k),
#   Mean_ER = c(metrics_H_SG$mean_ER, metrics_Sys$mean_ER, metrics_ZZ$mean_ER, metrics_ZZC$mean_ER),
#   Mean_se_ER = c(metrics_H_SG$mean_se_ER, metrics_Sys$mean_se_ER, metrics_ZZ$mean_se_ER, metrics_ZZC$mean_se_ER)
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
#     key = "hr",
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
