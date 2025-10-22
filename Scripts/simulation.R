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
wmu_number <- wmu_number_list[1]
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
points_within_FW_polygons <- st_join(points_sf, FW_Sys_G_design@region@region, join = st_within)

# Count the number of points in each polygon
FW_points_count <- points_within_FW_polygons %>%
  group_by(ID) %>%  # Replace `id` with the actual column name identifying polygons
  summarise(count = n()) %>%
  filter(!is.na(ID))  # Remove rows with NA in the id column

# View the result
print(FW_points_count)
FW_points_count$count

FW_density <- density
FW_density@density.surface[[1]] <- st_intersection(density@density.surface[[1]], FW_Sys_G_design@region@region)
FW_density@region.name <- FW_Sys_G_design@region@region.name
FW_density@strata.name <- FW_Sys_G_design@region@strata.name
FW_density@density.surface[[1]] <- FW_density@density.surface[[1]] %>%
  mutate(strata = ID) %>%
  select(-ID)

pop_desc_FW <- make.population.description(
  region = FW_Sys_G_design@region,
  density = FW_density,
  N = FW_points_count$count,
  fixed.N = TRUE
)

ddf_analyses_FW_2C <- make.ds.analysis(
  dfmodel = ~1,
  key = "hr",
  criteria = "AIC",
  truncation = 180,
  group.strata = data.frame(design.id = FW_Sys_G_design@region@strata.name, analysis.id = rep("A", length(FW_Sys_G_design@region@strata.name)))
)

ddf_analyses_FW_G <- make.ds.analysis(
  dfmodel = ~1,
  key = "hr",
  criteria = "AIC",
  truncation = 260,
  group.strata = data.frame(design.id = FW_Sys_G_design@region@strata.name, analysis.id = rep("A", length(FW_Sys_G_design@region@strata.name)))
)
FW_Sys_G_design@truncation <- 180
FW_ZZ_G_design@truncation <- 180
FW_Sys_2C_sim <- make.simulation(
  reps = SIM_REPS,
  design = FW_Sys_G_design,
  population.description = pop_desc_FW,
  detectability = detect_FW2,
  ds.analysis = ddf_analyses_FW_2C
)

FW_ZZ_2C_sim <- make.simulation(
  reps = SIM_REPS,
  design = FW_ZZ_G_design,
  population.description = pop_desc_FW,
  detectability = detect_FW2,
  ds.analysis = ddf_analyses_FW_2C
)
FW_Sys_G_design@truncation <- 260
FW_ZZ_G_design@truncation <- 260
FW_Sys_G_sim <- make.simulation(
  reps = SIM_REPS,
  design = FW_Sys_G_design,
  population.description = pop_desc_FW,
  detectability = detect_FWG,
  ds.analysis = ddf_analyses_FW_G
)

FW_ZZ_G_sim <- make.simulation(
  reps = SIM_REPS,
  design = FW_ZZ_G_design,
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
plot(QC_Sys_survey, region, main = "hallo!!!!!!!!!!!!!")
print(plot)


#create gif
# iterate_count <- 99
# for (i in 1:iterate_count) {
#   print(i)
#   QC_Sys_survey <- run.survey(QC_Sys_sim)
    
#   # Create a filename for saving the plot
#   filename <- here("Output", "Plots", "sim_gif", paste0("plot_", i, ".jpg"))
  
#   # open png object
#   png(filename, width = 1644, height = 682)
#   # Plot the results
#   plot(QC_Sys_survey, region)
  
  

#   # Save the plot as a JPEG file
#   dev.off()
# }


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


# input_path <- here("Output", "Simulation", paste0("cover-WMU", wmu_number, "-varEffort.RData")) # nothing can be removed
load(output_path)
ls()
# List all objects containing 'QC_Sys_design_'
loaded_objects <- ls(pattern = "QC_Sys_design_")
N_factor <- c(0.25, 0.5, 0.75, 1)

for (design_name in loaded_objects[5]) {
  design <- get(design_name)

  # Perform spatial join to count points within each polygon
  points_within_QC_polygons <- st_join(points_sf, design@region@region, join = st_within)

  # Count the number of points in each polygon
  QC_points_count <- points_within_QC_polygons %>%
    group_by(ID) %>%
    summarise(count = n()) %>%
    filter(!is.na(ID))

  id_df <- data.frame(ID = design@region@strata.name, count = 0)
  merged_df <- full_join(id_df, QC_points_count, by = "ID")
  merged_df$count <- ifelse(is.na(merged_df$count.y), merged_df$count.x, merged_df$count.y)
  QC_pop_count <- merged_df %>% select(ID, count)

  QC_Sys_density <- density
  QC_Sys_density@density.surface[[1]] <- st_intersection(density@density.surface[[1]], design@region@region)
  QC_Sys_density@region.name <- design@region@region.name
  QC_Sys_density@strata.name <- design@region@strata.name
  QC_Sys_density@density.surface[[1]] <- QC_Sys_density@density.surface[[1]] %>%
    mutate(strata = ID) %>%
    select(-ID)

  pop_desc_QC <- make.population.description(
    region = design@region,
    density = QC_Sys_density,
    N = QC_pop_count$count,
    fixed.N = TRUE
  )

  ddf_analyses_QC <- make.ds.analysis(
    dfmodel = ~1,
    key = "hr",
    criteria = "AIC",
    truncation = 260,
    group.strata = data.frame(design.id = design@region@strata.name, analysis.id = rep("A", length(design@region@strata.name)))
  )

  for (FACTOR in seq_along(N_factor)) {
    print(FACTOR)
    pop_desc_QC_density <- pop_desc_QC
    pop_desc_QC_density@N <- pop_desc_QC@N * N_factor[FACTOR]

    QC_Sys_sim_density <- make.simulation(
      reps = SIM_REPS,
      design = design,
      population.description = pop_desc_QC_density,
      detectability = detect_FWG,
      ds.analysis = ddf_analyses_QC
    )

    QC_Sys_sim_density <- run.simulation(simulation = QC_Sys_sim_density, run.parallel = TRUE, max.cores = 20)

    output_path <- here("Output", "Simulation", paste0(design_name, "-density_sim-WMU", wmu_number, "-D", N_factor[FACTOR], ".RData"))
    save(QC_Sys_sim_density, file = output_path)
  }
}

print('Done')