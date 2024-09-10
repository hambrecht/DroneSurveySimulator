# Load required packages
library(here)       # For relative paths
library(purrr)      # For functional programming
library(foreach)    # For parallel processing
library(doParallel) # For parallel processing
library(dsims)      # For density surface modelling
library(Distance)   # For distance sampling functions
library(dsm)        # For distance sampling models

# Define custom functions

#' Check for Required Columns in a Dataframe
#'
#' Checks if all specified columns are present in a dataframe. The check is case-insensitive.
#'
#' @param df A dataframe to check.
#' @param required_cols A character vector of column names that must be present in the dataframe.
#'
#' @return Throws an error if any required columns are missing.
#'
#' @examples
#' df <- data.frame(A = 1:5, b = 6:10, D = 11:15, E = 16:20)
#' required_columns <- c("A", "B", "C")
#' check_columns_present(df, required_columns)
#'
#' @seealso \link[base]{stop}
#' @export
check_columns_present <- function(df, required_cols) {
  # Convert column names to lowercase
  actual_cols <- tolower(colnames(df))
  required_cols_lower <- tolower(required_cols)
  
  # Identify missing columns
  missing_cols <- setdiff(required_cols_lower, actual_cols)
  
  # Stop execution if there are missing columns
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
}

possDS <- possibly(.f = Distance::ds, otherwise = NULL)

### Set Variables
# Set the number of cores to use for parallel processing
registerDoParallel(cores = parallel_cores)

# Define grid cell size (in meters)
GRID_SIZE <- 500
wum_number <-'528' # 501, 503, 512 517, 528

## Load and Check Data

# Load processed data
input_path <- here("Output", "PrepData", paste0("prepared_",wum_number,".RData"))
load(file = input_path)


# Define required columns for dataframes
segdata_required <- c("longitude", "latitude", "Transect.Label", "Sample.Label", "x", "y", "Effort")
distdata_required <- c("object", "size", "longitude", "latitude", "x", "y", "Effort")
obsdata_required <- c("object", "Sample.Label", "size", "distance", "Effort")

# Check presence of required columns
check_columns_present(segdata, segdata_required)
check_columns_present(distdata, distdata_required)
check_columns_present(obsdata, obsdata_required)

## Density Surface Template

# Create a density object for coordinate retrieval and abundance prediction
# Define survey region
# region <- make.region(
#   region.name = "study area",
#   shape = wmu
# )

# Define survey region with strata
region <- make.region(
  region.name = "study area with strata",
  shape = wmu_strata
)

# Create density surface
density <- dsims::make.density(
  region = region,
  x.space = GRID_SIZE
)

# Extract coordinates
coords <- sf::st_drop_geometry(density@density.surface[[1]][, c("x", "y")])


# Plot distance data
break_bins <- seq(from = 0, to = 0.6, by = 0.05)
hist(distdata$distance, main = "Moose line transects", xlab = "Distance (km)", breaks = break_bins)

max_truncation_distance <- round(max(distdata$distance),1)
truncation_distance_list <- c(max_truncation_distance, max_truncation_distance*.95, max_truncation_distance*.9, 0.5)
detfc_list <- list()
detfc_list_bin <- list()


for (truncation_distance in truncation_distance_list) {
  assign(paste0("detfc_hr_null_trunc_",truncation_distance), possDS(distdata, truncation_distance, key = "hr", adjustment = NULL))
  assign(paste0("detfc_hn_null_trunc_",truncation_distance), possDS(distdata, truncation_distance, key = "hn", adjustment = NULL))
  assign(paste0("detfc_hn_cos_trunc_",truncation_distance), possDS(distdata, truncation_distance, key = "hn", adjustment = "cos"))
  detfc_list <- list(
    paste0("detfc_hr_null_trunc_",truncation_distance),
    paste0("detfc_hn_null_trunc_",truncation_distance),
    paste0("detfc_hn_cos_trunc_",truncation_distance)
  )

  for ( i in seq(from = 0.08, to = 0.15, by = 0.01)){
    assign(paste0("detfc_hr_null_trunc_",truncation_distance,"bin_",i), possDS(distdata, truncation_distance, key = "hr", adjustment = NULL, cutpoints = seq(from = 0, to = truncation_distance, by = i)))
    assign(paste0("detfc_hn_null_trunc_",truncation_distance,"bin_",i), possDS(distdata, truncation_distance, key = "hn", adjustment = NULL, cutpoints = seq(from = 0, to = truncation_distance, by = i)))
    assign(paste0("detfc_hn_cos_trunc_",truncation_distance,"bin_",i), possDS(distdata, truncation_distance, key = "hn", adjustment = "cos", cutpoints = seq(from = 0, to = truncation_distance, by = i)))
    detfc_list_bin <- list(
      paste0("detfc_hr_null_trunc_",truncation_distance,"bin_",i),
      paste0("detfc_hn_null_trunc_",truncation_distance,"bin_",i),
      paste0("detfc_hn_cos_trunc_",truncation_distance,"bin_",i)
    )
  }

  
}



# clean lists
detfc_list <- compact(detfc_list)
detfc_list_bin <- compact(detfc_list_bin)
# compare models
do.call(summarize_ds_models, detfc_list)
do.call(summarize_ds_models, detfc_list_bin)

summarize_ds_models(detfc_list_bin[1], detfc_list_bin[2], detfc_list_bin[3])


## Estimate Detection Functions
# Store detection functions
detfc_list <- list()
# Null model (Hazard-rate)
detfc_hr_null <- Distance::ds(distdata, 0.5, key = "hr", adjustment = NULL)
detfc_hn_null <- Distance::ds(distdata, 0.5, key = "hn", adjustment = NULL, cutpoints = bins)
detfc_hn_null_cos <- Distance::ds(distdata, 0.5, key = "hn", adjustment = "cos", cutpoints = bins)
# summary(detfc_hr_null) #-2518.394
# summary(detfc_hn_null) #-2486.259
# summary(detfc_hn_null_cos) #-2536.719
par(mfrow = c(1, 3))
plot(detfc_hr_null, showpoints = FALSE, pl.den = 0, lwd = 2, breaks = bins, main = "Hazard-rate null model")
plot(detfc_hn_null, showpoints = FALSE, pl.den = 0, lwd = 2, main = "Half-normal null model")
plot(detfc_hn_null_cos, showpoints = FALSE, pl.den = 0, lwd = 2, main = "Half-normal null model")
ddf.gof(detfc_hr_null$ddf)
gof_ds(detfc_hr_null)
ddf.gof(detfc_hn_null$ddf)
ddf.gof(detfc_hn_null_cos$ddf)
par(mfrow = c(1, 1))
summarize_ds_models(detfc_hr_null, detfc_hn_null, detfc_hn_null_cos)
summarize_ds_models(detfc_hr_null_bins, detfc_hr_null)

# add to list
detfc_list[["Hazard-rate, null"]] <- detfc_hr_null
detfc_list[["Half-normal, null"]] <- detfc_hn_null
detfc_list[["Half-normal cos, null"]] <- detfc_hn_null_cos


detfc_hr_cc <- Distance::ds(distdata, 0.5, formula = ~as.factor(canopy_cover), key = "hr", adjustment = NULL)
detfc_hn_cc <- Distance::ds(distdata, 0.5, formula = ~as.factor(canopy_cover), key = "hn", adjustment = NULL)
detfc_hn_cc_cos <- Distance::ds(distdata, 0.5, formula = ~as.factor(canopy_cover), key = "hn", adjustment = "cos")
# summary(detfc_hr_cc) #-2443.976
# summary(detfc_hn_cc) #-2449.583
# summary(detfc_hn_cc_cos) #-2449.583
par(mfrow = c(1, 3))
plot(detfc_hr_null, showpoints = FALSE, pl.den = 0, lwd = 2, main = "Hazard-rate null model")
plot(detfc_hn_cc, showpoints = FALSE, pl.den = 0, lwd = 2, main = "Half-normal null model")
plot(detfc_hn_cc_cos, showpoints = FALSE, pl.den = 0, lwd = 2, main = "Half-normal null model")
ddf.gof(detfc_hr_null$ddf)
ddf.gof(detfc_hn_null$ddf)
ddf.gof(detfc_hn_null_cos$ddf)
par(mfrow = c(1, 1))


# add to list
detfc_list[["Hazard-rate, canopy cover"]] <- detfc_hr_cc
detfc_list[["Half-normal, canopy cover"]] <- detfc_hn_cc
detfc_list[["Half-normal cos, canopy cover "]] <- detfc_hn_cc_cos

table_ds_models <- summarize_ds_models(detfc_hr_null, detfc_hn_null, detfc_hn_null_cos, detfc_hr_cc, detfc_hn_cc, detfc_hn_cc_cos)
print(table_ds_models)
print(paste('Best model', table_ds_models[1,1]))

# add name to detection functions
for (i in 1:length(detfc_list)) {
  detfc_list[[i]]$name <- names(detfc_list)[i]
}


# Call summary and plot for each detection function after the loop
#for (model_name in names(detfc_list)) {
 # detfc <- detfc_list[[model_name]]
#  summary(detfc)
#  plot(detfc, showpoints = FALSE, pl.den = 0, lwd = 2, main = paste(model_name, "model"))
#  ddf.gof(detfc$ddf)
#}

## Fit and Analyse Distance Sampling Models
dsm_list <- list()

for (detfc in detfc_list) {
  # Uf detfc$name contains `null` then use quasipoisson, otherwise use Tweedie
    if (grepl("null", detfc$name)) {
      dsm_xy <- dsm::dsm(count ~ s(x, y), detfc, segdata, obsdata, method = "REML")
      dsm_list[[paste0('dsm_xy_',detfc$name)]] <- dsm_xy
    } else {
      dsm_est_model <- dsm(abundance.est ~ s(x, y), detfc, segdata, obsdata, method = "REML")
      dsm_list[[paste0('dsm_est_xy_',detfc$name)]] <- dsm_est_model
    }
}

# # Null model
# dsm.xy <- dsm::dsm(count ~ s(x, y), detfc.hr.null, segdata, obsdata, method = "REML")
# #summary(dsm.xy)
# #vis.gam(dsm.xy, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = GRID_SIZE)
# dsm_list[["null"]] <- dsm.xy
#
#
# # Models with covariates in the detection function
# for (model_name in names(models)) {
#   dsm_est_model <- dsm(abundance.est ~ s(x, y), detfc_list[[model_name]], segdata, obsdata, method = "REML")
#   dsm_list[[paste0('dsm.est.xy.',model_name)]] <- dsm_est_model
# }

# add name to dsm
for (i in 1:length(dsm_list)) {
  dsm_list[[i]]$name <- names(dsm_list)[i]
}


## Model Checking

# Call summary and plot for each detection function after the loop
#for (model_name in names(dsm_list)) {
#  dsm_model <- dsm_list[[model_name]]
#  summary(dsm_model)
#  plot(dsm_model, select = 2)
#  vis.gam(dsm_model, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = GRID_SIZE)
#}


# Check goodness of fit with Q-Q plots
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
for (model in dsm_list) {
 gam.check(model)
  # Add a title to the entire plotting area
  print(model$name)
  mtext(model$name, outer = TRUE, cex = 1.5)
}
par(mfrow = c(1, 1))


# Summarise model results
mod_results <- data.frame(
  "Model name" = names(dsm_list),
  "Description" = c(
    "Bivariate smooth of location, hazard-rate, quasipoisson",
    "Bivariate smooth of location, half-normal, quasipoisson",
    "Bivariate smooth of location, half-normal cos, quasipoisson",
    "Bivariate smooth of location, hazard-rate, quasipoisson, canopy height covariate in detection function",
    "Bivariate smooth of location, half-normal, quasipoisson, canopy height covariate in detection function",
    "Bivariate smooth of location, half-normal cos, quasipoisson, canopy height covariate in detection function"
  ),
  "Deviance explained" = sapply(
    dsm_list,
    function(x) paste0(round(summary(x)$dev.expl * 100, 2), "%")
  )
)
# wum_number <-'528' # 501, 503, 512 517, 528
# input_path <- here::here("Output", "DSM", paste0("dsm",wum_number,".RData"))
# load(file = input_path)
knitr::kable(mod_results, col.names = c("Model name", "Description", "Deviance explained"))

# save DSM list
output_path <- here("Output", "DSM", paste0("dsm",wum_number,".RData"))
save(dsm_list, detfc_list, file = output_path)

## Abundance Estimation

# Predict abundance using the null model
dsm_xy_pred <- predict(dsm_list[[1]], coords, dsm_list[[1]]$offset[1])

# Calculate total abundance over the survey area
total_abundance <- sum(dsm_xy_pred)

# Update density object with predicted values
density@density.surface[[1]]$density <- dsm_xy_pred

# Plot density surface
plot(density)
# plot(density@density.surface[[1]]['density'])

abudance_strata_list <- list()
for (strata_name in density@strata.name){
  density_strata <- subset(density@density.surface[[1]], strata == strata_name)
  abudance_strata_list[[strata_name]] <- as.numeric(sum(density_strata$density))
}

# Save processed data
output_path <- here("Output", "Density", paste0("density",wum_number,".RData"))
save(density, total_abundance, abudance_strata_list, region, file = output_path)

