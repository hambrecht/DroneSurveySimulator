# Load required packages
library(here)       # For relative paths
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

### Set Variables

# Define grid cell size (in meters)
GRID_SIZE <- 500
wum_number <-'501' # 501, 503, 512 517, 528

## Load and Check Data

# Load processed data
input_path <- here("Output", "PrepData", paste0("prepared",wum_number,".RData"))
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
region <- make.region(
  region.name = "study area",
  shape = wmu
)

# Create density surface
density <- dsims::make.density(
  region = region,
  x.space = GRID_SIZE
)

# Extract coordinates
coords <- sf::st_drop_geometry(density@density.surface[[1]][, c("x", "y")])

## Estimate Detection Functions
# Store detection functions
detfc_list <- list()
detfc_list[["null"]] <- detfc.hr.null
# Null model (Hazard-rate)
detfc.hr.null <- Distance::ds(distdata, max(distdata$distance), key = "hr", adjustment = NULL)
#summary(detfc.hr.null)
#par(mfrow = c(1, 2))
#plot(detfc.hr.null, showpoints = FALSE, pl.den = 0, lwd = 2, main = "Null model")
#ddf.gof(detfc.hr.null$ddf)

# Detection function models with covariates
models <- list(
  canopy_height = ~ as.factor(canopy_height),
  canopy_cover = ~ as.factor(canopy_cover),
  agb = ~ as.factor(agb),
  vol = ~ as.factor(vol)
)



for (model_name in names(models)) {
  detfc <- Distance::ds(distdata, max(distdata$distance),
                        formula = models[[model_name]],
                        key = "hr", adjustment = NULL)
  detfc_list[[model_name]] <- detfc
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

# Null model
dsm.xy <- dsm::dsm(count ~ s(x, y), detfc.hr.null, segdata, obsdata, method = "REML")
#summary(dsm.xy)
#vis.gam(dsm.xy, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = GRID_SIZE)
dsm_list[["null"]] <- dsm.xy

# Tweedie model
dsm.xy.tweedie <- dsm(count ~ s(x, y), detfc.hr.null, segdata, obsdata, family = tw(), method = "REML")
#summary(dsm.xy.tweedie)
#vis.gam(dsm.xy.tweedie, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = GRID_SIZE)
dsm_list[["dsm.xy.tweedie"]] <- dsm.xy.tweedie

# Models with covariates
for (model_name in names(models)) {
  dsm_model <- dsm(count ~ s(x, y, k = 10) + s(get(model_name), k = 10), detfc.hr.null, segdata, obsdata, method = "REML")
  dsm_list[[paste0('dsm.xy.', model_name)]] <- dsm_model 
}

# Models with covariates in the detection function
for (model_name in names(models)) {
  dsm_est_model <- dsm(abundance.est ~ s(x, y), detfc_list[[model_name]], segdata, obsdata, method = "REML")
  dsm_list[[paste0('dsm.est.xy.',model_name)]] <- dsm_est_model
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
#par(mfrow = c(2, 2))
#for (model in dsm_list) {
#  gam.check(model)
#}
#par(mfrow = c(1, 1))

# Randomised quantile residuals for Tweedie model
#rqgam_check(dsm.xy.tweedie)

# Check for autocorrelation
#for (model in dsm_list) {
#  dsm_cor(model, max.lag = 10, Segment.Label = "Sample.Label")
#}

## Model Selection

# Summarise model results
mod_results <- data.frame(
  "Model name" = c(
    "`dsm.xy`", "`dsm.xy.tweedie`", "`dsm.xy.height`", "`dsm.xy.cover`", "`dsm.xy.agb`", "`dsm.xy.vol`",
    "`dsm.est.xy.height`", "`dsm.est.xy.cover`", "`dsm.est.xy.agb`", "`dsm.est.xy.vol`"
  ),
  "Description" = c(
    "Bivariate smooth of location, quasipoisson",
    "Bivariate smooth of location, Tweedie, quasipoisson",
    "Bivariate smooth of location, smooth of height, quasipoisson",
    "Bivariate smooth of location, smooth of cover, quasipoisson",
    "Bivariate smooth of location, smooth of AGB, quasipoisson",
    "Bivariate smooth of location, smooth of Volume, quasipoisson",
    "Bivariate smooth of location, Tweedie, height covariate in detection function",
    "Bivariate smooth of location, Tweedie, cover covariate in detection function",
    "Bivariate smooth of location, Tweedie, AGB covariate in detection function",
    "Bivariate smooth of location, Tweedie, Volume covariate in detection function"
  ),
  "Deviance explained" = sapply(
    dsm_list,
    function(x) paste0(round(summary(x)$dev.expl * 100, 2), "%")
  )
)

knitr::kable(mod_results, col.names = c("Model name", "Description", "Deviance explained"))

# save DSM list
output_path <- here("Output", "DSM", paste0("dsm",wum_number,".RData"))
save(dsm_list, detfc_list, file = output_path)

## Abundance Estimation

# Predict abundance using the null model
dsm.xy.pred <- predict(dsm.xy, coords, dsm.xy$offset[1])

# Calculate total abundance over the survey area
total_abundance <- sum(dsm.xy.pred)

# Update density object with predicted values
density@density.surface[[1]]$density <- dsm.xy.pred

# Plot density surface
plot(density@density.surface[[1]]['density'])

# Save processed data
output_path <- here("Output", "Density", paste0("density",wum_number,".RData"))
save(density, total_abundance, region, file = output_path)

