# Load required packages
library(here)       # For relative paths
library(purrr)      # For functional programming
library(dplyr)      # For data manipulation
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
# Define grid cell size (in meters)
GRID_SIZE <- 500
wmu_number_list <- c('501', '503', '512', '517', '528')
# wmu_number <- wmu_number_list[1]
for (wmu_number in wmu_number_list){

  ## Load and Check Data

  # Load processed data
  input_path <- here("Output", "PrepData", paste0("prepared_",wmu_number,".RData"))
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
    shape = wmu
  )

  # Create density surface
  density <- dsims::make.density(
    region = region,
    x.space = GRID_SIZE
  )

  # Extract coordinates
  coords <- sf::st_drop_geometry(density@density.surface[[1]][, c("x", "y")])

  ## Distance Sampling Models
  truncation_distance_list <- c(round(max(distdata$distance),1), round(quantile(distdata$distance, 0.995),2))
  # detfc_list <- list()
  detfc_list_bin <- list()

  # Fixed bin size
  bin_size <- 0.1

  # Loop through truncation distances
  for (truncation_distance in truncation_distance_list) {
      detfc_list_bin[[paste0("detfc_hr_null_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, key = "hr", adjustment = NULL, cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
      detfc_list_bin[[paste0("detfc_hn_null_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, key = "hn", adjustment = NULL, cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
      detfc_list_bin[[paste0("detfc_hn_cos_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, key = "hn", adjustment = "cos", max_adjustments = 2,  monotonicity = 'strict', cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
      detfc_list_bin[[paste0("detfc_cc_hr_null_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, formula = ~as.factor(canopy_cover), key = "hr", adjustment = NULL, cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
      detfc_list_bin[[paste0("detfc_cc_hn_null_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, formula = ~as.factor(canopy_cover), key = "hn", adjustment = NULL, cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
      detfc_list_bin[[paste0("detfc_cc_hn_cos_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, formula = ~as.factor(canopy_cover), key = "hn", adjustment = "cos", max_adjustments = 2,  monotonicity = 'strict', cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
  }



  # clean lists
  # detfc_list_compact <- compact(detfc_list)
  detfc_list_bin_compact <- compact(detfc_list_bin)

# output_path <- here("Output", "DSM", paste0("detfc",wmu_number,".RData"))
# save(detfc_list_bin_compact, file = output_path)
# input_path <- here("Output", "DSM", paste0("detfc",wmu_number,".RData"))
# load(input_path)

  aic_values_bin <- purrr::map_dbl(detfc_list_bin_compact, ~ .x$ddf$criterion)
  # Extract CV values for each model
  cv_values_bin <- purrr::map_dbl(detfc_list_bin_compact, function(model) {
    model_summary <- summary(model)
    cv_value <- model_summary$ds$Nhat.se[1] / model_summary$ds$Nhat
    return(cv_value)
  })

  # Extract Chi-Square values for each model
  chi_square_values_bin <- purrr::map_dbl(detfc_list_bin_compact, function(model) {
    gof_test <- ddf.gof(model$ddf)
    chi_square_value <- gof_test$chisquare$chi1$chisq
    return(chi_square_value)
  })

  # Combine AIC and CV values into a dataframe
  model_metrics <- data.frame(
    Model = names(detfc_list_bin_compact),
    AIC = aic_values_bin,
    CV = cv_values_bin,
    Chi_Square = chi_square_values_bin
  )

  # Normalise AIC and CV (assuming lower values are better)
  model_metrics <- model_metrics %>%
    mutate(
      AIC_normalised = (AIC - min(AIC)) / (max(AIC) - min(AIC)),
      CV_normalised = (CV - min(CV)) / (max(CV) - min(CV)),
      Chi_Square_normalised = (Chi_Square - min(Chi_Square)) / (max(Chi_Square) - min(Chi_Square))
  )

  # Define weights for AIC, CV and Chi-Square (adjust these weights according to your preference)
  weight_AIC <- 1
  weight_CV <- 1
  weight_Chi_Square <- 1

  # Combine the normalized values using the weights
  model_metrics <- model_metrics %>%
    mutate(
      combined_score = weight_AIC * AIC_normalised + weight_CV * CV_normalised + weight_Chi_Square * Chi_Square_normalised
    ) %>%
    arrange(combined_score)

  row.names(model_metrics) <- NULL


  best_model_list <- list()
  best_model <- detfc_list_bin_compact[[model_metrics$Model[1]]]
  best_summary<- summary(best_model)
  best_cv <- (best_summary$ds$Nhat.se[1]/best_summary$ds$Nhat)
  best_model_list <- list(model = best_model, summary=best_summary, cv = best_cv)

  for (i in 1:length(detfc_list_bin_compact)) {
    detfc_list_bin_compact[[i]]$name <- names(detfc_list_bin_compact)[i]
  }

  ## Fit and Analyse Distance Sampling Models
  dsm_list <- list()

  for (detfc in detfc_list_bin_compact) {
    # Uf detfc$name contains `null` then use quasipoisson, otherwise use Tweedie
      if (grepl("cc", detfc$name)) {
        dsm_est_model <- dsm(abundance.est ~ s(x, y), detfc, segdata, obsdata, method = "REML")
        dsm_list[[paste0('dsm_est_xy_',detfc$name)]] <- dsm_est_model
      } else {
        dsm_xy <- dsm::dsm(count ~ s(x, y), detfc, segdata, obsdata, method = "REML")
        dsm_list[[paste0('dsm_xy_',detfc$name)]] <- dsm_xy
      }
  }

  # add name to dsm
  for (i in 1:length(dsm_list)) {
    dsm_list[[i]]$name <- names(dsm_list)[i]
  }

  # save DSM list
  output_path <- here("Output", "DSM", paste0("dsm",wmu_number,".RData"))
  save(dsm_list, detfc_list_bin_compact, best_model_list, file = output_path)
}

# Plot distance data
par(mfrow = c(1, 1))
break_bins <- seq(from = 0, to = 0.6, by = 0.05)
hist(distdata$distance, main = "Moose line transects", xlab = "Distance (km)", breaks = break_bins)
quantile(distdata$distance, 0.995)
max(distdata$distance)

# Plot the models
par(mfrow = c(2, 3))
for (model_name in model_metrics$Model) {
  model <- detfc_list_bin_compact[[model_name]]
  plot(model, showpoints = FALSE, pl.den = 0, lwd = 2, xlim = c(0,0.6), ylim = c(0,1), main = paste(model_name, "AIC:", round(model$ddf$criterion,2)))
  # ddf.gof(model$ddf, qq = TRUE, main = model_name)
}
par(mfrow = c(1, 1))

par(mfrow = c(1, 2))
plot(best_model, showpoints = FALSE, pl.den = 0, lwd = 2, ylim = c(0,1), main = names(sorted_aic_diff)[1])
test <- ddf.gof(best_model$ddf, asp=1, qq = TRUE)
test$chisquare$chi1$p
qqplot.ddf(best_model$ddf)
par(mfrow = c(1, 1))


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
    # "Description" = c(
    #   "Bivariate smooth of location, hazard-rate, quasipoisson",
    #   "Bivariate smooth of location, half-normal, quasipoisson",
    #   "Bivariate smooth of location, half-normal cos, quasipoisson",
    #   "Bivariate smooth of location, hazard-rate, quasipoisson, canopy height covariate in detection function",
    #   "Bivariate smooth of location, half-normal, quasipoisson, canopy height covariate in detection function",
    #   "Bivariate smooth of location, half-normal cos, quasipoisson, canopy height covariate in detection function"
    # ),
    "Deviance explained" = sapply(
      dsm_list,
      function(x) paste0(round(summary(x)$dev.expl * 100, 2), "%")
    )
  )
  # wmu_number <-'528' # 501, 503, 512 517, 528
  # input_path <- here::here("Output", "DSM", paste0("dsm",wmu_number,".RData"))
  # load(file = input_path)
  knitr::kable(mod_results, col.names = c("Model name", "Description", "Deviance explained"))


  ## Abundance Estimation

  # Predict abundance using the null model
  # Find the index of the sublist with name "Second"
index <- which(sapply(dsm_list, function(x) x$name) == "dsm_xy_detfc_hn_null_trunc_0.43_bin_0.1")

# Access the sublist

  sel_model <- dsm_list[[index]]
  dsm_xy_pred <- predict(sel_model, coords, sel_model$offset)

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
  output_path <- here("Output", "Density", paste0("density",wmu_number,".RData"))
  save(density, total_abundance, region, file = output_path)

