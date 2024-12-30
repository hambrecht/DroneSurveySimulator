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

wmu_number_list <- c('501', '503', '512', '517', '528')
wmu_number <- wmu_number_list[1]
# for (wmu_number in wmu_number_list){
  print(paste0("Processing WMU: ", wmu_number))

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

  


  ## Distance Sampling Models
  truncation_distance_list <- c(round(max(distdata$distance),1), round(quantile(distdata$distance, 0.995),2))
  # detfc_list <- list()
  detfc_list_bin <- list()

  # Fixed bin size
  bin_size <- 0.05

  # Loop through truncation distances
  for (truncation_distance in truncation_distance_list) {
      detfc_list_bin[[paste0("detfc_hr_null_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, key = "hr", adjustment = NULL, cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
      detfc_list_bin[[paste0("detfc_hn_null_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, key = "hn", adjustment = NULL, cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
      detfc_list_bin[[paste0("detfc_hn_cos_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, key = "hn", adjustment = "cos", max_adjustments = 2,  monotonicity = 'strict', cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
      detfc_list_bin[[paste0("detfc_cc_hr_null_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, formula = ~as.factor(canopy_cover), key = "hr", adjustment = NULL, cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
      detfc_list_bin[[paste0("detfc_cc_hn_null_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, formula = ~as.factor(canopy_cover), key = "hn", adjustment = NULL, cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
      # detfc_list_bin[[paste0("detfc_cc_hn_cos_trunc_", truncation_distance, "_bin_", bin_size)]] <- possDS(distdata, truncation_distance, formula = ~as.factor(canopy_cover), key = "hn", adjustment = "cos", max_adjustments = 2,  monotonicity = 'strict', cutpoints = seq(from = 0, to = truncation_distance, by = bin_size))
  }


  # clean lists
  # detfc_list_compact <- compact(detfc_list)
  detfc_list_bin_compact <- compact(detfc_list_bin)

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
  save(dsm_list, detfc_list_bin_compact, best_model_list, wmu, distdata, file = output_path)
# }



