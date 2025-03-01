library(here)
library(Distance) # for summarize_ds_models function
library(dsm) # for dsm_cor function
library(mrds)
library(mgcv)
library(dplyr)      # For data manipulation
library(ggplot2)
library(RColorBrewer)

# Ensure the output directory exists
dir.create(here("Data", "Plot"), showWarnings = FALSE)
output_dir <- here("Data", "Plot")
# Define a colour palette
colour_palette <- brewer.pal(5, "RdBu")

# define list with WUM numbers
wmu_number_list <- c('501','503', '512', '517', '528')
wmu_number <- wmu_number_list[4]
# loop through WUM numbers
# for(wmu_number in wmu_number_list){

  # Load processed data
  input_path <- here("Output", "DSM", paste0("dsm",wmu_number,".RData"))
  load(file = input_path)

  # Assign individal names to data frames
  # assign(paste0("dsm_list_",wmu_number), dsm_list)
  # assign(paste0("detfc_list_bin_compact",wmu_number), detfc_list_bin_compact)
  # rm(dsm_list, detfc_list_bin_compact)
# }

  # Plot distance data
  
  length(detfc_list_bin_compact)
# remove `detfc_cc_hn_null_trunc_0.6_bin_0.05` from the list in case of error in `cv_values_bin`
# detfc_list_bin_compact <- detfc_list_bin_compact[-grep("cc_hn_null_trunc_0.6_bin_0.05", detfc_list_bin_compact)]
# detfc_list_bin_compact <- detfc_list_bin_compact[-grep("cc_hn_null_trunc_0.55_bin_0.05", detfc_list_bin_compact)]
# detfc_list_bin_compact

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
  AIC = round(aic_values_bin,0),
  CV = round(cv_values_bin,3),
  Chi_Square = round(chi_square_values_bin,1)
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
weight_CV <- 0.5
weight_Chi_Square <- 0.5

# Combine the normalized values using the weights
model_metrics <- model_metrics %>%
  mutate(
    combined_score = weight_AIC * AIC_normalised + weight_CV * CV_normalised + weight_Chi_Square * Chi_Square_normalised
  ) %>%
  arrange(combined_score)

row.names(model_metrics) <- NULL

knitr::kable(model_metrics, digits = 3, escape = FALSE, caption = "Model comparison metrics") # %>%
  # kable_styling(full_width = FALSE, position = "center")





# Select items containing 'trunc_0.6'
detfc_list_trunc_0.6 <- grep("trunc_0.6", detfc_list_bin_compact)

# Select items not containing 'trunc_0.6'
detfc_list_trunc_99 <- grep("trunc_0.6", detfc_list_bin_compact, invert = TRUE)


# for (wmu_number in wmu_number_list){
  table_ds_models_trunc_0.6 <- do.call(summarize_ds_models, detfc_list_bin_compact[detfc_list_trunc_0.6])
  table_ds_models_trunc_99 <- do.call(summarize_ds_models, detfc_list_bin_compact[detfc_list_trunc_99])
  knitr::kable(table_ds_models_trunc_0.6[, -1], digits=3, escape=FALSE,
             caption="Model selection summary of the 600m truncation models") # %>%
  # kable_styling(full_width = FALSE, position = "center")
  knitr::kable(table_ds_models_trunc_99[, -1], digits=3, escape=FALSE,
              caption="Model selection summary of the 99.5 percentile truncation models") # %>%
  # kable_styling(full_width = FALSE, position = "center")
  # }




### DF
# select best model
best_model <- detfc_list_bin_compact[[model_metrics$Model[1]]]
print(paste0("The best model is: ", model_metrics$Model[1]))
# Plot the models
par(mfrow = c(4, 3)) # adjust based on the number of models +1
break_bins <- seq(from = 0, to = 0.6, by = 0.05)
hist(distdata$distance, main = "Moose distance from line transects", xlab = "Distance (km)", breaks = break_bins)
for (model_name in model_metrics$Model) {
  model <- detfc_list_bin_compact[[model_name]]
  # Check if the current model is the best model
  if (model_name == model_metrics$Model[1]) {
    ggplot(model, showpoints = FALSE, pl.den = 0, lwd =2, xlim = c(0,0.6), ylim = c(0,1), main = paste(model_name, "AIC:", round(model$ddf$criterion,2)))
    # Draw a red box around the plot
    usr <- par("usr")
    rect(usr[1], usr[3], usr[2], usr[4], border = "red", lwd = 2)
  } else {
    ggplot(model, showpoints = FALSE, pl.den = 0, lwd = 2, xlim = c(0,0.6), ylim = c(0,1), main = paste(model_name, "AIC:", round(model$ddf$criterion,2)))
  }
  # ddf.gof(model$ddf, qq = TRUE, main = model_name)
}
par(mfrow = c(1, 1))




par(mfrow = c(1, 2))
plot(best_model, showpoints = FALSE, pl.den = 0, lwd = 2, ylim = c(0,1), main = model_metrics$Model[1])
qqplot.ddf(best_model$ddf)
par(mfrow = c(1, 1))


## DSM

# Call summary and plot for each detection function after the loop
par(mfrow = c(4, 3))
for (model_name in names(dsm_list)) {
 dsm_model <- dsm_list[[model_name]]
 vis.gam(dsm_model, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = 500, main = model_name)
}
par(mfrow = c(1, 1))

# Check for autocorrelation
par(mfrow = c(4, 3))
for (i in 1:length(dsm_list)) {
  dsm_cor(dsm_list[[i]], max.lag = 6, Segment.Label = "Sample.Label", main = names(dsm_list)[[i]])
}
par(mfrow = c(1, 1))

# Check goodness of fit with Q-Q plots
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  model <- dsm_list[[7]]
  gam.check(model)
  # Add a title to the entire plotting area
  print(model$name)
  mtext(model$name, outer = TRUE, cex = 1.5)

par(mfrow = c(1, 1))




# Summarise model results
mod_results <- data.frame(
  "Deviance explained" = sapply(
    dsm_list,
    function(x) paste0(round(summary(x)$dev.expl * 100, 2), "%")
  )
)

knitr::kable(mod_results, col.names = c("Model name", "Deviance explained"))
