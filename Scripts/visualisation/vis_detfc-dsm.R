library(here)
library(mrds)
library(ggplot2)
library(RColorBrewer)

# Ensure the output directory exists
dir.create(here("Data", "Plot"), showWarnings = FALSE)

# Define a colour palette
colour_palette <- brewer.pal(5, "RdBu")

# define list with WUM numbers
wmu_number_list <- c('501','503', '512', '517', '528')

# loop through WUM numbers
for(wmu_number in wmu_number_list){
  # Load processed data
  input_path <- here("Output", "DSM", paste0("dsm",wmu_number,".RData"))
  load(file = input_path)

  # Assign individal names to data frames
  assign(paste0("dsm_list_",wmu_number), dsm_list)
  assign(paste0("detfc_list_",wmu_number), detfc_list)
  rm(dsm_list, detfc_list)
}

par(mfrow = c(2, 3))
for (wmu_number in wmu_number_list){
  for (i in length(paste0("detfc_list_",wmu_number))){
    plot(paste0("detfc_list_",wmu_number)[[i]] , showpoints = FALSE, pl.den = 0, lwd = 2, main = detfc_list501[[i]]$name)
    ddf.gof(paste0("detfc_list_",wmu_number)[[i]]$ddf)
  }
}
par(mfrow = c(1, 1))

for (wmu_number in wmu_number_list){
  table_ds_models <- do.call(summarize_ds_models, paste0("detfc_list_",wmu_number)[[i]])
  print(table_ds_models)
  print(paste('Best model', table_ds_models[1,1]))
}


## Model Checking

# Call summary and plot for each detection function after the loop
for (model_name in names(dsm_list)) {
  dsm_model <- dsm_list[[model_name]]
  summary(dsm_model)
  plot(dsm_model, select = 2)
  vis.gam(dsm_model, plot.type = "contour", view = c("x", "y"), asp = 1, type = "response", contour.col = "black", n.grid = GRID_SIZE)
}

# Check goodness of fit with Q-Q plots
par(mfrow = c(2, 2))
for (model in dsm_list) {
  gam.check(model)
}
par(mfrow = c(1, 1))

# Randomised quantile residuals for Tweedie model
rqgam_check(dsm.xy.tweedie)

# Check for autocorrelation
for (model in dsm_list) {
  dsm_cor(model, max.lag = 10, Segment.Label = "Sample.Label")
}

## Model Selection

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
knitr::kable(mod_results, col.names = c("Model name", "Description", "Deviance explained"))