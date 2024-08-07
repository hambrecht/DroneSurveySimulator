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
  assign(paste0("dsm_list",wmu_number), dsm_list)
  assign(paste0("detfc_list",wmu_number), detfc_list)
  rm(dsm_list, detfc_list)
}

detfc_list501$canopy_cover$ddf

summary(detfc_list501$canopy_cover)
par(mfrow = c(1, 2))
plot(detfc_list501$canopy_cover, showpoints = FALSE, pl.den = 0, lwd = 2, main = "Canopy cover model")
ddf.gof(detfc_list501$canopy_cover$ddf)

# Detection function models with covariates
models <- list(
  canopy_height = ~ as.factor(canopy_height),
  canopy_cover = ~ as.factor(canopy_cover),
  agb = ~ as.factor(agb),
  vol = ~ as.factor(vol)
)

# Call summary and plot for each detection function after the loop
for (model_name in names(detfc_list)) {
detfc <- detfc_list[[model_name]]
 summary(detfc)
 plot(detfc, showpoints = FALSE, pl.den = 0, lwd = 2, main = paste(model_name, "model"))
 ddf.gof(detfc$ddf)
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