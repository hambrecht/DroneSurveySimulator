# Load necessary libraries
library(here)
library(dsims)
library(knitr)

# Check if pbapply is installed
if (!requireNamespace("pbapply", quietly = TRUE)) {
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
#' @param FOV Numeric value representing the field of view of the camera in degrees. Default is 60 degrees.
#'
#' @return Numeric value representing the rounded image width in meters.
#' @throws Error if ALTITUDE or FOV are not positive.
#'
#' @examples
#' calculate_image_width(ALTITUDE = 100)
#' calculate_image_width(ALTITUDE = 150, FOV = 75)
#'
#' @export
calculate_image_width <- function(ALTITUDE, FOV = 25) {
  if (ALTITUDE <= 0 || FOV <= 0) stop("ALTITUDE and FOV must be positive numbers.")
  
  # Calculate and round the image width
  round(2 * ALTITUDE * tan((FOV * pi / 180) / 2), -1)
}

# Load density data
input_path <- here::here("Output", "Density", "density501.RData")
load(file = input_path)

# Define constants
ALTITUDE <- 120          # Height in meters
HFOV <- 25             # Horizontal FOV in degrees
IMAGE_WIDTH <- calculate_image_width(ALTITUDE, HFOV)
print(paste0("Half swath width is: ",IMAGE_WIDTH,"m"))


# Create population description
pop_desc <- make.population.description(
  region = region,
  density = density,
  N = total_abundance,
  fixed.N = TRUE
)

# Define and visualise detection function
detect_hn <- make.detectability(
  key.function = "hn",
  scale.param = 100,
  truncation = IMAGE_WIDTH
)
plot(detect_hn, pop_desc)

# Define and visualise uniform detection function
detect_uf <- make.detectability(
  key.function = "uf",
  scale.param = 0.9, # accounting for canopy cover
  truncation = IMAGE_WIDTH
)
plot(detect_uf, pop_desc)

# create coverage grid
cover <- make.coverage(region,
                       n.grid.points = 1000)
plot(region, cover)

# Define survey design
design <- make.design(
  region = region,
  transect.type = "line",
  design = "eszigzag",
  samplers = numeric(0),
  line.length = numeric(0),
  seg.length = numeric(0),
  effort.allocation = numeric(0),
  design.angle = 90,
  spacing = 1000,
  edge.protocol = "minus",
  seg.threshold = numeric(0),
  bounding.shape = "convex.hull",
  truncation = IMAGE_WIDTH,
  coverage.grid = NULL
)
default.design <- make.design(region = region,
                              transect.type = "line",
                              design = "systematic",
                              samplers = 20,
                              design.angle = c(155, 90),
                              edge.protocol = "minus",
                              truncation = 2000, 
                              coverage.grid = cover)
transects <- generate.transects(default.design)
plot(region, transects, lwd = 1.5, col = 4)


plot(region, transects, lwd = 1.5, col = 4)

# You can re-run the coverage simulation using the following code. Note, your
# results should vary slightly from mine, make sure you haven't set a seed!
design <- run.coverage(design, reps = 1000)
plot(design)

design <- make.design(region = region,
                      transect.type = "line",
                      design = "segmentedgrid",
                      spacing = c(1200, 1200, 1200), # segments seperated by 1.2km
                      seg.length = c(10000,10000, 10000), # segements of 10km
                      design.angle = c(0, 0, 0),
                      seg.threshold = c(10 ,10 ,10), # any segments less than 10% of the segment length (i.e. 1km) will be discarded.
                      edge.protocol = "minus",
                      truncation = 600,
                      coverage.grid = cover)
transects <- generate.transects(design)
plot(region, transects, lwd = 0.5, col = 4)

design <- run.coverage(design, reps = 1000)
plot(design)




# Visualise the survey design
transects <- generate.transects(design)
plot(region, transects)

# Define analysis models
ddf_analyses <- make.ds.analysis(
  dfmodel = ~1,
  key = "hn",
  criteria = "AIC",
  truncation = 600
)
ddf_analyses <- make.ds.analysis(
  dfmodel = list(~1, ~1),
  key = c("hn", "hr"),
  criteria = "AIC",
  truncation = 600
)

# Create and run the simulation
sim <- make.simulation(
  reps = 99,
  design = design,
  population.description = pop_desc,
  detectability = detect_uf,
  ds.analysis = ddf_analyses
)

survey <- run.survey(sim)
plot(survey, region)

# Run the full simulation
sim <- run.simulation(simulation = sim, run.parallel = F)

# Display results
summary(sim)
histogram.N.ests(sim, xlim=c(7500,9500))



## Testing truncation distances. Required with drone???
# Investigate truncation distances
truncation_distances <- c(
  calculate_image_width(100), calculate_image_width(200),
  calculate_image_width(300), calculate_image_width(400),
  calculate_image_width(500)
)

results_list <- vector("list", length(truncation_distances))
summary_list <- vector("list", length(truncation_distances))

for (i in seq_along(truncation_distances)) {
  cat(sprintf("\nRunning for truncation = %d", truncation_distances[i]))
  
  new_ds_analyses <- make.ds.analysis(
    dfmodel = list(~1, ~1),
    key = c("hn", "hr"),
    criteria = "AIC",
    truncation = truncation_distances[i]
  )
  
  sim@ds.analysis <- new_ds_analyses
  results_list[[i]] <- run.simulation(sim, run.parallel = F)
  summary_list[[i]] <- summary(results_list[[i]], description.summary = FALSE)
}

names(results_list) <- paste0("t", truncation_distances)
names(summary_list) <- paste0("t", truncation_distances)


# Extracting results statistics

N    <- unlist(lapply(summary_list, function(x){x@individuals$N$mean.Estimate}))
n    <- unlist(lapply(summary_list, function(x){x@individuals$summary$mean.n}))
se   <- unlist(lapply(summary_list, function(x){x@individuals$N$mean.se}))
sd_N <- unlist(lapply(summary_list, function(x){x@individuals$N$sd.of.means}))
bias <- unlist(lapply(summary_list, function(x){x@individuals$N$percent.bias}))
RMSE <- unlist(lapply(summary_list, function(x){x@individuals$N$RMSE}))
cov  <- unlist(lapply(summary_list, function(x){x@individuals$N$CI.coverage.prob}))

sim_data <- data.frame(trunc = truncation_distances,
                       n = round(n),
                       N = round(N),
                       se = round(se,2),
                       sd.N = round(sd_N,2),
                       bias = round(bias,2),
                       RMSE = round(RMSE,2),
                       cov = round(cov*100,1))

kable(sim_data,
      col.names = c("$Truncation$", "$mean\\ n$", "$mean\\ \\hat{N}$", "$mean\\ se$", "$SD(\\hat{N})$", "$\\% Bias$", "$RMSE$", "$\\%\\ CI\\ Coverage$"),
      row.names = FALSE,
      align = c('c', 'c', 'c', 'c', 'c', 'c', 'c', 'c'),
      caption = "Simulation Results for the simple half normal detection probability: The truncation distance, mean number of detections, mean estimated population size (N), mean standard error of $\\hat{N}$, the standard deviation of $\\hat{N}$, percentage bias, root mean squared error, percentage of times the true value of N was captured in the confidence intervals.",
      table.placement="!h",
      format = "simple")
