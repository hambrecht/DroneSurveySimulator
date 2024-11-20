# Load the required library
library(dsims)

# Create a detectability function
detect_uf <- make.detectability(
  key.function = "uf",
  scale.param = 0.8, 
  truncation = 2
)

# Define the outer boundary and holes for the polygons
outer <- matrix(c(0,0,15,0,15,10,0,10,0,0), ncol=2, byrow=TRUE)
hole1 <- matrix(c(2,2,2,3,3,3,3,2,2,2), ncol=2, byrow=TRUE)
hole2 <- matrix(c(5,5,5,6,7,6,8,5.5,7,5,5,5), ncol=2, byrow=TRUE)

# Create polygons with holes
pol1 <- sf::st_polygon(list(outer, hole1*1.5, hole2))
pol2 <- sf::st_polygon(list(outer + 15, hole2*1.5 + 12))
pol3 <- sf::st_polygon(list(outer + 30, hole2*2.5 + 20))

# Combine polygons into a simple feature collection
sfc <- sf::st_sfc(pol1, pol2, pol3)
strata.names <- c("SW", "central", "NE")
mp1 <- sf::st_sf(strata = strata.names, geom = sfc)

# Print the simple feature collection
print(mp1)

# Define the study area
area_m <- matrix(c(0,0,45,0,45,40,0,40,0,0), ncol=2, byrow=TRUE)
area <- sf::st_polygon(list(area_m))

# Plot the study area and the polygons
plot(area)
plot(mp1, add = TRUE)

# Create regions for the study area and subplots
study_region <- make.region(region.name = "study area", shape = area)
poly_region <- make.region(region.name = "sub plots", strata.name = strata.names, shape = mp1)

# Plot the regions
plot(study_region)
plot(poly_region)

# Create a flat density grid
density <- make.density(region = study_region, x.space = 1, constant = 1)

# Add hotspots to the density grid for spatial variability
density <- add.hotspot(object = density, centre = c(5, 5), sigma = 10, amplitude = 2)
density <- add.hotspot(object = density, centre = c(15, 25), sigma = 5, amplitude = 3)
density <- add.hotspot(object = density, centre = c(40, 15), sigma = 15, amplitude = -0.5)

# Plot the density grid
plot(density, study_region)

# Create a covariate list describing the distribution of cluster sizes
covariates <- list(cover = list(distribution = "ztruncpois", mean = 3))

# Create population descriptions for the study area and subplots
study_area_pop_desc_co <- make.population.description(region = study_region, density = density, covariates = covariates)
sub_plots_pop_desc_co <- make.population.description(region = poly_region, density = density, covariates = covariates)

# Create a systematic subplot design
subplots_design <- make.design(
  region = poly_region,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = numeric(0), # OR
  spacing = 3,
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 2
)

# Generate transects for the subplot design
subplots_transects <- generate.transects(subplots_design)

# Plot the transects on the subplot region
plot(poly_region, subplots_transects, lwd = 0.5, col = 4)

# Create detection function analyses
ddf_analyses <- make.ds.analysis(dfmodel = ~1, key = c("hn", "hr"), criteria = "AIC", truncation = 2)
ddf_analyses_sub <- make.ds.analysis(dfmodel = ~1, key = c("hn", "hr"), criteria = "AIC", truncation = 2, group.strata = data.frame(design.id = poly_region@strata.name, analysis.id = rep("A", length(poly_region@strata.name))))
ddf_analyses_co <- make.ds.analysis(dfmodel = list(~1, ~1, ~cover, ~cover), key = c("hn", "hr", "hn", "hr"), criteria = "AIC", truncation = 2)
ddf_analyses_sub_co <- make.ds.analysis(dfmodel = list(~1, ~1, ~cover, ~cover), key = c("hn", "hr", "hn", "hr"), criteria = "AIC", truncation = 2, group.strata = data.frame(design.id = poly_region@strata.name, analysis.id = rep("A", length(poly_region@strata.name))))

# Create a simulation for the subplots
sim_sub <- make.simulation(
  reps = 99,
  design = subplots_design,
  population.description = sub_plots_pop_desc_co, # study_area_pop_desc_co, sub_plots_pop_desc_co
  detectability = detect_uf,
  ds.analysis = ddf_analyses_sub_co # ddf_analyses, ddf_analyses_sub, ddf_analyses_co, ddf_analyses_sub_co
)

# Run the survey and summarize the results
summary(sim_sub, use.max.reps = TRUE, description.summary = FALSE) # here were no successful repetitions.

# Failing
sub_survey <- run.survey(sim_sub) # Error in sample.int(length(x), size, replace, prob) : incorrect number of probabilities

# Print the survey results
print(sub_survey)

# Plot the survey results
plot(sub_survey, subplots)