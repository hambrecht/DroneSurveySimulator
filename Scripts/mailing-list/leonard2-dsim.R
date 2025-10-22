##### By Eric:
# I'll start by describing the change to the second script, defining population and survey within poly_region first. This script works by adding a call to "run.simulation()" (rather than "run.survey()") near the bottom of the script.


library(dsims)
# Create Multi Strata Region
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

# Create regions for the study area and subplots
poly_region <- make.region(region.name = "sub plots", strata.name = strata.names, shape = mp1)

# Plot the regions
plot(poly_region)

# Create a systematic subplot design
poly_design <- make.design(
  region = poly_region,
  transect.type = "line",
  design = "systematic",
#  samplers = numeric(0), # OR
#  line.length = numeric(0), # OR
  spacing = 3,
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 2
)
# Generate transects for the subplot design
poly_transects <- generate.transects(poly_design)

# Plot the transects on the subplot region
plot(poly_region, poly_transects, lwd = 0.5, col = 4, covered.area=TRUE)

## Create Density
# Create a flat density grid
density_poly <- make.density(region = poly_region, x.space = 1, constant = 1)
density_poly@density.surface
# Add hotspots to the density grid for spatial variability
#density_poly <- add.hotspot(object = density_poly, centre = c(5, 5), sigma = 10, amplitude = 2)
#density_poly <- add.hotspot(object = density_poly, centre = c(15, 25), sigma = 5, amplitude = 3)
#density_poly <- add.hotspot(object = density_poly, centre = c(40, 15), sigma = 15, amplitude = -0.5)

# Plot the density grid
#plot(density_poly, poly_region)

# Create population descriptions for the study area and subplots to test both
poly_pop_desc <- make.population.description(region = poly_region, density = density_poly)

# Create a detectability function
detect_uf <- make.detectability(
  key.function = "uf",
  scale.param = 0.8, 
  truncation = 2
)
plot(detect_uf, poly_pop_desc)

# Create detection function analyses
ddf_analyses_poly <- make.ds.analysis(dfmodel = ~1, key = c("hn", "hr"), criteria = "AIC", truncation = 2, group.strata = data.frame(design.id = poly_region@strata.name, analysis.id = rep("A", length(poly_region@strata.name))))

# Create a simulation for the subplots
poly_sim <- make.simulation(
  reps = 99,
  design = poly_design,
  population.description = poly_pop_desc, 
  detectability = detect_uf,
  ds.analysis = ddf_analyses_poly 
)

#------------- this is the missing step
poly_sim.out <- run.simulation(poly_sim, run.parallel = TRUE)
summary(poly_sim.out)

# Run the survey and summarize the results
summary(poly_sim, use.max.reps = TRUE, description.summary = FALSE) # here were no successful repetitions.

# Failing
poly_survey <- run.survey(poly_sim) # Error in sample.int(length(x), size, replace, prob) : incorrect number of probabilities

# Print the survey results
print(poly_survey)

# Plot the survey results
plot(poly_survey, poly_region)
