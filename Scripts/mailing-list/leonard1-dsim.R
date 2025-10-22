##### By Eric:
# I'll start by describing the change to the second script, defining population and survey within poly_region first. This script works by adding a call to "run.simulation()" (rather than "run.survey()") near the bottom of the script.

# The first script, that you designed to directly address your extrapolation question, can be made to work (with the revisions I've included in the code), but the working code does not address the extrapolation issue.  By defining the survey within your subareas, the software cannot determine the detection probabilities of animals outside the subareas because the animals' distances from (non-existent) transects outside the subareas.

# Sadly, "dsims" cannot directly address the question that you want answered.




# Load the required library
library(dsims)
library(sf)
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

NOT_poly_design <- make.design(
  region = study_region, #poly_region,
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

plot(study_region, generate.transects(NOT_poly_design))

## Create Density
# Create a flat density grid
density <- make.density(region = study_region, x.space = 1, constant = 1)

# Add hotspots to the density grid for spatial variability
density <- add.hotspot(object = density, centre = c(5, 5), sigma = 10, amplitude = 2)
density <- add.hotspot(object = density, centre = c(15, 25), sigma = 5, amplitude = 3)
density <- add.hotspot(object = density, centre = c(40, 15), sigma = 15, amplitude = -0.5)

# Plot the density grid
plot(density, study_region)
density@density.surface[[1]] <- st_intersection(density@density.surface[[1]], poly_region@region)

# Create population descriptions for the study area and subplots to test both
study_area_pop_desc <- make.population.description(region = study_region, density = density)


# Create a detectability function
detect_uf <- make.detectability(
  key.function = "uf",
  scale.param = 0.8, 
  truncation = 2
)
plot(detect_uf, study_area_pop_desc)

# Create detection function analyses
ddf_analyses <- make.ds.analysis(dfmodel = ~1, key = c("hn", "hr"), criteria = "AIC", truncation = 2)

# Create a simulation for the subplots
sim <- make.simulation(
  reps = 99,
  design = NOT_poly_design, #poly_design,
  population.description = study_area_pop_desc, 
  detectability = detect_uf,
  ds.analysis = ddf_analyses
)

onlysamplesubsim <- make.simulation(
  reps = 99,
  design = poly_design,
  population.description = study_area_pop_desc, 
  detectability = detect_uf,
  ds.analysis = ddf_analyses
)



sim.out <- run.simulation(sim, run.parallel = TRUE)
onlysamplesubsim <- run.simulation(onlysamplesubsim, run.parallel = TRUE) #fails
# fails I believe because detection probabilities cannot be assigned to animals living
# outside "poly_region"

# Run the survey and summarize the results
summary(sim, use.max.reps = TRUE, description.summary = FALSE) # here were no successful repetitions.
summary(onlysamplesubsim, use.max.reps = TRUE, description.summary = FALSE) # here were no successful repetitions.

# Failing
survey <- run.survey(sim) # Error in sample.int(length(x), size, replace, prob) : incorrect number of probabilities

# Print the survey results
print(survey)

# Plot the survey results
plot(survey, study_region) # this works
