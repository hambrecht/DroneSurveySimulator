# install packages
# install.packages("terra")
# First, ensure you have a copy of the devtools package
#if (!nzchar(system.file(package = "devtools"))) install.packages("devtools")
# then ensure you have a copy of the dssd package:
# if (!nzchar(system.file(package = "dssd"))) devtools::install_github("DistanceDevelopment/dssd", build_vignettes = TRUE)
# finally install dsims from github:
# devtools::install_github("DistanceDevelopment/dsims", build_vignettes = TRUE)

library(devtools)
library(terra)
library(dssd)
library(dsims)
library(sf)
library(lwgeom)

# Define the height and field of view
altitude <- 100  # height in meters
fov <- 60      # field of view in degrees

# Calculate the width using the tangent function in meters
image_wdith <- 2 * altitude * tan((fov * (pi / 180)) / 2)
print(image_wdith)

# https://distancesampling.org/R/vignettes/mexico-analysis.html




# Distance example using dsims
# https://examples.distancesampling.org/dsims-truncation/dsims-examples.html
# load shapefile
wum_strata <- sf::st_read("D:\\WMU\\base_data\\WMU\\wmu_501_3400.shp")
strata <- sf::st_read(wum_strata)

# replace numbers with letters
strata$prop_dec <- as.character(strata$prop_dec)

# Check for invalid polygons
validity <- st_is_valid(strata)
invalid_polygons <- strata[!validity, ]
fixed_polygons <- st_make_valid(invalid_polygons)
strata[!validity, ] <- fixed_polygons
all_valid <- st_is_valid(strata)
if (all(all_valid)) {
  print("All polygons are now valid.")
} else {
  print("There are still some invalid polygons.")
}


# Create the survey region
region <- make.region(region.name = "study area",
                      units = "m",
                      strata.name = "A",
                      shape = wum_strata)

region@region$geometry <- st_make_valid(region@region$geometry)
# The plot function allows plotting in km or m.
plot(region)



## # Create the density surface
density <- make.density(region = region,
                        x.space = 600,
                        constant = 1,
                        density.surface = dsm.xy)


# Add a hotspot to the density surface, centre located at x = 15000, y = 4000 with
# a Gaussian decay parameter sigma = 1500. The value at the centre point will now
# be 1 (the current value of the density surface defined above) + 0.5 = 1.5
# eg.density <- add.hotspot(density, centre = c(340000,6000000), sigma = 10000, amplitude = 0.5)
# Add a lowspot to this new density surface (eg.density)
# eg.density <- add.hotspot(eg.density, centre = c(285000,5970000), sigma = 1000, amplitude = -0.25)
# Plot the density surface
plot(density, region)




# Create the population description, with a population size N = 200
pop.desc <- make.population.description(region = region,
                                        density = density ,
                                        N = rep(200, length(region@strata.name)),
                                        fixed.N = TRUE)


# Create the covariate list
covariate.list <- list()
# The population will be 50% males and 50% females
covariate.list$sex <- list(data.frame(level = c("female", "male"),
                                      prob = c(0.5,0.5)))


## # Create the population description, with a population size N = 200
## pop.desc.cov <- make.population.description(region = region,
##                                             density = density,
##                                             covariates = covariate.list,
##                                             N = 200)

# Create the population description, with a population size N = 200
pop.desc.cov <- make.population.description(region = region,
                                            density = density,
                                            covariates = covariate.list,
                                            N = 200)

covariate.list <- list()
covariate.list$size <- list(list(distribution = "poisson", lambda = 35))



# Make a simple half normal detection function with a scale parameter of 200
detect.hn <- make.detectability(key.function = "hn",
                                 scale.param = 200,
                                 truncation = image_wdith)
# We can now visualise these detection functions
plot(detect.hn, pop.desc)

# Create the covariate parameter list
cov.params <- list()
# Note the covariate parameters are supplied on the log scale
cov.params$sex = data.frame(level = c("female", "male"),
                            param = c(0, 1.5))

detect.cov <- make.detectability(key.function = "hn" ,
                                 scale.param = 120,
                                 cov.param = cov.params,
                                 truncation = 1000)



# We can now visualise these detection functions
plot(detect.cov, pop.desc.cov)



# Define the design
design <- make.design(region = region,
                      transect.type = "line",
                      design = "systematic",
                      spacing = 1000,
                      truncation = image_wdith)

parallel.design <- make.design(region = region,
                               design = "systematic",
                               spacing = 2500,
                               edge.protocol = "minus",
                               design.angle = 90,
                               truncation = image_wdith)

zigzag.design <- make.design(region = region,
                             design = "eszigzag",
                             spacing = 1000,
                             edge.protocol = "minus",
                             design.angle = 0,
                             bounding.shape = "convex.hull",
                             truncation = image_wdith)


transects <- generate.transects(design)
parallel.transects <- generate.transects(parallel.design)
zigzag.transects <- generate.transects(zigzag.design)
plot(region, transects)
plot(region, parallel.transects)
plot(region, zigzag.transects)


ddf.analyses <- make.ds.analysis(dfmodel = list(~1, ~1),
                                 key = c("hn", "hr"),
                                 criteria = "AIC",
                                 truncation = 600)



## sim <- make.simulation(reps = 999,
##                        region.obj = region,
##                        design.obj = design,
##                        detectability.obj = detect.hn,
##                        ddf.analyses.list = ddf.analyses)
##                        population.description.obj = pop.desc,
## # Produce simulation setup plots
## check.sim.setup(sim)

sim <- make.simulation(reps = 999,
                       design = design,
                       population.description = pop.desc,
                       detectability = detect.hn,
                       ds.analysis = ddf.analyses)
# Produce survey and plot it
survey <- run.survey(sim)
plot(survey, region)

sim <- run.simulation(simulation = sim, run.parallel = F)
# Display a summary of the simulation
summary(sim)
# Display a histogram of the estimates of abundance
histogram.N.ests(sim)



