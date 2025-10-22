library(dsims)
# library(knitr)
# library(RColorBrewer)
# library(dplyr)
library(sf)
# library(units)
# library(geosphere)


# Define the 10x10km study area
area_m <- matrix(c(0,0,10000,0,10000,10000,0,10000,0,0), ncol=2, byrow=TRUE)
area <- sf::st_polygon(list(area_m))

# Plot the study area and the polygons
# plot(area)

# Create regions for the study area and subplots
region <- make.region(region.name = "study area", shape = area, units = "m")
st_crs(region@region) <- 32633


# Define survey design
poly_dim <- c()
poly_dim$x_length <- 2000 # 2km x 2km plot

  
# Compute grid dimensions
n_side <- sqrt(16) # numer of subplots
if (n_side != floor(n_side)) stop("Number of plots must be a perfect square.")

# Calculate spacing between plots
total_plot_span <- n_side * poly_dim$x_length
total_gap_span <- 10000 - total_plot_span
gap <- total_gap_span / (n_side + 1)

# Create empty list to store polygons
plots <- list()

# Generate subplots
for (i in 0:(n_side - 1)) {
for (j in 0:(n_side - 1)) {
    xmin <- gap + i * (poly_dim$x_length + gap)
    ymin <- gap + j * (poly_dim$x_length + gap)
    xmax <- xmin + poly_dim$x_length
    ymax <- ymin + poly_dim$x_length
    plots[[length(plots) + 1]] <- st_polygon(list(rbind(
    c(xmin, ymin),
    c(xmin, ymax),
    c(xmax, ymax),
    c(xmax, ymin),
    c(xmin, ymin)
    )))
}
}

# Convert to sf object
subplots_sf <- st_sf(geometry = st_sfc(plots))
st_crs(subplots_sf) <- st_crs(region) # Set the same CRS as the region


# Create a subplot region using the generated polygons
subplots <- make.region(
region.name = "study area",
shape = subplots_sf
)

design <- make.design(
region = subplots,
transect.type = "line",
design = "systematic",
samplers = numeric(0),
line.length = numeric(0),
spacing = 80,
design.angle = 0,
edge.protocol = "minus",
truncation = 80
)

# create unique density
# Define a curved path (e.g., a semi-ellipse or arc)
centres <- list(
  c(5000, 3000),
  c(5789, 3108),
  c(6551, 3421),
  c(7236, 3906),
  c(7765, 4509),
  c(8039, 5165),
  c(7968, 5803),
  c(7506, 6354),
  c(6678, 6758),
  c(5591, 6972),
  c(4408, 6972),
  c(3321, 6758),
  c(2493, 6354),
  c(2031, 5803),
  c(1960, 5165),
  c(2234, 4509),
  c(2763, 3906),
  c(3448, 3421),
  c(4210, 3108),
  c(4999, 3000)
)


# Define and visualise uniform detection function
detectFun <- make.detectability(
  key.function = "uf",
  scale.param = 1, # accounting for canopy cover
  truncation = 80
)


# create design density
density <- make.density(region = design@region,
                    x.space = 500,
                    constant = 1)

# Add each hotspot along the curve
for (centre in centres) {
    density <- add.hotspot(object = density,
                            centre = centre,
                            sigma = 1500,
                            amplitude = 2)
}
density <- add.hotspot(object = density,
                    centre = c(5000, 8000),
                    sigma = 1500,
                    amplitude = -4)
plot(density, design@region, scale = 2)
lowDensity <- density
lowDensity@density.surface[[1]]$density <- lowDensity@density.surface[[1]]$density * 0.01 # fixes memory issues
plot(lowDensity, design@region, scale = 2)

ddfAnalyses <- make.ds.analysis(
                    dfmodel = ~1,
                    key = "hr",
                    criteria = "AIC",
                    truncation = 80,
                    group.strata = data.frame(design.id = design@region@strata.name, analysis.id = rep("A", length(design@region@strata.name)))
)

numIndividuals <- sample(10:30, length(design@region@strata.name), replace = TRUE)

# Create population description
popDesc <- make.population.description(
    region = design@region,
    density = density,
    N = numIndividuals,
    fixed.N = T
)

lowPopDesc <- make.population.description(
    region = design@region,
    density = lowDensity,
    N = numIndividuals,
    fixed.N = T
)

SIM_REPS <- 99

simHighDensity<- make.simulation(
    reps = SIM_REPS,
    design = design,
    population.description = popDesc,
    detectability = detectFun,
    ds.analysis = ddfAnalyses
)

simLowDensity <- make.simulation(
    reps = SIM_REPS,
    design = design,
    population.description = lowPopDesc,
    detectability = detectFun,
    ds.analysis = ddfAnalyses
)

simLowDensity <- run.simulation(simulation = simLowDensity, run.parallel = TRUE, max.cores = 20)
simHighDensity <- run.simulation(simulation = simHighDensity, run.parallel = TRUE, max.cores = 20)
