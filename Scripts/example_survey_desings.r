### Example 1: Creating density and population for study_region and run simulation

# Load the required library
library(dsims)

# Define the study area
area_m <- matrix(c(0,0,500,0,500,500,0,500,0,0), ncol=2, byrow=TRUE)
area <- sf::st_polygon(list(area_m))

# Plot the study area and the polygons
plot(area)

# Create regions for the study area and subplots
region <- make.region(region.name = "study area", shape = area)


## Helicopter design
heli_design <- make.design(
  region = region,
  transect.type = "line",
  design = "segmentedgrid",
  spacing = 12, # segments seperated by 1.2km
  seg.length = 100, # segements of 10km
  design.angle = 0, # align transect with north south
  seg.threshold = 10, # any segments less than 10% of the segment length (i.e. 1km) will be discarded.
  edge.protocol = "minus",
  truncation = 5, # IMAGE_WIDTH
)
heli_transects <- generate.transects(heli_design)

## Systematic design
sys_design <- make.design(
  region = region,
  transect.type = "line",
  design = "systematic",
  samplers = numeric(0), # OR
  line.length = heli_transects@line.length, # OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 5, # IMAGE_WIDTH
)
sys_transects <- generate.transects(sys_design)



## Random design
rnd_design <- make.design(
  region = region,
  transect.type = "line",
  design = "random",
  samplers = numeric(0), # OR
  line.length = heli_transects@line.length, # OR
  spacing = numeric(0),
  design.angle = 0,
  edge.protocol = "minus",
  truncation = 5, # IMAGE_WIDTH
)
rnd_transects <- generate.transects(rnd_design)



## Zigzag design
zigzag_design <- make.design(
  region = region,
  transect.type = "line",
  design = "eszigzag",
  samplers = numeric(0), # OR
  line.length = heli_transects@line.length, # OR
  spacing = numeric(0),
  design.angle = 90, # The design angle for the zigzag designs refers to the angle of a line which would run through the middle of each zigzag transect if the zigzags were to be generated within a rectangle. The design angle for zigzags should usually run along the longest dimension of the study region.
  edge.protocol = "minus",
  bounding.shape = "convex.hull", # rectangle or convex.hull. convex hull is generally more efficient.
  truncation = 5, # IMAGE_WIDTH
)
zigzag_transects <- generate.transects(zigzag_design)




## Zigzag with complementary line
zigzagcom_design <- make.design(
  region = region,
  transect.type = "line",
  design = "eszigzagcom", # eszigzag or eszigzagcom
  samplers = numeric(0), # OR
  line.length = heli_transects@line.length, # OR
  spacing = numeric(0),
  design.angle = 90,
  edge.protocol = "minus",
  bounding.shape = "convex.hull",
  truncation = 5, # IMAGE_WIDTH
)
zigzagcom_transects <- generate.transects(zigzagcom_design)


# Plot desings
par(mfrow = c(2, 3))
plot(region, heli_transects, lwd = 1, col = 4)
plot(region, sys_transects, lwd = 1, col = 4)
plot(region, rnd_transects, lwd = 1, col = 4)
plot(region, zigzag_transects, lwd = 1, col = 4)
plot(region, zigzagcom_transects, lwd = 1, col = 4)
par(mfrow = c(1, 1))
par(mfrow = c(2, 3))
plot(heli_design)
plot(sys_design)
plot(rnd_design)
plot(zigzag_design)
plot(zigzagcom_design)
par(mfrow = c(1, 1))