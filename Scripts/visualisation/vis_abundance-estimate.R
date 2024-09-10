library(here)
library(raster)
library(sf)
library(rayshader)
library(dssd)
library(dsims)



# load processed data
# define list with WMU numbers
wmu_number_list <- c('501','503', '512', '517', '528')

# loop through WMU numbers
for(wmu_number in wmu_number_list){
  # Load processed data
  input_path <- here("Output", "Density", paste0("density",wmu_number,".RData"))
  load(file = input_path)

  # Assign individual names to data frames
  assign(paste0("density_",wmu_number), density)
  assign(paste0("total_abundance_",wmu_number), total_abundance)
  assign(paste0("region_",wmu_number), region)
  rm(density, total_abundance, region)
}
## Abundance Estimation
# Plot density surface
density_list <- paste0("density_",wmu_number)
total_abundance_list <- paste0("total_abundance_",wmu_number)


for (wmu_number in wmu_number_list){
  density <- get(paste0("density_",wmu_number))
  total_abundance <- get(paste0("total_abundance_",wmu_number))
  plot(density)
}

plot(density_501)

# Rayshader
#And convert it to a matrix:
# Extract the geometry and density columns from the density.surface table
geometry <- density_501@density.surface[[1]]$geometry
density_values <- 1000* density_501@density.surface[[1]]$density

# Convert the geometry column to a spatial object
spatial_geometry <- st_as_sf(geometry)

# Add density values as an attribute to the spatial object
spatial_geometry$density <- density_values

# Convert sf object to Spatial (sp) object because rasterize in raster package works with sp objects
spatial_geometry_sp <- as(spatial_geometry, "Spatial")

# Create an empty raster object with the desired extent and resolution
raster_template <- raster(extent(spatial_geometry_sp), resolution = 10)

# Check CRS and set it to match spatial data if not already set
crs(raster_template) <- crs(spatial_geometry_sp)

# Rasterize the spatial object, assigning the density values to the raster cells
raster_geometry <- rasterize(spatial_geometry_sp, raster_template, field = "density")


# Extract geometry and density, convert to spatial object, and rasterize
library(sf)
library(raster)

geometry <- density_501@density.surface[[1]]$geometry
x <- density_501@density.surface[[1]]$x
y <- density_501@density.surface[[1]]$y
density <- 100 * density_501@density.surface[[1]]$density

# Create a SpatialPointsDataFrame with the original coordinates
coords <- data.frame(x = x, y = y)
data <- data.frame(density = density)
spatial_points <- SpatialPointsDataFrame(coords, data, proj4string = CRS("+proj=longlat +datum=WGS84"))

# Extract the extent of the SpatialPointsDataFrame
points_extent <- extent(spatial_points)

# Define the raster with the same extent as the SpatialPointsDataFrame
r <- raster(extent(points_extent), res = 500)

# Ensure the raster has the same CRS as the spatial points
crs(r) <- crs(spatial_points)

# Rasterize the spatial points
rasterized <- rasterize(spatial_points, r, field = "density", fun = mean, background = NA)

# Check and plot the raster
print(rasterized)
plot(rasterized)


# Convert the raster object to a matrix
elmat <- raster_to_matrix(rasterized)


#We use another one of rayshader's built-in textures:
elmat %>%
  sphere_shade(texture = "imhof1", colorintensity=0.1) %>% # “imhof1”, “imhof2”, “imhof3”, imhof4“,”desert“,”bw“,”unicorn”.
  plot_map()


elmat %>%
  height_shade() %>%
  # sphere_shade(texture = "imhof1", colorintensity=0.1) %>%
  # add_shadow(ray_shade(elmat, zscale = 3), 0.5) %>%
  add_shadow(ambient_shade(elmat), 0) %>%
  plot_3d(elmat, zscale = 10, fov = 0, theta = 135, zoom = 0.75, phi = 45, windowsize = c(1000, 800))
Sys.sleep(0.2)
render_scalebar(limits=c(0, 5, 10),label_unit = "km",position = "W", y=10,
                scale_length = c(0.33,1))
render_compass(position = "NE")
render_snapshot()
