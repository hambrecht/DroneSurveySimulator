# Required packages
# library(devtools)
# devtools::install_github("r-lidar/lidR") 
# devtools::install_github("r-lidar/lidRviewer") 
# devtools::install_github("ptompalski/lidRmetrics") 
# install.packages("RCSF") # for cloth simulation filter

# load packages
library(lidR) # main package
library(lidRviewer) # includes better plotting function
library(lidRmetrics) # includes function for canopy cover
library(terra) # for writing raster to file

#load ALS file from drone
las <- readALS("D:/Sync/2_Areas/MooseProject/Data/DL_development/Hinton2025/lidar/lidars/terra_las/cloud0.las")

# validate las file
# las_check(las)
las <- filter_duplicates(las)

# plot las file
# view(las)

# ground classification with cloth simulation filter
las <- classify_ground(las, algorithm = csf())
## plot classified point cloud for quality control
# view(las, color = "Classification", size = 3, bg = "white") 

# filter noise
# las <- classify_noise(las)
# las <- remove_noise(las)


# height normalization
nlas <- normalize_height(las, knnidw())
## plot histogram
# hist(filter_ground(nlas)$Z, breaks = seq(-0.7, 0.7, 0.01), main = "", xlab = "Elevation")

# canopy cover raster
ccRaster <- pixel_metrics(nlas, ~metrics_percabove(Z, zmin = 0), res = 1)
# Plot the first layer of the raster (safer than hardcoding the name)
plot(ccRaster[[1]])

# export raster
output_dir <- "Output"
if (!dir.exists(output_dir)) dir.create(output_dir)
writeRaster(ccRaster$pzabovemean, filename = file.path(output_dir, "ccraster.tif"), filetype = "GTiff", overwrite = TRUE)

# canopy cover hexagon
#ccHex <- hexagon_metrics(nlas, ~metrics_percabove(Z, zmin = 0), area = 10)
#plot(ccHex$pzabovemean)
#writeRaster(ccHex$pzabovemean, filename = file.path(output_dir, "cchex.tif"), filetype = "GTiff", overwrite = TRUE)
#writeRaster(ccHex$pzabovemean, filename = "cchex.tif", filetype = "GTiff", overwrite = TRUE)

# retrieve canopy cover values for horse locations
# load horse locations
horseLocations <- vect("D:/Sync/2_Areas/MooseProject/Data/DL_development/Hinton2025/GPS/HorseLocations.shp")
# create 3m buffer around horse locations
horseLocationsBuffer <- buffer(horseLocations, width = 3)
# extract mean values from raster cells within the buffer
valsBuffer <- extract(ccRaster$pzabove2, horseLocationsBuffer, fun= "mean", method = "simple", na.rm = TRUE)
# add names to the locations
rownames(valsBuffer) <- values(horseLocations)$name
print(valsBuffer)

# save values to csv
write.csv(valsBuffer, file = file.path(output_dir, "horseCC.csv"), row.names = TRUE)
