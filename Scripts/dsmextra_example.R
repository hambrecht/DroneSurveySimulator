if (!require("remotes")) install.packages("remotes")
remotes::install_github("hambrecht/dsmextra")


#'---------------------------------------------
# Other required libraries
#'---------------------------------------------
library(dsmextra)     # Extrapolation toolkit for ecological models
library(raster)       # Geographic data analysis and modeling
library(tidyverse)    # Packages for data science
library(magrittr)     # Pipe operator
library(smoothr)

#'--------------------------------------------------------------------
# Set tibble options
#'--------------------------------------------------------------------
options(tibble.width = Inf) # All tibble columns shown
options(pillar.neg = FALSE) # No colouring negative numbers
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)

#'--------------------------------------------------------------------
# Set knitr options
#'--------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#'---------------------------------------------
# Load and extract the data
#'---------------------------------------------
data("spermwhales")
segs <- spermwhales$segs
predgrid <- spermwhales$predgrid

#'---------------------------------------------
# Quick inspection
#'---------------------------------------------
knitr::kable(head(segs[,c(1, 3:12)]), format = "pandoc")

knitr::kable(head(predgrid), format = "pandoc")

#'---------------------------------------------
# Create an outline of the study area boundaries
#'---------------------------------------------
study_area <- predgrid[, c("x", "y")]
study_area$value <- 1
study_area <- raster::rasterFromXYZ(study_area, crs = sp::CRS("+proj=aea +lat_1=38 +lat_2=30 +lat_0=34 +lon_0=-73 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
study_area <- raster::rasterToPolygons(study_area, dissolve = TRUE)
study_area <- sp::spTransform(study_area, CRSobj = sp::proj4string(transects))
study_area <- smoothr::smooth(study_area, method = "ksmooth", smoothness = 5)

#'---------------------------------------------
# Produce a simple plot
#'---------------------------------------------
plot(study_area, col = "lightskyblue1") # Region boundary
plot(transects, add = TRUE, col = "skyblue3") # Survey tracks
maps::map("world", fill = TRUE, col = "grey",
          xlim = range(obs$coords.x1),
          ylim = range(obs$coords.x2), add = TRUE)
pts <- obs # Sightings
coordinates(pts) <- ~ coords.x1 + coords.x2
axis(1); axis(2); box()
points(pts, pch = 16)

#'---------------------------------------------
# Define projected coordinate system
#'---------------------------------------------
aftt_crs <- sp::CRS("+proj=aea +lat_1=38 +lat_2=30 +lat_0=34 +lon_0=-73 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

#'---------------------------------------------
# Define environmental covariates of interest
#'---------------------------------------------
covariates.spermwhale <- c("Depth", "SST", "NPP", "DistToCAS", "EKE")

spermwhale.extrapolation <- compute_extrapolation(samples = segs,
                                                  covariate.names = covariates.spermwhale,
                                                  prediction.grid = predgrid,
                                                  coordinate.system = aftt_crs)

# Number of cells subject to univariate extrapolation (see below for definition)
predgrid %>%
  dplyr::filter(!dplyr::between(Depth, min(segs$Depth), max(segs$Depth)) |
                  !dplyr::between(SST, min(segs$SST), max(segs$SST)) |
                  !dplyr::between(NPP, min(segs$NPP), max(segs$NPP)) |
                  !dplyr::between(DistToCAS, min(segs$DistToCAS), max(segs$DistToCAS)) |
                  !dplyr::between(EKE, min(segs$EKE), max(segs$EKE))) %>%
  nrow()
str(spermwhale.extrapolation, 2)

compare_covariates(extrapolation.type = "both",
                   extrapolation.object = spermwhale.extrapolation,
                   n.covariates = NULL,
                   create.plots = TRUE,
                   display.percent = TRUE)

#'---------------------------------------------
# Calculate Gower's distances and %N
#'---------------------------------------------
spermwhale.nearby <- compute_nearby(samples = segs,
                                    prediction.grid = predgrid,
                                    coordinate.system = aftt_crs,
                                    covariate.names = covariates.spermwhale,
                                    nearby = 1)

#'---------------------------------------------
# Rename coordinates and convert to SpatialPointsdf
#'---------------------------------------------
obs.sp <- obs %>%
  dplyr::rename(., x = coords.x1, y = coords.x2) %>%
  sp::SpatialPointsDataFrame(coords = cbind(.$x, .$y), data = ., proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>%
  sp::spTransform(., CRSobj = aftt_crs)

map_extrapolation(map.type = "extrapolation",
                  extrapolation.object = spermwhale.extrapolation,
                  sightings = obs.sp,
                  tracks = transects)

map_extrapolation(map.type = "mic",
                  extrapolation.object = spermwhale.extrapolation,
                  sightings = obs.sp,
                  tracks = transects)


map_extrapolation(map.type = "nearby",
                  extrapolation.object = spermwhale.nearby,
                  sightings = obs.sp,
                  tracks = transects)


#'---------------------------------------------
# One function to rule them all, one funtion to find them, 
# One function to bring them all, and in the darkness bind them.
#'---------------------------------------------
spermwhale.analysis <- extrapolation_analysis(samples = segs,
                                              covariate.names = covariates.spermwhale,
                                              prediction.grid = predgrid,
                                              coordinate.system = aftt_crs,
                                              compare.covariates = TRUE,
                                              compare.extrapolation.type = "both",
                                              compare.n.covariates = NULL,
                                              compare.create.plots = TRUE,
                                              compare.display.percent = TRUE,
                                              nearby.compute = TRUE,
                                              nearby.nearby = 1,
                                              map.generate = TRUE,
                                              map.sightings = obs.sp,
                                              map.tracks = transects)

