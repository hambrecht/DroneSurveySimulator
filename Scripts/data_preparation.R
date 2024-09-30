# Load necessary libraries
library(here)    # Relative paths
library(sf)      # Simple Features for R
library(dplyr)   # Data manipulation
library(tidyr)   # Data tidying
library(units)   # Units for spatial data
library(purrr)   # Functional programming tools

# Define functions

#' Load and Transform Spatial Data
#'
#' This function reads spatial data from specified file paths and transforms
#' it to the given coordinate reference system (CRS). It returns a list of
#' spatial data frames.
#'
#' @param moose_path A character string specifying the file path to the moose 
#'   location data.
#' @param transects_path A character string specifying the file path to the
#'   transect data.
#' @param sbfi_path A character string specifying the file path to the SBFI 
#'   data.
#' @param wmu_path A character string specifying the file path to the WMU 
#'   data.
#' @param crs A numeric or character string specifying the CRS to transform
#'   the data to.
#' @param transects2_path Optional. A character string specifying the file path to 
#'   a second transect dataset. If provided, it will be joined with the first transect.
#' @return A named list containing the transformed spatial data frames for moose,
#'   transects, sbfi, and wmu.
#' @export
load_spatial_data <- function(moose_path, sbfi_path, wmu_path, crs, transects_path, transects2_path = NULL) {
  # Helper function to safely read spatial data
  safe_st_read <- function(path, layer = NULL) {
    if (file.exists(path)) {
      if (!is.null(layer)) {
        sf::st_read(path, layer = layer)
      } else {
        sf::st_read(path)
      }
    } else {
      stop(paste("File not found:", path))
    }
  }
  
  moose <- safe_st_read(moose_path) %>%
    st_transform(crs = crs) %>%
    select(Latitude, Longitude, name)
  
  sbfi <- safe_st_read(sbfi_path) %>%
    st_transform(crs = crs)
  
  wmu <- safe_st_read(wmu_path) %>%
    st_transform(crs = crs)
 
  transects <- safe_st_read(transects_path, layer = "tracks") %>%
    st_transform(crs = crs) %>%
    select(1)
  
  if (!is.null(transects2_path)) {
    transects2 <- safe_st_read(transects2_path, layer = "tracks") %>%
      st_transform(crs = crs) %>%
      select(1)
    transects <- rbind(transects, transects2)
  }
  
  list(moose = moose, transects = transects, sbfi = sbfi, wmu = wmu)
}


#' Split LINESTRING into Equal-Length Segments
#'
#' This function divides a LINESTRING geometry into equal-length segments. 
#' Each segment is represented as a LINESTRING.
#'
#' @param linestring A LINESTRING object from the sf package.
#' @return An sf object containing the segments as LINESTRING geometries.
#' @export
split_into_segments <- function(linestring) {
  total_length <- st_length(linestring)
  total_length <- if (!inherits(total_length, "units")) set_units(total_length, "m") else total_length
  
  num_segments <- ceiling(as.numeric(set_units(total_length, "km")))
  equal_segment_length <- set_units(total_length, "km") / num_segments
  segment_lengths <- rep(equal_segment_length, num_segments)
  segment_lengths[num_segments] <- set_units(total_length, "km") - sum(segment_lengths[1:(num_segments - 1)])
  
  points <- st_line_sample(linestring, sample = seq(0, 1, length.out = num_segments + 1)) %>%
    st_cast("POINT")
  
  segments <- map2(points[-length(points)], points[-1], ~ {
    segment <- st_sfc(st_linestring(x = st_coordinates(c(.x, .y))), crs = st_crs(linestring))
    st_sf(geometry = segment)
  })
  
  do.call(rbind, segments)
}

# Define file paths

data_paths <- list(
  "501" = list(
    moose_path = "D:\\WMU\\survey_data\\501_moose_locations.shp",
    sbfi_path = "D:\\WMU\\base_data\\CA_Forest_Satellite_Based_Inventory_2020\\clipped\\sbfi_501.shp",
    wmu_path = "D:\\WMU\\base_data\\WMU\\wmu_501_3400.shp",
    transects_path = "D:\\WMU\\survey_data\\WMU 501 (2018-2019)\\WMU501_transects_2018.gpx",
    transects2_path = NULL
  ),
  "503" = list(
    moose_path = "D:\\WMU\\survey_data\\503_moose_locations.shp",
    sbfi_path = "D:\\WMU\\base_data\\CA_Forest_Satellite_Based_Inventory_2020\\clipped\\sbfi_503.shp",
    wmu_path = "D:\\WMU\\base_data\\WMU\\wmu_503_3400.shp",
    transects_path = "D:\\WMU\\survey_data\\WMU 503 (2021-2022)\\wmu503_transects_ph1.gpx",
    transects2_path = NULL
  ),
  "512" = list(
    moose_path = "D:\\WMU\\survey_data\\512_moose_locations.shp",
    sbfi_path = "D:\\WMU\\base_data\\CA_Forest_Satellite_Based_Inventory_2020\\clipped\\sbfi_512.shp",
    wmu_path = "D:\\WMU\\base_data\\WMU\\wmu_512_3400.shp",
    transects_path = "D:\\WMU\\survey_data\\WMU 512 (2019-2020)\\WMU512_transects_EW_block1.gpx",
    transects2_path = "D:\\WMU\\survey_data\\WMU 512 (2019-2020)\\WMU512_transects_EW_block2.gpx"
  ),
  "517" = list(
    moose_path = "D:\\WMU\\survey_data\\517_moose_locations.shp",
    sbfi_path = "D:\\WMU\\base_data\\CA_Forest_Satellite_Based_Inventory_2020\\clipped\\sbfi_517.shp",
    wmu_path = "D:\\WMU\\base_data\\WMU\\wmu_517_3400.shp",
    transects_path = "D:\\WMU\\survey_data\\WMU 517 (2018-2019)\\WMU517_Transects_D4.gpx",
    transects2_path = "D:\\WMU\\survey_data\\WMU 517 (2018-2019)\\WMU517_Transects_D123.gpx"
  ),
  "528" = list(
    moose_path = "D:\\WMU\\survey_data\\528_moose_locations.shp",
    sbfi_path = "D:\\WMU\\base_data\\CA_Forest_Satellite_Based_Inventory_2020\\clipped\\sbfi_528.shp",
    wmu_path = "D:\\WMU\\base_data\\WMU\\wmu_528_3400.shp",
    transects_path = "D:\\WMU\\survey_data\\WMU 528 (2018-2019)\\WMU528_Transects_D4.gpx",
    transects2_path = "D:\\WMU\\survey_data\\WMU 528 (2018-2019)\\WMU528_Transects_D123.gpx"
  )
)

# Define WMU numbers and loop through
wmu_number_list <- c('501', '503', '512', '517', '528')


for (wmu_number in wmu_number_list){
  print(paste("Processing WMU", wmu_number))

  # Load data
  data <- load_spatial_data(
    moose_path = data_paths[[wmu_number]]$moose_path,
    sbfi_path = data_paths[[wmu_number]]$sbfi_path,
    wmu_path = data_paths[[wmu_number]]$wmu_path,
    crs = 3400,
    transects_path = data_paths[[wmu_number]]$transects_path,
    transects2_path = data_paths[[wmu_number]]$transects2_path
  )

  moose <- data$moose
  transects <- data$transects
  sbfi <- data$sbfi
  wmu <- data$wmu


  # Rename columns in sbfi for consistency
  colnames(sbfi) <- c(
    "OBJECTID", "ID", "TILE", "AREA_HA", "PERIMETER_M", "JURISDICTION",
    "ECOZONE", "ECOPROVINCE", "ECOREGION", "MANAGEMENT", "LC_WATER",
    "LC_SNOW_ICE", "LC_ROCK_RUBBLE", "LC_EXPOSED_BARREN", "LC_BRYOIDS",
    "LC_SHRUBS", "LC_WETLAND", "LC_WETLAND_TREED", "LC_HERBS",
    "LC_CONIFEROUS", "LC_BROADLEAF", "LC_MIXEDWOOD", "LC_TREED",
    "LC_FAO_FOREST", "LC_WETLAND_VEGETATION", "DISTURB_FIRE_PERC",
    "DISTURB_FIRE_YEAR", "DISTURB_FIRE_MAGNITUDE_MIN", "DISTURB_FIRE_MAGNITUDE_MAX",
    "DISTURB_FIRE_MAGNITUDE_AVG", "DISTURB_FIRE_MAGNITUDE_SD", "DISTURB_FIRE_MAGNITUDE_MEDIAN",
    "DISTURB_HARVEST_PERC", "DISTURB_HARVEST_YEAR", "RECOVERY_FIRE_MIN",
    "RECOVERY_FIRE_MAX", "RECOVERY_FIRE_AVG", "RECOVERY_FIRE_SD",
    "RECOVERY_FIRE_MEDIAN", "RECOVERY_HARVEST_MIN", "RECOVERY_HARVEST_MAX",
    "RECOVERY_HARVEST_AVG", "RECOVERY_HARVEST_SD", "RECOVERY_HARVEST_MEDIAN",
    "AGE_MIN", "AGE_MAX", "AGE_AVG", "AGE_SD", "AGE_MEDIAN",
    "AGE_0_10", "AGE_10_20", "AGE_20_30", "AGE_30_40", "AGE_40_50",
    "AGE_50_60", "AGE_60_70", "AGE_70_80", "AGE_80_90", "AGE_90_100",
    "AGE_100_110", "AGE_110_120", "AGE_120_130", "AGE_130_140",
    "AGE_140_150", "AGE_GT_150", "STRUCTURE_CANOPY_HEIGHT_MIN",
    "STRUCTURE_CANOPY_HEIGHT_MAX", "STRUCTURE_CANOPY_HEIGHT_AVG",
    "STRUCTURE_CANOPY_HEIGHT_SD", "STRUCTURE_CANOPY_HEIGHT_MEDIAN",
    "STRUCTURE_CANOPY_COVER_MIN", "STRUCTURE_CANOPY_COVER_MAX",
    "STRUCTURE_CANOPY_COVER_AVG", "STRUCTURE_CANOPY_COVER_SD",
    "STRUCTURE_CANOPY_COVER_MEDIAN", "STRUCTURE_LOREYS_HEIGHT_MIN",
    "STRUCTURE_LOREYS_HEIGHT_MAX", "STRUCTURE_LOREYS_HEIGHT_AVG",
    "STRUCTURE_LOREYS_HEIGHT_SD", "STRUCTURE_LOREYS_HEIGHT_MEDIAN",
    "STRUCTURE_BASAL_AREA_MIN", "STRUCTURE_BASAL_AREA_MAX",
    "STRUCTURE_BASAL_AREA_AVG", "STRUCTURE_BASAL_AREA_SD",
    "STRUCTURE_BASAL_AREA_MEDIAN", "STRUCTURE_BASAL_AREA_TOTAL",
    "STRUCTURE_AGB_MIN", "STRUCTURE_AGB_MAX", "STRUCTURE_AGB_AVG",
    "STRUCTURE_AGB_SD", "STRUCTURE_AGB_MEDIAN", "STRUCTURE_AGB_TOTAL",
    "STRUCTURE_VOLUME_MIN", "STRUCTURE_VOLUME_MAX", "STRUCTURE_VOLUME_AVG",
    "STRUCTURE_VOLUME_SD", "STRUCTURE_VOLUME_MEDIAN", "STRUCTURE_VOLUME_TOTAL",
    "SPECIES_NUMBER", "SPECIES_1", "SPECIES_1_PERC", "SPECIES_2",
    "SPECIES_2_PERC", "SPECIES_3", "SPECIES_3_PERC", "SPECIES_4",
    "SPECIES_4_PERC", "SPECIES_5", "SPECIES_5_PERC", "SPECIES_CONIFEROUS_PERC",
    "SPECIES_CML_1", "SPECIES_CML_1_PERC", "SPECIES_CML_2",
    "SPECIES_CML_2_PERC", "SPECIES_CML_3", "SPECIES_CML_3_PERC",
    "SPECIES_CML_4", "SPECIES_CML_4_PERC", "SPECIES_CML_5",
    "SPECIES_CML_5_PERC", "SPECIES_CML_CONIFEROUS_PERC",
    "SPECIES_CML_ASSEMBLAGES", "SPECIES_CML_ASSEMBLAGES_PERC",
    "SYMB_LAND_BASE_LEVEL", "SYMB_LAND_COVER_LEVEL", "SYMB_VEGETATION_LEVEL",
    "SYMB_DISTURBANCE", "SYMB_RECOVERY", "SYMB_AGE", "Shape_Length",
    "Shape_Area", "layer", "path", "geometry"
  )

  # Select only the OBJECTID column from wmu
  wmu <- wmu[, "OBJECTID"]

  # Perform spatial join between moose points and sbfi polygons
  #
  # This step joins the moose points with the sbfi polygons to enrich the moose data
  # with vegetation structure metrics from sbfi. Missing values are replaced with zero.
  joined <- st_join(moose, sbfi) %>%
    select(
      STRUCTURE_CANOPY_HEIGHT_MEDIAN,
      STRUCTURE_CANOPY_COVER_MEDIAN,
      STRUCTURE_AGB_MEDIAN,
      STRUCTURE_VOLUME_MEDIAN
    ) %>%
    replace_na(list(
      STRUCTURE_CANOPY_HEIGHT_MEDIAN = 0,
      STRUCTURE_CANOPY_COVER_MEDIAN = 0,
      STRUCTURE_AGB_MEDIAN = 0,
      STRUCTURE_VOLUME_MEDIAN = 0
    ))

  # Add vegetation structure metrics to the moose data
  #
  # The metrics from the joined data are added to the moose dataset, with some rounding
  # for efficiency and consistency.
  moose <- moose %>%
    mutate(
      canopy_height = round(joined$STRUCTURE_CANOPY_HEIGHT_MEDIAN, 0),
      canopy_cover = round(joined$STRUCTURE_CANOPY_COVER_MEDIAN, 0),
      agb = round(joined$STRUCTURE_AGB_MEDIAN, -1),
      vol = round(joined$STRUCTURE_VOLUME_MEDIAN, -1)
    )

  # Split transects into equal-length segments
  #
  # Transects are converted to LINESTRING geometries and split into segments using
  # the `split_into_segments` function. The segments are then combined into a single
  # sf object.
  transects_segments <- transects %>%
    st_cast("LINESTRING") %>%
    st_geometry() %>%
    map(split_into_segments) %>%
    bind_rows() %>%
    st_sf()

  # Set CRS for the segments
  st_crs(transects_segments) <- st_crs(transects)

  # Add sample labels to the segments
  transects_segments$Sample.Label <- row_number(transects_segments)

  # Calculate distances between moose points and transect segments
  #
  # Compute the distance between each moose point and each segment of the transects.
  distances <- st_distance(moose, transects_segments)

  # Identify the closest segment for each moose point
  #
  # For each moose point, find the nearest transect segment and the distance to it.
  closest_segments <- tibble(
    moose_id = seq_len(nrow(moose)),
    Sample.Label = map_int(seq_len(nrow(moose)), ~ which.min(distances[., ])),
    distance = map_dbl(seq_len(nrow(moose)), ~ min(distances[., ]))
  )

  # Merge closest segment information with the moose data
  #
  # Append the closest segment and distance information to the moose dataset.
  moose <- moose %>%
    mutate(
      object = row_number(),
      size = 1,
      Effort = 10,
      distance = closest_segments$distance / 1000
    ) %>%
    left_join(closest_segments, by = c("object" = "moose_id")) %>%
    rename(
      Transect.Label = name,
      Sample.Label = Sample.Label,
      distance = distance.x
    ) %>%
    select(-distance.y)

  # Filter out moose points with distance greater than 600 meters
  #
  # Remove moose points that are more than 600 meters from the nearest transect segment.
  moose <- moose %>%
    filter(distance <= 0.6)

  # Extract x and y coordinates from the moose sf object
  coords <- st_coordinates(moose)

  # Prepare data for density surface modelling
  #
  # Convert moose data to a dataframe for modelling, including x and y coordinates.
  segdata <- moose %>%
    select(
      Latitude, Longitude, Effort, Transect.Label, Sample.Label,
      canopy_height, canopy_cover, agb, vol
    ) %>%
    st_drop_geometry() %>%
    as.data.frame() %>%
    mutate(
      x = coords[, 1],
      y = coords[, 2]
    )

  # Prepare distance data for density surface modelling
  #
  # Create a dataframe for distance data, including x and y coordinates.
  distdata <- moose %>%
    select(
      object, Latitude, Longitude, distance, Effort, size,
      canopy_height, canopy_cover, agb, vol
    ) %>%
    st_drop_geometry() %>%
    mutate(
      detected = 1,
      x = coords[, 1],
      y = coords[, 2]
    )

  # Prepare observation data
  #
  # Create a dataframe with observation data.
  obsdata <- moose %>%
    select(object, distance, Effort, Sample.Label, size) %>%
    st_drop_geometry() %>%
    as.data.frame()

  # Save the processed data
  #
  # Save the prepared datasets to a file for further use.
  output_path <- here("Output", "PrepData", paste0("prepared_", wmu_number, ".RData"))
  save(segdata, distdata, obsdata, wmu, file = output_path)
}