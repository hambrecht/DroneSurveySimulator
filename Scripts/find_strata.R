# Load the terra package
library(terra)

# Step 1: Create two test rasters of size 1000x1000 with random values between 0 and 1
# set.seed(123)  # For reproducibility
r1 <- rast(nrows = 1000, ncols = 1000, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
r2 <- rast(nrows = 1000, ncols = 1000, xmin = 0, xmax = 100, ymin = 0, ymax = 100)

# Assign random values between 0 and 1 to the rasters
values(r1) <- runif(ncell(r1), 0, 1)
values(r2) <- runif(ncell(r2), 0, 1)

# Step 2: Multiply the two rasters
combined <- r1 * r2
plot(combined)

# Step 3: Apply a focal smoothing (mean over a 3x3 window)
combined_smoothed99 <- focal(combined, w = matrix(1, 99, 99), fun = mean)
plot(combined_smoothed99)
combined_smoothed99 <- focal(combined_smoothed99, w = matrix(1, 99, 99), fun = mean)
plot(combined_smoothed99)

combined_smoothed9 <- focal(combined, w = matrix(1, 9, 9), fun = mean)
plot(combined_smoothed9)
combined_smoothed9 <- focal(combined_smoothed9, w = matrix(1, 9, 9), fun = mean)
plot(combined_smoothed9)


# Step 4: Dynamically define the breaks based on the value range of combined_smoothed
min_value <- min(values(combined_smoothed99), na.rm = TRUE) # Get the minimum value
max_value <- max(values(combined_smoothed99), na.rm = TRUE) # Get the maximum value

# Define breaks (e.g., 4 equal intervals across the value range)
breaks <- seq(min_value, max_value, length.out = 5) # Creates 4 categories

# Define categories (e.g., 1, 2, 3, 4 for the 4 intervals)
categories <- 1:(length(breaks) - 1)

combined_categorized <- classify(combined_smoothed99, rcl = cbind(breaks[-length(breaks)], breaks[-1], categories))
plot(combined_categorized)

# Step 5: Group contiguous areas (patches) for each category separately
patches_list <- lapply(categories, function(cat) {
  cat_raster <- combined_categorized == cat
  cat_patches <- patches(cat_raster, directions = 8, zeroAsNA = TRUE)
  cat_patches <- mask(cat_patches, cat_raster) # Mask to retain only the patches of the current category
  cat_patches
})

# Combine the patches from all categories
combined_patches <- do.call(merge, patches_list)
plot(combined_patches)
# head(combined_patches[])
# value <- terra::extract(combined_patches, cbind(x, y))




# Step 6: Convert the patches raster to polygons and preserve the value attribute
patches_polygons <- as.polygons(combined_patches, dissolve = FALSE)
# patches_polygons$value <- combined_patches[]

# Plot the patches polygons
plot(patches_polygons)

# Dissolve polygons by category to create one homogeneous area per category
homogenous_areas <- aggregate(patches_polygons, by = "value")

# Step 7: Optionally, merge small patches with larger neighboring patches to reduce the number of patches
# Define a threshold for the minimum patch size (e.g., 1000 square units)
min_patch_size <- 1000
large_patches <- homogenous_areas[area(homogenous_areas) > min_patch_size, ]
small_patches <- homogenous_areas[area(homogenous_areas) <= min_patch_size, ]

# Merge small patches with neighboring large patches
for (i in 1:nrow(small_patches)) {
  small_patch <- small_patches[i, ]
  neighbors <- st_touches(small_patch, large_patches, sparse = FALSE)
  if (any(neighbors)) {
    large_patch_index <- which(neighbors)[1]
    large_patches[large_patch_index, ] <- st_union(large_patches[large_patch_index, ], small_patch)
  }
}

# Combine the large patches back into a single object
reduced_patches <- rbind(large_patches, small_patches[!st_touches(small_patches, large_patches, sparse = FALSE), ])

# Step 8: Plot the results
plot(reduced_patches, col = c("lightblue", "lightgreen", "orange", "red"))

# Optionally, save the output as a file (e.g., GeoPackage format)
# writeVector(reduced_patches, "reduced_patches.gpkg", overwrite=TRUE)

