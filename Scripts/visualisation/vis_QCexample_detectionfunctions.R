# Load necessary library
library(grDevices)
library(RColorBrewer)
library(here)

# Function to compute detectability values
compute_detectability <- function(key.function, scale.param, shape.param, truncation) {
  distances <- seq(0, truncation, length.out = 100)
  if (key.function == "hn") {
    detectability <- exp(-distances^2 / (2 * scale.param^2))
  } else if (key.function == "uf") {
    detectability <- rep(scale.param, length(distances))
  } else {
    stop("Unsupported key function")
  }
  data.frame(distance = distances, detectability = detectability)
}

# Generate data for each detection function
detect_H_data <- compute_detectability("hn", 140, NULL, 500)
detect_FW2_data <- compute_detectability("hn", 170, NULL, 110)
detect_FWG_data <- compute_detectability("hn", 170, NULL, 155)
detect_uf_data <- compute_detectability("uf", 1, NULL, 40)

# Define custom colors
custom_colors <- c("orange", "#1f78b4")

# Set the output file path using the here package
output_path <- here("Output", "Plots", "quadcopterDetectability_plot.tiff")

# Open a TIFF device with square dimensions
tiff(filename = output_path, width = 1500, height = 1500, units = "px", res = 300)

# Set graphical parameters for A4 size
par(mar = c(5, 5, 3, 2) + 0.1, cex = 0.9, cex.lab = 1.2, cex.axis = 1, lwd = 2, bty = "l", family = "Arial")

# Plot the data
plot(detect_FWG_data$distance, detect_FWG_data$detectability, type = "l", col = custom_colors[2], lwd = 5,
     xlim = c(0, 170), ylim = c(0, 1), xlab = "Distance (m)", ylab = "Detection probability")

# Plot the additional line
lines(detect_uf_data$distance, detect_uf_data$detectability, col = custom_colors[1], lwd = 5)

# Add a legend
legend("bottomright", legend = c("NADIR","Gimbal"), col = custom_colors, lty = c(1, 1), lwd = 4,
       cex = 1, box.lwd = 1, box.col = "NA", bg = "white")

# Close the TIFF device
dev.off()
