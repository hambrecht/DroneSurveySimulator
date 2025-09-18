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
detect_H_data <- compute_detectability("hn", 140, NULL, 600)
detect_FW2_data <- compute_detectability("hn", 170, NULL, 180)
detect_FWG_data <- compute_detectability("hn", 170, NULL, 260)
detect_uf_data <- compute_detectability("uf", 1, NULL, 40)
# detect_uf_data <- compute_detectability("hn", 170, NULL, 50)

# Define custom colors
custom_colors <- c("#a6cee3", "#1f78b4", "#1f78b4", "#b2df8a")

# Set the output file path using the here package
output_path <- here("Output", "Plots", "detectability_plot2.tiff")

# Open a TIFF device
tiff(filename = output_path, width = 2250, height = 1125, units = "px", res = 300)

# Set graphical parameters for A4 size
par(mar = c(5, 5, 4, 2) + 0.1, cex = 0.9, cex.lab = 0.9, cex.axis = 0.9, lwd = 1, bty = "n", family = "Arial")

# Plot the data
# plot(detect_H_data$distance, detect_H_data$detectability, type = "l", col = custom_colors[4], lwd = 7, xlim = c(0, 600), ylim = c(0, 1), xlab = "Distance", ylab = "Detectability")
plot(detect_FWG_data$distance, detect_FWG_data$detectability, type = "l", col = custom_colors[3], lwd = 7, xlim = c(0, 260), ylim = c(0, 1), xlab = "Distance", ylab = "Detectability")
# lines(detect_FWG_data$distance, detect_FWG_data$detectability, col = custom_colors[3], lwd = 7, lty = 2) # Dotted line
# lines(detect_FW2_data$distance, detect_FW2_data$detectability, col = 'black', lwd = 7) # Shadow
# lines(detect_FW2_data$distance, detect_FW2_data$detectability, col = custom_colors[2], lwd = 5)
lines(detect_uf_data$distance, detect_uf_data$detectability, col = custom_colors[1], lwd = 7)

# Add a legend
# legend("topright", legend = c("NADIR", "2 Cameras", "Gimbal", "Helicopter"), col = custom_colors, lty = c(1, 1, 2, 1), lwd = 2, cex = 0.9, box.lwd = 1, box.col = "NA", bg = "white")
legend("topright", legend = c("NADIR", "Gimbal"), col = custom_colors, lty = c(1, 1), lwd = 2, cex = 0.9, box.lwd = 1, box.col = "NA", bg = "white")


# Close the TIFF device
dev.off()
