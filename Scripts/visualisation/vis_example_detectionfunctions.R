library(RColorBrewer)
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
detect_uf_data <- compute_detectability("uf", 1, NULL, 50)
# detect_uf_data <- compute_detectability("hn", 170, NULL, 50)

# Define custom colors
custom_colors <- c("#1f78b4", "#1f78b4", "#a6cee3", "#b2df8a")

# Set graphical parameters for A4 size
par(mar = c(5, 5, 4, 2) + 0.1, cex = 1.2, cex.lab = 1.5, cex.axis = 1.2, lwd = 2, bty = "n")

# Plot the data
plot(detect_FW2_data$distance, detect_FW2_data$detectability, type = "l", col = custom_colors[1], lwd = 3, xlim = c(0, 600), ylim = c(0, 1), xlab = "Distance", ylab = "Detectability")
lines(detect_FWG_data$distance, detect_FWG_data$detectability, col = custom_colors[2], lwd = 3, lty = 2) # Dotted line
lines(detect_uf_data$distance, detect_uf_data$detectability, col = custom_colors[3], lwd = 3)
lines(detect_H_data$distance, detect_H_data$detectability, col = custom_colors[4], lwd = 3)

# Add a legend
legend("topright", legend = c("2 Cameras", "Gimbal", "NADIR", "Helicopter"), col = custom_colors, lty = c(1, 1, 2, 1), lwd = 3, cex = 1.5, box.lwd = 2, box.col = "NA", bg = "white")



