library(imager)
library(ggplot2)
library(dplyr)
##################################################################

# Read arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
working_dir <- args[1]
output_dir <- args[2]
main_img_path <- args[3]
x_grids <- as.integer(args[4])
y_grids <- as.integer(args[5])
Labels <- args[6]

# Set working directory
setwd(working_dir)

# Define legend image path
legend_img_path <- "legend.png"

# Define output filenames
legend_colors_output <- file.path(output_dir, paste0(Labels, "_legend_colors_repeat.png"))
rgb_percent_repeat_output <- file.path(output_dir, paste0(Labels, "_repeat.png"))
rgb_percent_txt_output <- file.path(output_dir, paste0(Labels, "_rgb_percent.txt"))


##################################################################
# Extract colors and their corresponding values from the legend
legend_img <- load.image(legend_img_path)
plot(legend_img)

legend_array <- as.array(legend_img)

# Image size and middle
legend_height <- dim(legend_array)[2]  # Height
legend_width <- dim(legend_array)[1]   # Width
middle_column <- round(legend_width / 2)

# Cut height into 100 indices
n_points <- 101
y_indices <- round(seq(1, legend_height, length.out = n_points))  # Sample indices in height range

# Extract normalized RGB and multiply by 255
legend_colors_df <- data.frame(matrix(0, 101, 3))

for (i in 1:100) {
  legend_colors_df[i, ] <- legend_array[middle_column, y_indices[i], 1, 1:3] * 255 
}

colnames(legend_colors_df) <- c("R", "G", "B")

# Convert RGB to hex colors
legend_rgb <- apply(legend_colors_df, 1, function(color) {
  rgb(color[1], color[2], color[3], maxColorValue = 255)  # RGB range is 0-255
})

# Plot and check extracted legend colors
png(legend_colors_output, width = 1200, height = 300)
barplot(rep(1, 101), col = legend_rgb, border = NA,
        main = "Extracted Colors from Legend", 
        axes = FALSE, space = 0, ylim = c(0, 1))
dev.off()

# Generate a table for color and percentage
legend_colors_df_hex <- cbind(legend_colors_df, legend_rgb)
legend_colors_df_hex$percentage <- (rev(seq(0:100)) - 1)

# Extract colors from each spot in the main image
main_img <- load.image(main_img_path)
plot(main_img)

# Calculate the grid size
x_step <- dim(main_img)[1] / x_grids  # Width direction
y_step <- dim(main_img)[2] / y_grids  # Height direction

# Create an empty data frame to store RGB values of each grid
rgb_values_center <- data.frame(x_grid = integer(0), y_grid = integer(0), R = numeric(0), G = numeric(0), B = numeric(0))

for (i in 1:x_grids) {
  for (j in 1:y_grids) {
    x_min <- round(1 + (i - 1) * x_step)
    x_max <- round(i * x_step)
    y_min <- round(1 + (j - 1) * y_step)
    y_max <- round(j * y_step)

    x_center <- round((x_min + x_max) / 2)
    y_center <- round((y_min + y_max) / 2)
    
    center_rgb <- main_img[x_center, y_center, 1:3] * 255  # Multiply by 255 to convert to 0-255 RGB
    rgb_values_center <- rbind(rgb_values_center, data.frame(x_grid = i, y_grid = j, R = center_rgb[1], G = center_rgb[2], B = center_rgb[3]))
  }
}

# Convert RGB values to colors
rgb_values_center$color <- rgb(rgb_values_center$R / 255, rgb_values_center$G / 255, rgb_values_center$B / 255)

# Plot RGB colors as tiles
p <- ggplot(rgb_values_center, aes(x = x_grid, y = -1 * y_grid, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  theme_minimal() +
  labs(title = "RGB Colors from Center Points") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

ggsave(rgb_percent_repeat_output, plot = p, width = 600 / 40, height = 140 / 30, dpi = 300)

# Function to convert RGB to HEX
rgb_to_hex <- function(R, G, B) {
  sprintf("#%02X%02X%02X", R, G, B)
}

# Create HEX values for each RGB value in the center points
rgb_values_center$hex <- mapply(rgb_to_hex, rgb_values_center$R, rgb_values_center$G, rgb_values_center$B)

# Calculate the Euclidean distance between RGB values and legend's RGB values
calc_distance <- function(rgb1, rgb2) {
  sqrt(sum((rgb1 - rgb2)^2))
}

# Function to assign percentage based on closest matching HEX
assign_percentage <- function(hex_color, legend_df) {
  legend_rgb <- legend_df$legend_rgb
  legend_percentage <- legend_df$percentage
  
  distances <- sapply(legend_rgb, function(legend_hex) {
    legend_rgb_values <- col2rgb(legend_hex)[, 1]
    input_rgb_values <- col2rgb(hex_color)[, 1]
    calc_distance(input_rgb_values, legend_rgb_values)
  })
  
  closest_match_index <- which.min(distances)
  return(legend_percentage[closest_match_index])
}

# Apply the percentage assignment based on the closest legend color
rgb_values_center$percentage <- sapply(rgb_values_center$hex, assign_percentage, legend_df = legend_colors_df_hex)

write.table(rgb_values_center, rgb_percent_txt_output, sep = "\t", quote = F, row.names = F, col.names = T)

