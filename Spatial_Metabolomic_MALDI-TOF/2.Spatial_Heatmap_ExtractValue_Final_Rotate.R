##################################################################
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
rgb_percent_rotate_output <- file.path(output_dir, paste0(Labels, "_rotate.png"))

####################################################################################################

######
######
######
######
###### legend colors

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

######
######
######
######
###### read image

#main_img_path <- "Normal_Bio2_GAci_195.png"
main_img <- load.image(main_img_path)
plot(main_img)

#x_grids <- 40
#y_grids <- 57

###### repeat image
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
  #labs(title = "RGB Colors from Center Points") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))


ggsave(rgb_percent_repeat_output, plot = p, width = 210 / 30 , height = 400 / 40, dpi = 300)


###### image rotation 

x_step <- dim(main_img)[1] / x_grids  # Width direction
y_step <- dim(main_img)[2] / y_grids  # Height direction
distance <- sqrt((x_step)^2 + (y_step)^2)
print(x_grids)
print(y_grids)
print(distance)

binary_image <- threshold(main_img, "50%") 
plot(binary_image, main = "binary_image")


coords <- as.data.frame(which(binary_image > 0, arr.ind = TRUE))
fit <- prcomp(coords) 
angle <- atan2(fit$rotation[2, 1], fit$rotation[1, 1]) * 180 / pi  


main_img_rotated <- imrotate(main_img, angle = 180 - angle)
plot(main_img_rotated)

###### signal report

binary_image <- threshold(main_img_rotated, "50%")  
plot(binary_image, main = "binary_image")

coords <- as.data.frame(which(binary_image > 0, arr.ind = TRUE))
x_range <- range(coords$dim1)  
y_range <- range(coords$dim2)  


x_signal_count <- diff(x_range) + 1 
y_signal_count <- diff(y_range) + 1 



###### image with grids

png(rgb_percent_rotate_output, width = 800, height = 600)

distance <- sqrt((x_step)^2 + (y_step)^2)

x_coords <- seq(from = x_range[1], to = x_range[2], by = distance)
y_coords <- seq(from = y_range[1], to = y_range[2], by = distance)

plot(main_img_rotated, main = "")

for (x in x_coords) {
  abline(v = x, col = "white", lwd = 1) 
}

for (y in y_coords) {
  abline(h = y, col = "white", lwd = 1) 
}


dev.off()


###### representative colors for each new grid

x_step <- distance
y_step <- distance

grid_colors <- matrix(nrow = ((x_range[2] - x_range[1])/distance), ncol = ((y_range[2] - y_range[1])/distance))


for (i in 1:nrow(grid_colors)) {
  for (j in 1:ncol(grid_colors)) {
    x_start <- x_coords[i] 
    x_end <-  x_coords[(i +1)] +1
    y_start <- y_coords[j] 
    y_end <- y_coords[(j +1)] + 1
    
    grid_region <- main_img_rotated[x_start:x_end, y_start:y_end, , drop = FALSE]

    grid_pixels <- rgb(grid_region[,,1], grid_region[,,2], grid_region[,,3], maxColorValue = 1)

    dominant_color <- names(sort(table(grid_pixels), decreasing = TRUE))[1]
    second <- names(sort(table(grid_pixels), decreasing = TRUE))[2]
    if (dominant_color == "#000000" & !is.na(second) & table(grid_pixels)[second]!= 1){
    	grid_colors[i, j] <- second
    }	else{
  	    grid_colors[i, j] <- dominant_color
    }
  }
}

grid_colors_final <- t(grid_colors)


print(paste0("n_X after rotation = " , dim(grid_colors_final)[1]))
print(paste0("n_Y after rotation = " , dim(grid_colors_final)[2]))

###### convert back to RGB

rgb_df <- data.frame(
  x = integer(0),            
  y = integer(0),            
  hex_color = character(0),  
  R = integer(0),            
  G = integer(0),            
  B = integer(0)             
)

for (i in 1:ncol(grid_colors_final)) {  
  for (j in 1:nrow(grid_colors_final)) {  
    hex_color <- grid_colors_final[j, i]  
    
    rgb_values <- col2rgb(hex_color)  
    
    rgb_df <- rbind(rgb_df, data.frame(
      x = i,                        
      y = j,                        
      hex_color = hex_color,        
      R = rgb_values[1],            
      G = rgb_values[2],            
      B = rgb_values[3]             
    ))
  }
}


###### convert RGB into percentage based on legend


color_distance <- function(rgb1, rgb2) {
  sqrt(sum((rgb1 - rgb2)^2))
}

rgb_df$percentage <- NA  

for (i in 1:nrow(rgb_df)) {
  current_rgb <- c(rgb_df$R[i], rgb_df$G[i], rgb_df$B[i])  
  
  min_distance <- Inf
  best_percentage <- NA
  
  for (j in 1:nrow(legend_colors_df_hex)) {
    legend_rgb <- c(legend_colors_df_hex$R[j], legend_colors_df_hex$G[j], legend_colors_df_hex$B[j])
    
    distance <- color_distance(current_rgb, legend_rgb)
    
    if (distance < min_distance) {
      min_distance <- distance
      best_percentage <- legend_colors_df_hex$percentage[j]
    }
  }
  
  rgb_df$percentage[i] <- best_percentage
}


write.table(rgb_df, rgb_percent_txt_output, sep = "\t", quote = F, row.names = F, col.names = T)


