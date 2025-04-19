
######

# Only features that satisfy both conditions were retained: 
# (1) their maximum intensity is at least 10 times higher than the background across all regions from all samples 
# (2) their average intensity in each region is greater than their corresponding standard deviations
# feature-based
# 12 bios

setwd("/Users/joweihsieh/Dropbox/YCL/tension_wood/spatial_metabolome/20241007/")
library(readxl)
library(ggrepel)
library(ggplot2)
library(tidyr)

############################################# Step 1. filtering

# Define the function to process metabolome data
process_metabolome_data <- function(sample_list_path, file_path, fold_cutoff, rep_n) {
  
  ########################
  # Sample List Preparation
  ########################
  sample_list <- read_excel(sample_list_path)

  # Clean sample names
  sample_list$Sample_name_clean <- gsub("LCM_Meta_", "", sample_list$`Sample name`)

  # Create Group1 and Group2 based on Condition
  sample_list$Group1 <- sapply(sample_list$Condition, function(x) {
    if (grepl("Layer", x)) {
      # Extract the first three components if "Layer" is present
      paste(strsplit(x, " ")[[1]][1:3], collapse = " ")
    } else {
      x  # Return original value if "Layer" is not present
    }
  })

  sample_list$Group2 <- sapply(sample_list$Condition, function(x) {
    if (grepl("Layer", x)) {
      # Extract the first two components if "Layer" is present
      paste(strsplit(x, " ")[[1]][1:2], collapse = " ")
    } else {
      x  # Return original value if "Layer" is not present
    }
  })

  sample_list <- sample_list[!grepl("Blank",sample_list$Sample_name_clean),]
  ########################
  # Read Data and Apply Filters
  ########################
  # Read the data file
  data <- read.csv(file_path, header = T, sep = "\t")

  modified_data <- gsub(".*_(NEG|POS)_", "\\1_", colnames(data))  
  modified_data <- gsub("\\.", "-", modified_data)      
  colnames(data) <- modified_data

  data_Neg <- data[data$mode =="negative", c(grepl("NEG", colnames(data)))]
  data_Pos <- data[data$mode =="positive", c(grepl("POS", colnames(data)))]


  modified_data_Neg <- gsub("^NEG_", "", colnames(data_Neg))
  modified_data_Pos <- gsub("^POS_", "", colnames(data_Pos))

  colnames(data_Neg) <- modified_data_Neg
  colnames(data_Pos) <- modified_data_Pos


  data_Neg_Pos <- rbind(data_Neg, data_Pos)
  data_Neg_Pos_2 <- cbind(data[,c(1:2)], data_Neg_Pos)

  # fill in zero for NA
  data_Neg_Pos_2[is.na(data_Neg_Pos_2)] <- 0

  # Create Mode_ID column
  data_Neg_Pos_2$Mode_ID <- paste0(data_Neg_Pos_2$mode, "_", data_Neg_Pos_2$mass_binned)
  colname_names <- colnames(data_Neg_Pos_2)

  # mean of blank
  data_Neg_Pos_2$Blank <- apply(data_Neg_Pos_2[,grepl("Blank", colnames(data_Neg_Pos_2))] , 1, mean)

  # Find maximum values across selected columns and filter based on fold cutoff
  max_values <- apply(data_Neg_Pos_2[, colname_names[colname_names %in% sample_list$Sample_name_clean]], 1, max)
  filter_condition <- max_values >= fold_cutoff * data_Neg_Pos_2$Blank
  data_rm <- data_Neg_Pos_2[filter_condition, ]
  

  # Select desired columns
  matching_columns <- colnames(data_rm) %in% sample_list$Sample_name_clean
  data_selected <- data_rm[, c("Mode_ID", "Blank", names(data_rm)[matching_columns])]

  ########################
  # Calculate Mean and Standard Deviation
  ########################
  data_subset <- data_selected[, 3:ncol(data_selected)]

  # Modify the column names for NB and TB
  colnames_selected <- gsub("NB(\\d+)_(L\\d+)-.*", "N\\2", colnames(data_selected)[3:ncol(data_selected)])  # For NB
  colnames_selected <- gsub("TB(\\d+)_(L\\d+)-.*", "T\\2", colnames_selected)  # For TB


  # Initialize result storage
  mean_sd_list <- list()

  # Calculate means and standard deviations
  for (i in 1:length(unique(colnames_selected))) {
    treatments <- unique(colnames_selected)[i]

    group <- data_subset[, which(colnames_selected==treatments)]

    # Calculate mean and SD
    mean_vals <- rowMeans(group, na.rm = TRUE)
    sd_vals <- apply(group, 1, sd, na.rm = TRUE)

    # Store results in list
    mean_sd_list[[paste0("Mean_", treatments)]] <- mean_vals
    mean_sd_list[[paste0("SD_", treatments)]] <- sd_vals
  }

  # Convert list to data frame
  mean_sd_results <- as.data.frame(mean_sd_list)
  rownames(mean_sd_results) <- data_selected$Mode_ID

  ########################
  # Compare Mean and SD
  ########################
  result_labels <- data.frame(Metabolite = data_selected$Mode_ID)

  # Loop through mean and SD columns in pairs
  for (i in seq(1, ncol(mean_sd_results), by = 2)) {
    mean_col <- mean_sd_results[, i]
    sd_col <- mean_sd_results[, i + 1]

    # Create new column for comparison
    column_name <- paste0("KorR_", unique(colnames_selected)[[((i + 1) / 2)]])  # Naming the new column
    print(column_name)
    result_labels[[column_name]] <- ifelse(mean_col > sd_col, "KEEP", "REMOVED")
  }

  # Determine overall status based on individual results
  result_labels$Overall_Status <- ifelse(rowSums(result_labels[, -1] == "REMOVED") > 0, "REMOVED", "KEEP")

  # Return the results
  return(list(data_selected = data_selected, mean_sd_results = mean_sd_results, result_labels = result_labels))
}


result_rep12_10 <- process_metabolome_data(
   sample_list_path = "/Users/joweihsieh/Dropbox/YCL/tension_wood/spatial_metabolome/20241004/20240923_LCM spatial metabolome sample list_9PM.xlsx",
   file_path = "/Users/joweihsieh/Dropbox/YCL/tension_wood/spatial_metabolome/20241007/ms1_combined_LCM_spatial.tsv",
   fold_cutoff = 10,
   rep_n = 12
)

############################################# Step 2. get PC

# Define the function for PCA analysis, plotting PCA scatter, and boxplot of PC1
perform_pca_analysis <- function(data_selected, mean_sd_results, result_labels, output_csv, output_path = "pca_plot.png", boxplot_path = "boxplot_PC1.png") {
  
  # Filter mean_sd_results for KEEP metabolites  
  data_selected_KEEP <- data_selected[data_selected$Mode_ID %in% result_labels[result_labels$Overall_Status == "KEEP", "Metabolite"], ]
  write.table(data_selected_KEEP, output_csv, sep = "\t", quote = F, row.names = F, col.names = T)

  
  
  # Prepare data for PCA
  data_for_pca <- data_selected_KEEP[,c(3:ncol(data_selected_KEEP))]
  extracted_strings <- gsub("NB(\\d+)_(L\\d+)-.*", "N\\2", colnames(data_for_pca))
  extracted_strings <- gsub("TB(\\d+)_(L\\d+)-.*", "T\\2", extracted_strings)

  
  # Perform PCA
  pca_result <- prcomp(t(data_for_pca), scale. = TRUE)  # Added scaling
  
  # Convert PCA results to data frame
  pca_data <- as.data.frame(pca_result$x)
  
  # Add categories to the PCA data
  pca_data$category <- extracted_strings
  
  # Define color values based on categories
  unique_categories <- unique(pca_data$category)  # Get unique categories
  color_values <- setNames(c("black", "#8D8D92", "lightgrey"), unique_categories[1:3])
  
  if (length(unique_categories) > 3) {
    additional_colors <- c("darkred", "#CB769E", "#F2D7CF")
    color_values <- c(color_values, setNames(additional_colors[1:(length(unique_categories) - 3)], unique_categories[4:length(unique_categories)]))
  }
  
  # Create the PCA scatter plot
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = category)) +
    geom_point(size = 6) +
    scale_color_manual(values = color_values) +  # Use the defined colors
    xlab(paste("PC1 - ", round(summary(pca_result)$importance[2, 1] * 100, 2), "%")) +
    ylab(paste("PC2 - ", round(summary(pca_result)$importance[2, 2] * 100, 2), "%")) +
    ggtitle(paste0 ("n = ",nrow(data_selected_KEEP))) +  # Set the plot title
    theme_minimal() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      plot.title = element_text(size = 20, hjust = 0.5),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
  
  # Display the PCA scatter plot
  print(pca_plot)
  
  # Save the PCA scatter plot
  ggsave(output_path, plot = pca_plot, width = 10, height = 8, dpi = 300)
  
  # Create a boxplot for PC1 based on category
  boxplot_PC1 <- ggplot(pca_data, aes(x = category, y = PC1, fill = category)) +
    geom_boxplot() +
    scale_fill_manual(values = color_values) +  # Use the same color scheme as the PCA scatter plot
    xlab("Category") +
    ylab("PC1") +
    ggtitle("Boxplot of PC1 by Category") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      plot.title = element_text(size = 20, hjust = 0.5),
      legend.position = "none"  # Hide legend for the boxplot
    )
  
  # Display the boxplot
  print(boxplot_PC1)
  
  # Save the boxplot
  ggsave(boxplot_path, plot = boxplot_PC1, width = 10, height = 8, dpi = 300)
  
  NL_group <- pca_data[grepl("NL", pca_data$category), , drop = FALSE]
  TL_group <- pca_data[grepl("TL", pca_data$category), , drop = FALSE]
  
  
  return(pca_data)
}


pca_data_rep12_10 <- perform_pca_analysis(result_rep12_10[[1]], result_rep12_10[[2]], result_rep12_10[[3]], output_csv = "LCM_metabolite_12rep_features_10x_0418.txt", output_path = "pca_plot_rep12_10x_0418.png", boxplot_path = "boxplot_PC1_rep12_10x_0418.png")


############################################# Step 3. distance calculation
######### dist 10 2 vs 3 (based on PC)


pc_data <- pca_data_rep12_10[, grepl("PC", colnames(pca_data_rep12_10))]

euclidean_distances <- data.frame(as.matrix(dist(pc_data, method = "euclidean")))
colnames(euclidean_distances) <- rownames(euclidean_distances)
dist_string <- gsub("NB(\\d+)_(L\\d+)-.*", "N\\2", colnames(euclidean_distances))
dist_string <- gsub("TB(\\d+)_(L\\d+)-.*", "T\\2", dist_string)

  

euclidean_distances_NL1_NL2 <- euclidean_distances[which(dist_string=="NL1"), which(dist_string=="NL2")]
euclidean_distances_NL1_NL2$mean_dist <- rowMeans(euclidean_distances_NL1_NL2)
euclidean_distances_NL1_NL2$category <- paste0("NL1_NL2")


euclidean_distances_NL2_NL3 <- euclidean_distances[which(dist_string=="NL2"), which(dist_string=="NL3")]
euclidean_distances_NL2_NL3$mean_dist <- rowMeans(euclidean_distances_NL2_NL3)
euclidean_distances_NL2_NL3$category <- paste0("NL2_NL3")


euclidean_distances_TL1_TL2 <- euclidean_distances[which(dist_string=="TL1"), which(dist_string=="TL2")]
euclidean_distances_TL1_TL2$mean_dist <- rowMeans(euclidean_distances_TL1_TL2)
euclidean_distances_TL1_TL2$category <- paste0("TL1_TL2")


euclidean_distances_TL2_TL3 <- euclidean_distances[which(dist_string=="TL2"), which(dist_string=="TL3")]
euclidean_distances_TL2_TL3$mean_dist <- rowMeans(euclidean_distances_TL2_TL3)
euclidean_distances_TL2_TL3$category <- paste0("TL2_TL3")


euclidean_distances_NL1_NL2$ID  <- rownames(euclidean_distances_NL1_NL2)
euclidean_distances_NL2_NL3$ID  <- rownames(euclidean_distances_NL2_NL3)

euclidean_distances_TL1_TL2$ID  <- rownames(euclidean_distances_TL1_TL2)
euclidean_distances_TL2_TL3$ID  <- rownames(euclidean_distances_TL2_TL3)


euclidean_distances_NL1_NL2_NL3 <- rbind(euclidean_distances_NL1_NL2[,c("ID","mean_dist", "category")], euclidean_distances_NL2_NL3[,c("ID","mean_dist", "category")])
euclidean_distances_TL1_TL2_TL3 <- rbind(euclidean_distances_TL1_TL2[,c("ID","mean_dist", "category")], euclidean_distances_TL2_TL3[,c("ID","mean_dist", "category")])

euclidean_distances_NL123_TL123 <- rbind(euclidean_distances_NL1_NL2_NL3, euclidean_distances_TL1_TL2_TL3)


t.test(euclidean_distances_NL1_NL2$mean_dist, euclidean_distances_NL2_NL3$mean_dist)$p.value
t.test(euclidean_distances_TL1_TL2$mean_dist, euclidean_distances_TL2_TL3$mean_dist)$p.value

t.test(euclidean_distances_NL1_NL2$mean_dist, euclidean_distances_TL1_TL2$mean_dist)$p.value
t.test(euclidean_distances_NL2_NL3$mean_dist, euclidean_distances_TL2_TL3$mean_dist)$p.value
  


euclidean_distances_NL123_TL123$category <- factor(euclidean_distances_NL123_TL123$category, 
                                                   levels = c("NL1_NL2", "TL1_TL2", "NL2_NL3", "TL2_TL3"))

write.table(euclidean_distances_NL123_TL123, "boxplot_distance_rep12_10x_12_23_values.txt", sep = "\t", quote = F, row.names = F, col.names = T)


############################################# Step 4. plotting

# Create the boxplot
ggplot(euclidean_distances_NL123_TL123, aes(x = category, y = mean_dist)) +  # Use fill instead of color
    geom_boxplot(color = "black") +  # Optionally set the outline color
    #scale_fill_manual(values = color_values) +  # Use the defined fill colors
    labs(title = "",
         x = "",
         y = "Distance between two layers") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      plot.title = element_text(size = 20, hjust = 0.5),
      legend.position = "none"  # Hide legend for the boxplot
    )

ggsave("boxplot_distance_rep12_10x_12_23.png", width = 10, height = 8, dpi = 300)


