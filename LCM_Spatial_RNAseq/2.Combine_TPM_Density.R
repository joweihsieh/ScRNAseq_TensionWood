######## combine
sample_names <- c("Normal_Bio1_L1", "Normal_Bio1_L2", "Normal_Bio1_L3", 
    "Normal_Bio2_L1", "Normal_Bio2_L2", "Normal_Bio2_L3",
    "Normal_Bio3_L1", "Normal_Bio3_L2", "Normal_Bio3_L3",
     "Tension_Bio1_L1", "Tension_Bio1_L2", "Tension_Bio1_L3",
  	"Tension_Bio2_L1", "Tension_Bio2_L2", "Tension_Bio2_L3",
  	"Tension_Bio3_L1", "Tension_Bio3_L2", "Tension_Bio3_L3")


######## combine count table
file_paths1 <- list.files("/home/woodydrylab/FileShare/LCM_spatial_RNAseq/TNNR24PRR02378_2xRNAseq_lib/01.Raw_read/LCMSpTrNW", pattern = "*counts.tsv", full.names = TRUE)
file_paths2 <- list.files("/home/woodydrylab/FileShare/LCM_spatial_RNAseq/TNNR24PRR02378_2xRNAseq_lib/01.Raw_read/LCMSpTrTW", pattern = "*counts.tsv", full.names = TRUE)

file_paths <- c(file_paths1, file_paths2)

combined_data <- NULL

for (file_path in file_paths) {
  data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sample_name <- gsub("_counts.tsv$", "", basename(file_path))
  colnames(data)[2] <- sample_name
  
  if (is.null(combined_data)) {
    combined_data <- data
  } else {
    combined_data <- merge(combined_data, data, by = "gene", all = TRUE) 
  }
}

combined_data[is.na(combined_data)] <- 0

output_path <- "/home/woodydrylab/FileShare/LCM_spatial_RNAseq/TNNR24PRR02378_2xRNAseq_lib/01.Raw_read/count_table.txt"
write.table(combined_data, output_path, sep = "\t", row.names = FALSE, quote = FALSE)

######## normalization into TPM

library(tidyr)
library(ggplot2)
library(DESeq2)
library(rtracklayer)
library(dplyr)
library(tibble)  # Load tibble package

#### get gene length
# Path to the GTF file
gtf_file <- "/home/woodydrylab/FileShare/LCM_spatial_RNAseq/TNNR24PRR02378_2xRNAseq_lib/01.Raw_read/genome/Ptrichocarpa_533_v4.1.gene.gtf"
gtf_data <- import(gtf_file, format = "gtf")
genes <- gtf_data[gtf_data$type == "transcript"]

gene_info <- as.data.frame(genes) %>%
  select(seqnames, gene_id = gene_id, start, end)

gene_info <- gene_info %>%
  mutate(length = end - start + 1)


gene_info <- gene_info %>%
  group_by(gene_id) %>%
  slice_max(length, with_ties = FALSE) %>% 
  ungroup() 

write.csv(gene_info, "/home/woodydrylab/FileShare/LCM_spatial_RNAseq/TNNR24PRR02378_2xRNAseq_lib/01.Raw_read/genome/gene_lengths.csv", row.names = FALSE)


#
rownames(combined_data) <- combined_data$gene
combined_data <- combined_data[, -1]  
raw_counts <- combined_data

raw_counts <- raw_counts %>%
  mutate(gene_id = rownames(raw_counts)) %>%
  left_join(gene_info, by = "gene_id") %>%
  select(-seqnames, -start, -end) 


raw_counts <- na.omit(raw_counts)

calculate_tpm <- function(counts, lengths) {
  rpk <- counts / (lengths / 1000)  
  scale_factor <- colSums(rpk) / 1e6  
  tpm <- t(t(rpk) / scale_factor) 
  return(tpm)
}

counts <- raw_counts[, 1:ncol(combined_data)]  
gene_lengths <- raw_counts$length

tpm_matrix <- calculate_tpm(as.matrix(counts), gene_lengths)
rownames(tpm_matrix) <- raw_counts$gene_id
write.table(tpm_matrix, "/home/woodydrylab/FileShare/LCM_spatial_RNAseq/TNNR24PRR02378_2xRNAseq_lib/01.Raw_read/tpm_table.txt", row.names = TRUE, col.names = T, quote = F, sep = "\t")

######## how many genes are dectected
non_zero_counts <- apply(tpm_matrix, 2, function(x) sum(x != 0))
all_zero_rows <- apply(tpm_matrix, 1, function(x) all(x == 0))

if (any(all_zero_rows)) {
  print("Some rows contain only zeros across all columns:")
  print(rownames(tpm_matrix)[all_zero_rows])
} else {
  print("No rows contain only zeros across all columns.")
}

######## filtering

tpm_df <- data.frame(tpm_matrix)  

Treatment_Layers <- c("NL1", "NL2", "NL3", 
    "NL1", "NL2", "NL3",
    "NL1", "NL2", "NL3",
    "TL1", "TL2", "TL3",
  	"TL1", "TL2", "TL3",
  	"TL1", "TL2", "TL3")

Treatment_Layers_type <- unique(Treatment_Layers)

for (i in 1:length(Treatment_Layers_type)) {
  mean_col <- paste0("mean_", Treatment_Layers_type[i])
  SD_col <- paste0("SD_", Treatment_Layers_type[i])
  keep_col <- paste0("Keep_", Treatment_Layers_type[i])  
  
  sub <- tpm_df[, which(Treatment_Layers == Treatment_Layers_type[i])]
  tpm_df[[mean_col]] <- apply(sub, 1, function(x) mean(x, na.rm = TRUE))
  tpm_df[[SD_col]] <- apply(sub, 1, function(x) sd(x, na.rm = TRUE))
  
  tpm_df[[keep_col]] <- ifelse(tpm_df[[mean_col]] > tpm_df[[SD_col]], "Keep", "Discard")
}

keep_columns <- grep("^Keep_", colnames(tpm_df), value = TRUE)  
tpm_df$Keep_All <- apply(tpm_df[, keep_columns], 1, function(x) all(x == "Keep"))  

filtered_tpm_df <- tpm_df[tpm_df$Keep_All == TRUE, 1:length(Treatment_Layers)]
#325



######## SD distribution

sd_columns <- tpm_df[, grep("SD", colnames(tpm_df))]
filtered_sd_20_tmp_df <- tpm_df[apply(sd_columns, 1, function(x) all(x <= 20)), 1:length(Treatment_Layers)]

########


data_for_pca <- tpm_matrix
#11294
Treatment_Layers <- c("NL1", "NL2", "NL3", 
    "NL1", "NL2", "NL3",
    "NL1", "NL2", "NL3",
    "TL1", "TL2", "TL3",
  	"TL1", "TL2", "TL3",
  	"TL1", "TL2", "TL3")


dist_string <- Treatment_Layers

perform_pca_analysis <- function(data_for_pca, Treatment_Layers, output_path, boxplot_path) {
  
  
  pca_result <- prcomp(t(data_for_pca), scale. = TRUE)  # Added scaling
  pca_data <- as.data.frame(pca_result$x)
  
  # Add categories to the PCA data
  pca_data$category <- Treatment_Layers
  
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
    #ggtitle(paste0 ("n = ",nrow(data_selected_KEEP))) +  # Set the plot title
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
  
  
  # Perform t-test between NL and TL groups as a whole
  combined_NL <- NL_group$PC1
  combined_TL <- TL_group$PC1
  t_test_NL_TL <- t.test(combined_NL, combined_TL)
  cat("T-test result between NL group and TL group:\n")
  print(t_test_NL_TL$p.value)
  
  return(pca_data)
}


pca_data_SD_20 <- perform_pca_analysis(filtered_sd_20_tmp_df, Treatment_Layers, output_path = "pca_plot_SD_20.png", boxplot_path = "boxplot_PC1_SD_20.png")



######### density plot

compute_density_plot <- function(pca_data, group_labels, selected_keys, output_value, output_file_density) {

  # Step 1: Prepare the PCA data (select only PC columns)
  pc_data_clean <- pca_data[, grepl("PC", colnames(pca_data))]

  # Step 2: Calculate Euclidean distances
  euclidean_distances <- as.matrix(dist(pc_data_clean, method = "euclidean"))


  # Step 4: Subset distances for specified groups
  distance_list <- list(
    NL1_NL2 = euclidean_distances[which(group_labels == "NL1"), which(group_labels == "NL2")],
    NL2_NL3 = euclidean_distances[which(group_labels == "NL2"), which(group_labels == "NL3")],
    TL1_TL2 = euclidean_distances[which(group_labels == "TL1"), which(group_labels == "TL2")],
    TL2_TL3 = euclidean_distances[which(group_labels == "TL2"), which(group_labels == "TL3")]
  )


  selected_distances <- distance_list[selected_keys]
  #print(selected_distances)
  print(selected_keys[2])
  ks_test_result <- ks.test(selected_distances[[selected_keys[1]]], selected_distances[[selected_keys[2]]])
  print(ks_test_result)

  t_test_result <- t.test(selected_distances[[selected_keys[1]]], selected_distances[[selected_keys[2]]])
  print(t_test_result)

  # Step 2: Combine selected distances into a single data frame
  combined_df <- do.call(rbind, lapply(names(selected_distances), function(category) {
    data.frame(distance = as.vector(selected_distances[[category]]), category = category)
  }))

  write.csv(combined_df, output_value, row.names = FALSE, quote = F)

  hist_data_list <- lapply(selected_distances, function(dist) {
    # Generate histogram with breaks (bins)
    hist(dist, plot = FALSE)
  })
  

  # Step 5: Combine histogram data into a data frame
  combined_hist_df <- do.call(rbind, lapply(names(hist_data_list), function(category) {
    hist_data <- hist_data_list[[category]]
    data.frame(
      category = category,
      count = hist_data$counts,
      mid = hist_data$mids
    )
  }))


 # Step 6: Create a density plot
  p_density <- ggplot(combined_df, aes(x = distance, fill = category, color = category)) +
    geom_density(alpha = 0.4) +  # Adjust transparency
    labs(title = "",
         x = "Euclidean Distance",
         y = "Density") +
    scale_fill_manual(values = c("grey", "orange")) +
    scale_color_manual(values = c("grey", "orange")) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 15),
      plot.title = element_text(size = 20, hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      legend.position = "right",
      panel.grid = element_blank()  # Remove gridlines
  ) +
  coord_cartesian(ylim = c(0, 0.05)) +  # Replace <MAX_Y_VALUE> with the desired maximum y-axis value
  scale_x_continuous(limits = c(0, 160), breaks = seq(0, 160, by = 10))  # Customize limits and breaks


  # Save the density plot to a file
  ggsave(output_file_density, plot = p_density, width = 12, height = 8, dpi = 300)

  return(combined_df)
}


df_density_NT12 <- compute_density_plot(pca_data_SD_20, dist_string, c("NL1_NL2", "TL1_TL2"), "values_for_density_12_SD_20.txt" , "density_distance_12_SD_20.png")
df_density_NT23 <- compute_density_plot(pca_data_SD_20, dist_string, c("NL2_NL3", "TL2_TL3"), "values_for_density_23_SD_20.txt" , "density_distance_23_SD_20.png")

