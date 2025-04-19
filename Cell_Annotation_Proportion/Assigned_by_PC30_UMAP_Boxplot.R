library(Seurat)
library(magrittr)
library(dplyr)
library(Matrix)
library(tools)
library(ggplot2)
library(ggpubr)

setwd("/home/f06b22037/SSD2/JW/1136project_SingleCell/JW_customized/UMAP_recolor/PC30/all")


###################################################### Step 1. extract PC1 ~ PC30

get_object_PC30 <- function(integration_RDS) {

  combined_object <- readRDS(integration_RDS)
  combined_object_PC30 <- combined_object@reductions$pca@cell.embeddings
  return(combined_object_PC30)
}


input_path <- "/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Multi_species_analysis/"
integration_RDS <- paste0(input_path,"/all_data_rds/integration_PtrTensionWoodALL_13samples_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds")
combined_object_PC30 <- get_object_PC30(integration_RDS) 
write.table(combined_object_PC30, "integration_PtrTensionWoodALL_13samples_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base_PC30.txt", 
    sep="\t", quote=F, row.names=T, col.names=T)


###################################################### Step 2. assign colors by the nearest cells of Normal Bio1 (Tung) and plot UMAP
extract_pairsubset <- function(combined_object_PC30, plotting_csv, referencename, Sp2name) {
  # comepare Sp2 to refsp1 to find out the closet cells and its colors in refsp1
  combined_object_PC30_df <- combined_object_PC30[grepl(paste0(referencename,"_"), rownames(combined_object_PC30))|grepl(paste0(Sp2name,"_"), rownames(combined_object_PC30)),]  

  Dist_PC30 <- as.matrix(dist(combined_object_PC30_df))

  refsp1_Sp2 <- Dist_PC30[grepl(paste0(referencename,"_"), rownames(Dist_PC30)), grepl(paste0(Sp2name,"_"), colnames(Dist_PC30))]  
  refsp1_Sp2_compare_matrix <- data.frame(matrix(0, ncol(refsp1_Sp2), 2))
  colnames(refsp1_Sp2_compare_matrix) <- c("sp2_name", "refsp1_name")
  
  for (i in 1:ncol(refsp1_Sp2)){
    refsp1_Sp2_compare_matrix[i, 1] = colnames(refsp1_Sp2)[i]
    refsp1_Sp2_compare_matrix[i, 2] = rownames(refsp1_Sp2)[order(refsp1_Sp2[, i])[1]]
  }
  

  # extract UMAP of refsp1 and Sp2
  input_MS_plotting_csv <- read.csv(plotting_csv)
  refsp1_UMAP <- input_MS_plotting_csv[input_MS_plotting_csv$Sample == referencename, c("Barcode", "UMAP.1", "UMAP.2", "SS_Cluster", "SS_Color")]
  refsp1_UMAP$refsp1_name <- paste0("geneUMI_", referencename, "_", refsp1_UMAP$Barcode)
  Sp2_UMAP <- input_MS_plotting_csv[input_MS_plotting_csv$Sample == Sp2name, c("Barcode", "UMAP.1", "UMAP.2")]
  Sp2_UMAP$Barcode <- paste0("geneUMI_", Sp2name, "_", Sp2_UMAP$Barcode)
  

  # combine
  refsp1_Sp2_compare_matrix_color <- merge(refsp1_Sp2_compare_matrix, refsp1_UMAP[, c("refsp1_name", "SS_Cluster", "SS_Color")], by = "refsp1_name", all.x = TRUE)
  colnames(refsp1_Sp2_compare_matrix_color)[2] <- c("Barcode")
  Sp2_UMAP_Cluster <- merge(Sp2_UMAP, refsp1_Sp2_compare_matrix_color, by = "Barcode") 
  
  return(list(Sp2_UMAP_Cluster = Sp2_UMAP_Cluster, refsp1_UMAP = refsp1_UMAP))
}

######

combined_object_PC30 <- read.table("integration_PtrTensionWoodALL_13samples_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base_PC30.txt", sep="\t")
plotting_csv <- paste0(input_path, "all_plotting_tables_addSS/plotting_PtrTensionWoodALL_13samples_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30.csv")

###### color list of Ptr
Refsp1_Control <- read.csv("/home/f06b22037/DiskArray_f06b22037/SSD2/RK/1136project_SingleCell/results/Single_species_analysis/all_plotting_tables/plotting_TenX_Ptr.csv")
Refsp1_Control$refsp1_name <- paste0("geneUMI_TenX_Ptr_",Refsp1_Control$Barcode)

referencename <- "TenX_Ptr"
Sp2name <- c("TenX_PtrVertical01","TenX_PtrVertical02","TenX_PtrVertical03","TenX_PtrVertical04",
            "TenX_PtrTW", "TenX_PtrTension02", "TenX_PtrTension03","TenX_PtrTension04",
                "TenX_PtrOW", "TenX_PtrOpposite02", "TenX_PtrOpposite03", "TenX_PtrOpposite04")

######  purple - change to darker one
Refsp1_Control2 <- Refsp1_Control
Refsp1_Control2$Color <- ifelse(Refsp1_Control2$Color=="#9467BD", "#6D289D", Refsp1_Control2$Color) 
Refsp1_Control <- Refsp1_Control2


###### Plot UMAP with assigned colors
Refsp1_sp2_results <- list()
Refsp1_sp2_UMAP_Cluster_Color <- list()
for (i in 1:length(Sp2name)){
    results <- extract_pairsubset(combined_object_PC30, plotting_csv, referencename, Sp2name[i])
    Refsp1_sp2_results[[i]] <- results
    # recolor sp2
    results2 <- merge(results$Sp2_UMAP_Cluster[,1:4], Refsp1_Control[,c("refsp1_name","Cluster","Color")], by = "refsp1_name")
    Refsp1_sp2_UMAP_Cluster_Color[[i]] <- results2
    

    # plotting UMAP
    output_without_margin = FALSE
    output_names <- paste0("UMAP_integration_dist_assign_cluster_", Sp2name[i],".png")
    png(filename = output_names, width = 2800, height = 2000, res = 400);

    plot(
        x = Refsp1_sp2_UMAP_Cluster_Color[[i]]$UMAP.1,
        y = Refsp1_sp2_UMAP_Cluster_Color[[i]]$UMAP.2,
        col = Refsp1_sp2_UMAP_Cluster_Color[[i]]$Color,
        pch = 20,
        cex = 0.2,
        axes = !output_without_margin, las = 1,
        xlim = c(-12,12),
        ylim = c(-10,10),
        ylab = "UMAP.2",
        xlab = "UMAP.1"
        )
    dev.off()
}

Refsp1_sp2_UMAP_Cluster_Color_all <- do.call(rbind,Refsp1_sp2_UMAP_Cluster_Color)

write.table(Refsp1_sp2_UMAP_Cluster_Color_all, "integration_PtrTensionWoodALL_13samples_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base_PC30_assingedColors.csv", sep=",", quote=F, row.names=T, col.names=T)

###### re-Plot UMAP for Normal Bio1 (Tung)

Refsp1_UMAP <- Refsp1_sp2_results[[1]]$refsp1_UMAP
Refsp1_UMAP_Cluster_Color <- merge(Refsp1_UMAP[,c("refsp1_name","UMAP.1","UMAP.2")], Refsp1_Control[,c("refsp1_name","Cluster","Color")], by = "refsp1_name")


png(filename = "UMAP_integration_cluster_Ptr.png", width = 2800, height = 2000, res = 400);

plot(
    x = Refsp1_UMAP_Cluster_Color$UMAP.1,
    y = Refsp1_UMAP_Cluster_Color$UMAP.2,
    col = Refsp1_UMAP_Cluster_Color$Color,
    pch = 20,
    cex = 0.2,
    axes = !output_without_margin, las = 1,
    xlim = c(-12,12),
    ylim = c(-10,10),
    ylab = "UMAP.2",
    xlab = "UMAP.1"
    )
dev.off()

###################################################### Step 2. proportion calculation and boxplot

###### calculation
Count <- list()

Refsp1_sp2_UMAP_Cluster_Color[[13]] <- Refsp1_UMAP_Cluster_Color

Names <- c("PtrVertical01","PtrVertical02","PtrVertical03","PtrVertical04",
            "PtrTW", "PtrTension02", "PtrTension03","PtrTension04",
                "PtrOW", "PtrOpposite02", "PtrOpposite03", "PtrOpposite04","Ptr")



for (i in 1:length(Names)){

    Count[[i]] <- data.frame(table(Refsp1_sp2_UMAP_Cluster_Color[[i]]$Cluster))
    Count[[i]]$percent <- Count[[i]]$Freq/sum(Count[[i]]$Freq)*100
    Count[[i]]$Sample <- Names[i]
}



percent_output <- do.call(rbind, Count)
write.table(percent_output, "integration_PtrTensionWoodALL_13samples_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base_PC30_proportion.txt", sep="\t", quote=F, row.names=T, col.names=T)


######
library(ggpubr)
library(dplyr)

input_file <- "/Users/joweihsieh/Dropbox/YCL/tension_wood/20231115/results/integration/all/integration_PtrTensionWoodALL_13samples_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base_PC30_proportion.txt"
data <- read.table(input_file, header = TRUE)


combined_data <- data %>% filter(Var1 %in% 1:8)
ordered_var1 <- c(6, 4, 2, 1, 7, 3, 5, 8)

plot_comparison <- function(data, type1, type2, colors, filename) {
  df <- data %>%
    filter(Type %in% c(type1, type2)) %>%
    mutate(
      Var1 = factor(Var1, levels = ordered_var1),
      Type = factor(Type, levels = c(type1, type2))
    )
  
  p <- ggboxplot(df, x = "Type", y = "percent",
                 color = "Type", palette = colors,
                 ylab = "Percentage", xlab = "Type") +
    facet_wrap(~ Var1, scales = "free", ncol = 5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 46))
  
  print(p)
  
  ggsave(filename, plot = p, width = 30, height = 15, units = "cm", dpi = 300)
}

plot_comparison(
  data = combined_data,
  type1 = "PtrVertical",
  type2 = "PtrTension",
  colors = c("black", "black"),
  filename = "Boxplot_combined_Boxplot_1_to_8_vert_tens.png"
)

# Vertical vs Opposite
plot_comparison(
  data = combined_data,
  type1 = "PtrVertical",
  type2 = "PtrOpposite",
  colors = c("black", "black"),
  filename = "Boxplot_combined_Boxplot_1_to_8_vert_op.png"
)



######################################################
######################################################
######################################################  ray and fusiform - separate vertical and opposite
# === Load and preprocess data ===
library(dplyr)
library(ggpubr)

# Read data
input_file <- "/Users/joweihsieh/Dropbox/YCL/tension_wood/20231115/results/integration/all/integration_PtrTensionWoodALL_13samples_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base_PC30_proportion.txt"
data <- read.table(input_file, header = TRUE)

# Remove unwanted clusters
filtered_data <- data[!data$Var1 %in% c(9, 10), ]

# Calculate total Freq per Sample
sample_totals <- filtered_data %>%
  group_by(Sample) %>%
  summarise(total = sum(Freq))

# Add percentage to table
filtered_data$percent2 <- mapply(function(i, sample) {
  freq <- filtered_data[i, "Freq"]
  total <- sample_totals$total[sample_totals$Sample == sample]
  freq / total * 100
}, seq_len(nrow(filtered_data)), filtered_data$Sample)


# === Define lineages ===
ray_ids <- c(3, 5, 8)
fusiform_ids <- c(1, 2, 4, 6, 7)
early_ids <- c(6, 4)
late_ids <- c(1, 2, 7)

# Flexible lineage mapping function
assign_lineage <- function(var_id, lineage_type) {
  if (lineage_type == "ray_fusiform") {
    if (var_id %in% ray_ids) return("ray")
    if (var_id %in% fusiform_ids) return("fusiform")
  } else if (lineage_type == "early_late") {
    if (var_id %in% early_ids) return("early")
    if (var_id %in% late_ids) return("late")
  }
  return("others")
}

# Main lineage analysis function
analyze_lineage <- function(df, lineage_type, contrast1, contrast2, plot_filename_prefix, test_alt) {
  df$lineage <- sapply(df$Var1, assign_lineage, lineage_type = lineage_type)
  df <- df[df$lineage != "others", ]
  df$Sample2 <- paste0(df$Sample, "_", df$lineage)

  # Sum Freq by Sample2
  sum_data <- df %>% group_by(Sample2) %>% summarise(amount = sum(Freq))
  sum_all <- df %>% group_by(Sample) %>% summarise(total = sum(Freq))

  # Match Sample order (rotate if needed)
  sum_all <- sum_all[c(2:13, 1), ]

  # Calculate percentage
  sum_data$percent <- mapply(function(i, sample) {
    amt <- sum_data$amount[i]
    tot <- sum_all$total[ceiling(i/2)]
    amt / tot * 100
  }, seq_len(nrow(sum_data)), sum_data$Sample2)

  sum_data$lineage <- rep(if (lineage_type == "ray_fusiform") c("fusiform", "ray") else c("early", "late"), 13)
  sum_data$Treatment <- factor(c(rep("PtrOpposite", 8), rep("PtrTension", 8), rep("PtrVertical", 10)),
                                levels = c("PtrVertical", "PtrTension", "PtrOpposite"))

  # Plot and test
  for (contrast in list(c(contrast1, contrast2))) {
    subset_data <- sum_data[sum_data$Treatment %in% contrast, ]
    p <- ggboxplot(subset_data, x = "Treatment", y = "percent",
                   #color = "Treatment", palette = c("black", "blue", "darkgrey"),
                   order = contrast, ylab = "Percentage", xlab = "") +
      facet_wrap(~ lineage, scales = "free") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      coord_cartesian(ylim = c(0, 100))

    ggsave(paste0(plot_filename_prefix, paste(contrast, collapse = "_"), ".png"),
           plot = p, width = 15, height = 15, units = "cm", dpi = 300)

    for (lin in unique(subset_data$lineage)) {
      test_data <- subset_data[subset_data$lineage == lin, ]
      result <- pairwise.t.test(test_data$percent, test_data$Treatment, p.adj = "none", alternative = test_alt)
      print(paste("Lineage:", lin))
      print(result$p.value)
    }
  }
}

# === Run lineage analyses ===
analyze_lineage(filtered_data, "ray_fusiform", "PtrVertical", "PtrTension", "Boxplot_combined_Boxplot_ray_fusiform_", "greater")
analyze_lineage(filtered_data, "ray_fusiform", "PtrVertical", "PtrOpposite", "Boxplot_combined_Boxplot_ray_fusiform_", "greater")
analyze_lineage(filtered_data, "early_late", "PtrVertical", "PtrTension", "Boxplot_combined_Boxplot_early_late_fusiform_", "less")
analyze_lineage(filtered_data, "early_late", "PtrVertical", "PtrOpposite", "Boxplot_combined_Boxplot_early_late_fusiform_", "less")
