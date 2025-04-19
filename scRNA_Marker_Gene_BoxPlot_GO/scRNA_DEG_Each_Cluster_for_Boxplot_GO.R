
####### A gene is defined as a marker gene if it meets the following criteria:

#It must be a marker in either the tension or normal condition (i.e., showing either the highest or lowest expression compared to other clusters within that condition).
#The fold change compared to other clusters within the same condition must be ≥ 2.
#The corresponding cell cluster must show a p-value < 0.05 when comparing tension vs. normal.

####### 

setwd("/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/all_UMI_tables")

library(scales)
library(reshape)
library(ggplot2)
library(tidyverse)


files <- list.files(pattern = paste0("ray_fusiform.txt$"))
files <- files[c(1,11:13,10,5,2:4,9,6:8)]
Bios_labels <- c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio1", "Bio2","Bio3", "Bio4", "Bio1", "Bio2","Bio3", "Bio4")

combined_df <- data.frame()

condition <- c("VerticalBio1",
"VerticalBio2",
"VerticalBio3",
"VerticalBio4",
"VerticalBio5",
"OppositeBio1",
"OppositeBio2",
"OppositeBio3",
"OppositeBio4",
"TensionBio1",
"TensionBio2",
"TensionBio3",
"TensionBio4")

type <- c(rep("Vertical", 5), rep("Opposite", 4), rep("Tension", 4))

for (i in 1:length(files)) {
  data <- read.table(files[i], header = TRUE)  
  data$condition <- condition[i]
  data$type <- type[i]
  combined_df <- rbind(combined_df, data)
}

#write.table(combined_df, "/home/woodydrylab/FileShare/marker_gene_tension/plots_most/ray_fusiform_Tension_Normal.txt", sep = "\t", quote = F, row.names = F, col.names = T)



gene_ID <- unique(combined_df$gene_id)
names <- gene_ID


########
########

Bio_all_df_mean <- list()
all_results <- data.frame(gene_id = character(),
                          color = character(),
                          Vertical_mean = numeric(),
                          Tension_mean = numeric(),
                          p_value = numeric(),
                          stringsAsFactors = FALSE)

colors <- c("red", "orange", "yellow", "pink", "purple", "green", "brown", "blue")


fold_cutoff <- 2 #upregulation
fold_cutoff2 <- 0.5 #downregulation

########
for (i in seq_along(gene_ID)) {
  Bio_all_df_mean[[i]] <- combined_df[combined_df$gene_id %in% gene_ID[i], ]

  print(names[i])

  for (color in colors) {
    color_column <- paste0(color, "_mean")

    if (color_column %in% colnames(Bio_all_df_mean[[i]])) {
      p_values <- pairwise.t.test(
        Bio_all_df_mean[[i]][[color_column]], 
        Bio_all_df_mean[[i]]$type, 
        p.adj = "none", 
        alternative = "two.sided"
      )$p.value

      p_value <- p_values["Vertical", "Tension"]

      # Calculate means for other color columns excluding the current color
      vertical_other_color <- Bio_all_df_mean[[i]][Bio_all_df_mean[[i]]$type == "Vertical", 
                                                    !colnames(Bio_all_df_mean[[i]]) %in% color_column]
      vertical_other_color_mean_cols <- grep("mean", colnames(vertical_other_color), value = TRUE)
      vertical_other_color_mean_values <- colMeans(vertical_other_color[, vertical_other_color_mean_cols], na.rm = TRUE)
      vertical_other_color_overall_mean <- mean(vertical_other_color_mean_values, na.rm = TRUE)

      tension_other_color <- Bio_all_df_mean[[i]][Bio_all_df_mean[[i]]$type == "Tension", 
                                                   !colnames(Bio_all_df_mean[[i]]) %in% color_column]
      tension_other_color_mean_cols <- grep("mean", colnames(tension_other_color), value = TRUE)
      tension_other_color_mean_values <- colMeans(tension_other_color[, tension_other_color_mean_cols], na.rm = TRUE)
      tension_other_color_overall_mean <- mean(tension_other_color_mean_values, na.rm = TRUE)

      vertical_mean <- mean(Bio_all_df_mean[[i]][[color_column]][Bio_all_df_mean[[i]]$type == "Vertical"], na.rm = TRUE)
      tension_mean <- mean(Bio_all_df_mean[[i]][[color_column]][Bio_all_df_mean[[i]]$type == "Tension"], na.rm = TRUE)

      vertical_fold_to_other <- (vertical_mean + 1e-5) / (vertical_other_color_overall_mean + 1e-5)
      tension_fold_to_other <- (tension_mean + 1e-5) / (tension_other_color_overall_mean + 1e-5)

      # Determine marker status
      marker_status <- if (
        !is.na(p_value) && p_value < 0.05 &&
        (
          (vertical_fold_to_other >= fold_cutoff && all(vertical_mean > vertical_other_color_mean_values)) ||
          (tension_fold_to_other >= fold_cutoff && all(tension_mean > tension_other_color_mean_values)) ||
          (vertical_fold_to_other <= fold_cutoff2 && all(vertical_mean < vertical_other_color_mean_values)) ||
          (tension_fold_to_other <= fold_cutoff2 && all(tension_mean < tension_other_color_mean_values))
        )
      ) {
        "Marker"
      } else {
        "non-Marker"
      }

      # Store the results
      all_results <- rbind(all_results, data.frame(
        gene_id = gene_ID[i],
        color = color,
        Vertical_mean = vertical_mean,
        Vertical_other_mean = vertical_other_color_overall_mean,
        Vertical_fold_to_other = vertical_fold_to_other,
        Tension_mean = tension_mean,
        Tension_other_mean = tension_other_color_overall_mean,
        Tension_fold_to_other = tension_fold_to_other,
        p_value = ifelse(is.na(p_value), NA, p_value),
        MarkerGene = marker_status
      ))

      # Plot significant markers
      if (marker_status == "Marker") {
        print(paste("P-value:", p_value))

        color_dir <- file.path("/home/woodydrylab/FileShare/marker_gene_tension/plots_most", color)
        if (!dir.exists(color_dir)) dir.create(color_dir, recursive = TRUE)

        plot_name <- file.path(color_dir, paste0("box_ray_fusiform_exp_", names[i], "_", color, "_20241227.png"))

        long_df <- reshape2::melt(Bio_all_df_mean[[i]], id.vars = c("gene_id", "type", "condition"))
        long_df$value <- as.numeric(long_df$value)
        long_df$type2 <- paste0(long_df$type, "_", long_df$variable)

        long_df <- long_df[!grepl("Opposite", long_df$type2), ]
        long_df$type2 <- factor(long_df$type2, levels = c(
          "Vertical_orange_mean", "Vertical_yellow_mean", "Vertical_pink_mean", "Vertical_purple_mean", 
          "Vertical_green_mean", "Vertical_brown_mean", "Vertical_blue_mean", "Vertical_red_mean", 
          "Tension_orange_mean", "Tension_yellow_mean", "Tension_pink_mean", "Tension_purple_mean", 
          "Tension_green_mean", "Tension_brown_mean", "Tension_blue_mean", "Tension_red_mean"
        ))

        plot_object <- ggplot(long_df, aes(x = type2, y = value, fill = variable)) +
          geom_boxplot(width = 0.5) +
          geom_point(color = "black", size = 3) +
          labs(
            title = paste(names[i], "-", color, "\n", paste("p-value:", round(p_value, 4))),
            x = "Conditions",
            y = "Expression level"
          ) +
          scale_fill_manual(values = c(
            "#FF7F0F", "#F8E71C", "#E377C2", "#6D289D", 
            "#2AA02A", "#8C564B", "#1F77B4", "#D62728"
          )) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line()
          )

        ggsave(filename = plot_name, plot = plot_object, width = 30, height = 15, units = "cm", dpi = 400)
      } else {
        print(paste("P-value not significant or conditions not met for", names[i], "and color", color, "- skipping plot."))
      }
    } else {
      print(paste("Color column", color_column, "not found for gene", gene_ID[i]))
    }
  }
}
###########################

all_results_mk <- all_results[all_results$MarkerGene == "Marker", ]
color_gene_counts <- table(all_results_mk$color)

write.table(all_results_mk, "/home/woodydrylab/FileShare/marker_gene_tension/plots_most/marker_gene_list.txt", sep = "\t", quote = F, row.names = F, col.names = T)


########################### name conversion (v4 -> v3)

mk <- read.csv("/Users/joweihsieh/Dropbox/YCL/tension_wood/scRNA_marker_gene/marker_gene_list.txt", sep = "\t", header = T)
df <- read.csv("/Users/joweihsieh/Dropbox/YCL/tension_wood/Module1_GO/Ptrichocarpa_533_v4.1.synonym_20241227.txt", sep = "\t", header = T)

mk$v4_1 <- sub("\\.v4.1", "", mk$gene_id)
mk_re <- merge(mk, df, by = "v4_1", all.x = T)

write.table(mk_re, "/Users/joweihsieh/Dropbox/YCL/tension_wood/scRNA_marker_gene/marker_gene_list_names.txt", sep = "\t", quote = F, row.names = F, col.names = T)

####################################################################################################################################### co-expression network
########################### collecting individual and mean values for each marker

reshaped_data_color <- list()

for (j in 1:length(colors)) {
  subset_color <- all_results_mk[all_results_mk$color == colors[j], ]
  orders <- which(gene_ID %in% subset_color$gene_id)
  color_column <- paste0(colors[j], "_mean")
  
  reshaped_data_list <- list()
  
  for (k in seq_along(orders)) {
    subset_data <- Bio_all_df_mean[[orders[k]]] %>%
      filter(type %in% c("Vertical", "Tension")) %>%
      select(gene_id, condition, color_column, type)
    
    reshaped_data <- subset_data %>%
      select(gene_id, condition, color_column) %>%
      pivot_wider(names_from = condition, values_from = color_column) %>%
      rowwise() %>%
      mutate(
        Vertical_mean = mean(c_across(starts_with("VerticalBio")), na.rm = TRUE),
        Tension_mean = mean(c_across(starts_with("TensionBio")), na.rm = TRUE)
      )
    
    reshaped_data_list[[k]] <- reshaped_data
  }
  
  reshaped_data_color[[j]] <- do.call(rbind, reshaped_data_list)
  reshaped_data_color[[j]]$Type <- colors[j]
}

reshaped_data_color_all <- do.call(rbind, reshaped_data_color)

write.table(reshaped_data_color_all, "/home/woodydrylab/FileShare/marker_gene_tension/plots_most/marker_gene_list_indi_values.txt", sep = "\t", quote = F, row.names = F, col.names = T)


