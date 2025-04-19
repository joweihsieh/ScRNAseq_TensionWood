setwd("/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/all_UMI_tables")

library(fields)
library(ggplot2)
library(magrittr)
library(dplyr)

library(scales)
library(reshape)
library(ggplot2)


compute_cluster_means <- function(umi_file, sample_prefix, integ_clu, output_file, excluded_clusters = c(9, 10)) {
  message("Reading: ", umi_file)
  umi_data <- read.csv(umi_file)
  umi_data$Barcode <- paste0(sample_prefix, "_", umi_data$Barcode)
  
  # 合併整合的 cluster 資訊
  merged <- merge(umi_data, integ_clu, by = "Barcode")
  filtered <- merged[!merged$Cluster %in% excluded_clusters, ]
  
  filtered[, 2:(ncol(filtered) - 1)] <- lapply(filtered[, 2:(ncol(filtered) - 1)], function(x) as.numeric(as.character(x)))
  
  gene_ids <- colnames(filtered)[2:(ncol(filtered) - 1)]
  cluster_ids <- c(3, 5, 8, 6, 4, 2, 1, 7)
  cluster_names <- c("orange", "yellow", "pink", "purple", "green", "brown", "blue", "red")
  
  df <- data.frame(gene_id = gene_ids)
  
  for (i in seq_along(cluster_ids)) {
    cluster <- cluster_ids[i]
    cname <- paste0(cluster_names[i], "_mean")
    means <- sapply(2:(ncol(filtered) - 1), function(gene_col) {
      mean(filtered[filtered$Cluster == cluster, gene_col], na.rm = TRUE)
    })
    df[[cname]] <- means
  }
  
  write.table(df, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
}

integ_clu_path = "/home/woodydrylab/DiskArray/f06b22037/SSD2/JW/1136project_SingleCell/results/JW_customized/UMAP_recolor/PC30/all/integration_PtrTensionWoodALL_13samples_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base_PC30_assingedColors.csv"
Integ_Clu <- read.csv(integ_clu_path)

compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrVertical01.csv",
  sample_prefix = "geneUMI_TenX_PtrVertical01",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrVertical01_ray_fusiform.txt"
)

compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrTension02.csv",
  sample_prefix = "geneUMI_TenX_PtrTension02",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrTension02_ray_fusiform.txt"
)

compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrOpposite02.csv",
  sample_prefix = "geneUMI_TenX_PtrOpposite02",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrOpposite02_ray_fusiform.txt"
)


compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrVertical02.csv",
  sample_prefix = "geneUMI_TenX_PtrVertical02",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrVertical02_ray_fusiform.txt"
)

compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrOW.csv",
  sample_prefix = "geneUMI_TenX_PtrOW",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrOW_ray_fusiform.txt"
)

compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrTW.csv",
  sample_prefix = "geneUMI_TenX_PtrTW",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrTW_ray_fusiform.txt"
)

compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrVertical03.csv",
  sample_prefix = "geneUMI_TenX_PtrVertical03",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrVertical03_ray_fusiform.txt"
)

compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrTension03.csv",
  sample_prefix = "geneUMI_TenX_PtrTension03",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrTension03_ray_fusiform.txt"
)

compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrOpposite03.csv",
  sample_prefix = "geneUMI_TenX_PtrOpposite03",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrOpposite03_ray_fusiform.txt"
)

compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrVertical04.csv",
  sample_prefix = "geneUMI_TenX_PtrVertical04",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrVertical04_ray_fusiform.txt"
)

compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrTension04.csv",
  sample_prefix = "geneUMI_TenX_PtrTension04",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrTension04_ray_fusiform.txt"
)

compute_cluster_means(
  umi_file = "geneUMI_TenX_PtrOpposite04.csv",
  sample_prefix = "geneUMI_TenX_PtrOpposite04",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_PtrOpposite04_ray_fusiform.txt"
)


compute_cluster_means(
  umi_file = "geneUMI_TenX_Ptr.csv",
  sample_prefix = "geneUMI_TenX_Ptr",
  integ_clu = Integ_Clu,
  output_file = "geneUMI_TenX_Ptr_ray_fusiform.txt"
)


################################################################################
################################################################################
################################################################################ boxplot


files <- list.files(pattern = paste0("ray_fusiform.txt$"))
files <- files[c(1,11:13,10,5,2:4,9,6:8)]
Bios_labels <- c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio1", "Bio2", "Bio3", "Bio4", "Bio1", "Bio2", "Bio3", "Bio4")
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

################ tension wood marker genes

names <- c("PtFLA6_12", "PtrLBD39", "PtrLBD22", "PtrCOB3", "PtrCOB11", 
  "PopFLA1", "PopFLA2", "PopFLA3", "PopFLA4", "PopFLA5", "PopFLA6", "PopFLA7",
  "PopFLA8","PopFLA9","PopFLA10","PopFLA11", "PopFLA12","PopFLA13","PopFLA14","PopFLA15",
  "PHBMT1", "PtNCED3", "DDE2", "PtBES1",
  "MYB52", "PtrCOMT2", "PtACO6", "GASA14")

gene_ID <- c("Potri.013G151366.v4.1","Potri.005G097800.v4.1", "Potri.007G066700.v4.1", "Potri.004G117200.v4.1", "Potri.015G060100.v4.1",
  "Potri.019G121200.v4.1", "Potri.013G151432.v4.1", "Potri.013G151500.v4.1", "Potri.013G014200.v4.1", "Potri.019G121100.v4.1", "Potri.013G151300.v4.1", "Potri.012G014966.v4.1",
  "Potri.009G012200.v4.1","Potri.004G210600.v4.1","Potri.009G012100.v4.1","Potri.016G088700.v4.1", "Potri.014G162900.v4.1","Potri.006G129200.v4.1","Potri.001G320800.v4.1","Potri.015G129400.v4.1",
  "Potri.001G448000.v4.1","Potri.001G393800.v4.1", "Potri.002G130700.v4.1","Potri.014G041600.v4.1",
  "Potri.007G134500.v4.1","Potri.012G006400.v4.1", "Potri.004G003000.v4.1", "Potri.012G076700.v4.1")



Bio_all_df_mean <- list()
Fibers <- list()
for (i in 1:length(gene_ID)){
  name <- paste0("box_ray_fusiform_exp_", names[i], "_20240717.allbio.png")

  Bio_all_df_mean[[i]] = combined_df[combined_df$gene_id%in%gene_ID[i],]
  long_df <- melt(Bio_all_df_mean[[i]], id.vars = c("gene_id","type","condition"))
  long_df$value <- as.numeric(long_df$value)
  long_df$type2 <- paste0(long_df$type,"_",long_df$variable)

  long_df$type2 <- factor(long_df$type2, levels = c("Vertical_orange_mean", "Vertical_yellow_mean", "Vertical_pink_mean", "Vertical_purple_mean", "Vertical_green_mean", "Vertical_brown_mean", "Vertical_blue_mean", "Vertical_red_mean", 
                  "Opposite_orange_mean", "Opposite_yellow_mean", "Opposite_pink_mean", "Opposite_purple_mean", "Opposite_green_mean", "Opposite_brown_mean", "Opposite_blue_mean", "Opposite_red_mean", 
                  "Tension_orange_mean", "Tension_yellow_mean", "Tension_pink_mean", "Tension_purple_mean", "Tension_green_mean", "Tension_brown_mean", "Tension_blue_mean", "Tension_red_mean"))
 
  plot_object <- ggplot(long_df, aes(x = type2, y = value, fill = variable)) +
  geom_boxplot(width = 0.5) +  
  geom_point(color = "black", size = 3) +  
  labs(title = names[i],
       x = "Conditions",
       y = "Expression level") +  
  scale_fill_manual(values = c("#FF7F0F","#F8E71C","#E377C2","#6D289D","#2AA02A","#8C564B", "#1F77B4", "#D62728")) + 
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

  print(plot_object)
  ggsave(filename = name,
       plot = plot_object,
       width = 30, height = 15, units = "cm", dpi = 400)


  Fibers[[i]]  <- Bio_all_df_mean[[i]] 
  print(names[i])
   
  print(pairwise.t.test(Fibers[[i]]$red_mean, Fibers[[i]]$type, p.adj = "none", alternative = "two.sided"))

}

Bio_all_df_mean_all <- do.call(rbind, Bio_all_df_mean)

write.table(Bio_all_df_mean_all, "ray_fusiform_TW_markers_original_exp_level_20240717.txt", quote = F, sep="\t", row.names = F, col.names = T)

################ SND and VND


gene_ID <- 
c("Potri.015G127400.v4.1", "Potri.012G126500.v4.1","Potri.003G113000.v4.1", "Potri.001G120000.v4.1", "Potri.007G014400.v4.1", "Potri.005G116800.v4.1", 
  "Potri.011G153300.v4.1", "Potri.001G448400.v4.1", "Potri.014G104800.v4.1", "Potri.002G178700.v4.1")

#"Potri.010G170000.v4.1"
names <- c("VND6-A1", "VND6-A2", "VND6-B1", "VND6-B2", "VND6-C1", "VND6-C2",
  "SND1-A1", "SND1-A2", "SND1-B1", "SND1-B2")

# exclude ray
fusiform <- combined_df[,-c(2:4)]


Bio_all_df_mean <- list()
Fibers <- list()

for (i in 1:length(gene_ID)){
  name <- paste0("box_ray_fusiform_exp_", names[i], "_20250418.allbio.png")

  Bio_all_df_mean[[i]] = fusiform[fusiform$gene_id %in% gene_ID[i],]
  long_df <- melt(Bio_all_df_mean[[i]], id.vars = c("gene_id","type","condition"))
  long_df$value <- as.numeric(long_df$value)
  long_df$type2 <- paste0(long_df$type,"_",long_df$variable)

  long_df$type2 <- factor(long_df$type2, levels = c("Vertical_orange_mean", "Vertical_yellow_mean", "Vertical_pink_mean", "Vertical_purple_mean", "Vertical_green_mean", "Vertical_brown_mean", "Vertical_blue_mean", "Vertical_red_mean", 
                  "Opposite_orange_mean", "Opposite_yellow_mean", "Opposite_pink_mean", "Opposite_purple_mean", "Opposite_green_mean", "Opposite_brown_mean", "Opposite_blue_mean", "Opposite_red_mean", 
                  "Tension_orange_mean", "Tension_yellow_mean", "Tension_pink_mean", "Tension_purple_mean", "Tension_green_mean", "Tension_brown_mean", "Tension_blue_mean", "Tension_red_mean"))
 
  plot_object <- ggplot(long_df, aes(x = type2, y = value, fill = variable)) +
  geom_boxplot(width = 0.5) +  
  geom_point(color = "black", size = 3) +  
  labs(title = names[i],
       x = "Conditions",
       y = "Expression level") +  
  scale_fill_manual(values = c("#6D289D","#2AA02A","#8C564B", "#1F77B4", "#D62728")) + 
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

  print(plot_object)
  ggsave(filename = name,
       plot = plot_object,
       width = 30, height = 15, units = "cm", dpi = 400)


  Fibers[[i]]  <- Bio_all_df_mean[[i]] 
  print(names[i])
   
  print(pairwise.t.test(Fibers[[i]]$red_mean, Fibers[[i]]$type, p.adj = "none", alternative = "two.sided"))

}

Bio_all_df_mean_all <- do.call(rbind, Bio_all_df_mean)


