setwd("/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/all_UMI_tables")

library(magrittr)
library(dplyr)
library(reshape2)
library(readxl)

library(ggplot2)
library(RColorBrewer)
library(gplots)
library(tidyr)
library(class)
library(cluster)
library(dplyr)

################################################################## Step 1. load fold change of UMI between brown and either blue or red
Bio1 <- read.table("geneUMI_TenX_PtrVertical01_PtrTension02_PtrOpposite02_brown_blue_red.txt", header = T)
Bio2 <- read.table("geneUMI_TenX_PtrVertical02_PtrTW_PtrOW_brown_blue_red.txt", header = T)
Bio3 <- read.table("geneUMI_TenX_PtrVertical03_PtrTension03_PtrOpposite03_brown_blue_red.txt", header = T)
Bio4 <- read.table("geneUMI_TenX_PtrVertical04_PtrTension04_PtrOpposite04_brown_blue_red.txt", header = T)
Bio5 <- read.table("geneUMI_TenX_Ptr_brown_blue_red.txt", header = T)
rownames(Bio5) <- Bio5$gene_id

Bio5 <- Bio5[,c("blue_brown_fc","red_brown_fc")]
colnames(Bio5) <- c("blue_brown_fc_vertical","red_brown_fc_vertical")

################################################################## Step 2. load WGCNA output - pattern filtering

table2=read.table("WGCNA_results/WGCNA_modules_mean_log_expression_value_colors_gene_id.txt",sep="\t",header=T)
colnames(table2) <- c("row_names","blue_brown_fc_vertical","red_brown_fc_vertical","blue_brown_fc_tension","red_brown_fc_tension","blue_brown_fc_opposite","red_brown_fc_opposite","color_clu","gene_id")
table2 <- table2[,c("row_names","blue_brown_fc_vertical","blue_brown_fc_opposite","blue_brown_fc_tension","red_brown_fc_vertical","red_brown_fc_opposite","red_brown_fc_tension","color_clu","gene_id")]

colors2 = "blue"

i = colors2

condition_data <- list()
pairwise_result<- list()

for (i in colors2){
  module <- i
  data2_ori <- table2[table2$color_clu == module, c(1:7)]
  name <- paste0("box_WGCNA_mean_module_log_expression_", i, "_con_20240422.allbio.png")

  #Bio_all_df = matrix(0, 5, 6)
  Bio1_df <- data.frame(Bio1[rownames(Bio1) %in% data2_ori$row_names,])
  Bio2_df <- data.frame(Bio2[rownames(Bio2) %in% data2_ori$row_names,])
  Bio3_df <- data.frame(Bio3[rownames(Bio3) %in% data2_ori$row_names,])
  Bio4_df <- data.frame(Bio4[rownames(Bio4) %in% data2_ori$row_names,])
  Bio5_df <- data.frame(Bio5[rownames(Bio5) %in% data2_ori$row_names,])

  condition <- colnames(Bio1_df)
  colnames(Bio1_df) <- paste0(colnames(Bio1_df), "_bio1")
  colnames(Bio2_df) <- paste0(colnames(Bio2_df), "_bio2")
  colnames(Bio3_df) <- paste0(colnames(Bio3_df), "_bio3")
  colnames(Bio4_df) <- paste0(colnames(Bio4_df), "_bio4")
  colnames(Bio5_df) <- paste0(colnames(Bio5_df), "_bio5")

  Bio1_df$condition <- rownames(Bio1_df)
  Bio2_df$condition <- rownames(Bio2_df)
  Bio3_df$condition <- rownames(Bio3_df)
  Bio4_df$condition <- rownames(Bio4_df)
  Bio5_df$condition <- rownames(Bio5_df)

  Bio_df_12 <- merge(Bio1_df, Bio2_df,  by = "condition")
  Bio_df_123 <- merge(Bio_df_12, Bio3_df,  by = "condition")
  Bio_df_1234 <- merge(Bio_df_123, Bio4_df,  by = "condition")
  Bio_df_all <- merge(Bio_df_1234, Bio5_df,  by = "condition")
  rownames(Bio_df_all) <- Bio_df_all$condition
  for (j in 1:length(condition)){
  	condition_data[[j]] <- Bio_df_all[,grepl(condition[j], colnames(Bio_df_all))]
  	condition_data[[j]][,paste0(condition[j],"_mean")] <- apply(condition_data[[j]], 1 , function(x) mean(x))
  }

	Bio_df_all_mean <- do.call(cbind, condition_data)
	
  # fiber: tension > op and normal?
  Bio_df_all_mean$tension_en <- ifelse(Bio_df_all_mean$red_brown_fc_tension_mean > Bio_df_all_mean$red_brown_fc_opposite_mean & 
    Bio_df_all_mean$red_brown_fc_tension_mean > Bio_df_all_mean$red_brown_fc_vertical_mean, "yes" ,"no")
	
  # fiber changes in opposite > vessel changes?
  Bio_df_all_mean$vessel_low <- ifelse(Bio_df_all_mean$red_brown_fc_opposite_mean > Bio_df_all_mean$blue_brown_fc_tension_mean & 
    Bio_df_all_mean$red_brown_fc_opposite_mean > Bio_df_all_mean$blue_brown_fc_vertical_mean & 
    Bio_df_all_mean$red_brown_fc_opposite_mean > Bio_df_all_mean$blue_brown_fc_opposite_mean, "yes" ,"no")

  # fiber changes in vertical > vessel changes?
	Bio_df_all_mean$vessel_low2 <- ifelse(Bio_df_all_mean$red_brown_fc_vertical_mean > Bio_df_all_mean$blue_brown_fc_tension_mean & 
    Bio_df_all_mean$red_brown_fc_vertical_mean > Bio_df_all_mean$blue_brown_fc_vertical_mean & 
    Bio_df_all_mean$red_brown_fc_vertical_mean > Bio_df_all_mean$blue_brown_fc_opposite_mean, 
    "yes" ,"no")




	test <- Bio_df_all_mean[(Bio_df_all_mean$tension_en == "yes" & Bio_df_all_mean$vessel_low == "yes"  
    & Bio_df_all_mean$vessel_low2 == "yes"), grepl("mean",colnames(Bio_df_all_mean))]

  # higer than brown in vertical and opposite (those lower than brown are excluded)
	test2 <- test[test$red_brown_fc_vertical_mean > 0 & test$red_brown_fc_opposite_mean > 0,]
  # difference of changes in red/brown between tension and veritcal > between vertical and opposite
  test2_ft <- test2[(test2$red_brown_fc_tension_mean - test2$red_brown_fc_vertical_mean) > (test2$red_brown_fc_vertical_mean - test2$red_brown_fc_opposite_mean),]
	test2 <- test2_ft


  ### extract 
	Bio_df_all_mean_ft <- Bio_df_all_mean[rownames(Bio_df_all_mean) %in% rownames(test2),]
	Bio_df_all_mean_ft$Gene.ID <- rownames(Bio_df_all_mean_ft)

	test2_mean <- data.frame(Expression = apply(test2, 2 , function(x) mean(x)))
	test2$Gene.ID <- rownames(test2)
	long_df <- melt(test2, id.vars = "Gene.ID")
  long_df$variable <- factor(long_df$variable, levels = c("blue_brown_fc_vertical_mean", "blue_brown_fc_opposite_mean", "blue_brown_fc_tension_mean", "red_brown_fc_vertical_mean", "red_brown_fc_opposite_mean", "red_brown_fc_tension_mean"))



  #pairwise_result[[i]] <- pairwise.t.test(long_df$value, long_df$variable, p.adj = "none", alternative = "greater")
 
  plot_object <- ggplot(long_df, aes(x = variable, y = value)) +
    geom_boxplot(fill = "white", width = 0.5) +  
    geom_point(color = "black", size = 3) +  
    labs(title = paste0("selection in Module1: #", nrow(test2)),
         x = "Conditions",
         y = "Expression changes") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      panel.grid.major = element_line(color = "gray", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.line = element_line()
    )

  print(plot_object)
  ggsave(filename = name,
         plot = plot_object,
         width = 25, height = 15, units = "cm", dpi = 300)
}

################################################################## Step 3. t-test - 160 genes retained

CompleteTable <- read_excel("Ptrichocarpa_v4.1_complete table_with best hit in Egr, Lch_v1.0.xlsx")
CompleteTable2 <- data.frame(CompleteTable[,c(1,7,8)])


vertical <- c(7:11)
tension <- c(18:21)
opposite <- c(28:31)

Bio_df_all_mean_ft2 = Bio_df_all_mean_ft
Bio_df_all_mean_ft$v_t_p = apply(Bio_df_all_mean_ft[, c(vertical, tension)], 1, function(y){t.test(y[1:5], y[6:10], alternative="less")$p.value}) #sig
Bio_df_all_mean_ft$t_o_p = apply(Bio_df_all_mean_ft[, c(tension, opposite)], 1, function(y){t.test(y[1:4], y[5:9], alternative="greater")$p.value}) #sig
Bio_df_all_mean_ft$v_o_p = apply(Bio_df_all_mean_ft[, c(vertical, opposite)], 1, function(y){t.test(y[1:5], y[6:10], alternative="two.sided")$p.value}) #non
Bio_df_all_mean_ft$b_r_p = apply(Bio_df_all_mean_ft[, c(1:5,13:16, 23:26, 7:11, 18:21,28:31)], 1, function(y){t.test(y[1:13], y[14:26], alternative="less")$p.value}) #sig

Bio_df_all_mean_ft_anno <- merge(Bio_df_all_mean_ft, CompleteTable2, by = "Gene.ID", all.x = T)
write.table(Bio_df_all_mean_ft_anno, "WGCNA_modules_mean_indi_blue_filter_gene_id_annotation_20240501_p.txt", quote = F, sep = "\t", row.names =F, col.names = T)


# filter pvalue

Bio_df_all_mean_ft_anno_ttest <- Bio_df_all_mean_ft_anno[Bio_df_all_mean_ft_anno$v_t_p < 0.05 & Bio_df_all_mean_ft_anno$t_o_p < 0.05 & Bio_df_all_mean_ft_anno$v_o_p > 0.05 & Bio_df_all_mean_ft_anno$b_r_p < 0.05,]


################################################################## rename and write files


Table <- Bio_df_all_mean_ft_anno_ttest

colnames(Table)[1:33] <- c("Gene.ID", "blue_brown_fc_vertical_Bio5", "blue_brown_fc_vertical_Bio2",
"blue_brown_fc_vertical_Bio3", "blue_brown_fc_vertical_Bio4", "blue_brown_fc_vertical_Bio1",
"blue_brown_fc_vertical_mean", "red_brown_fc_vertical_Bio5",  "red_brown_fc_vertical_Bio2", 
"red_brown_fc_vertical_Bio3",  "red_brown_fc_vertical_Bio4",  "red_brown_fc_vertical_Bio1", 
"red_brown_fc_vertical_mean",  "blue_brown_fc_tension_Bio2",  "blue_brown_fc_tension_Bio1", 
"blue_brown_fc_tension_Bio3",  "blue_brown_fc_tension_Bio4",  "blue_brown_fc_tension_mean", 
"red_brown_fc_tension_Bio2",   "red_brown_fc_tension_Bio1",   "red_brown_fc_tension_Bio3",  
"red_brown_fc_tension_Bio4",   "red_brown_fc_tension_mean",   "blue_brown_fc_opposite_Bio2",
"blue_brown_fc_opposite_Bio1", "blue_brown_fc_opposite_Bio3", "blue_brown_fc_opposite_Bio4",
"blue_brown_fc_opposite_mean", "red_brown_fc_opposite_Bio2",  "red_brown_fc_opposite_Bio1", 
"red_brown_fc_opposite_Bio3",  "red_brown_fc_opposite_Bio4",  "red_brown_fc_opposite_mean")


Final_Table <- Table[,c("Gene.ID","blue_brown_fc_vertical_Bio1","blue_brown_fc_vertical_Bio2", "blue_brown_fc_vertical_Bio3", "blue_brown_fc_vertical_Bio4", "blue_brown_fc_vertical_Bio5", "blue_brown_fc_vertical_mean",
    "blue_brown_fc_opposite_Bio1","blue_brown_fc_opposite_Bio2", "blue_brown_fc_opposite_Bio3", "blue_brown_fc_opposite_Bio4","blue_brown_fc_opposite_mean",
    "blue_brown_fc_tension_Bio1","blue_brown_fc_tension_Bio2", "blue_brown_fc_tension_Bio3", "blue_brown_fc_tension_Bio4", "blue_brown_fc_tension_mean",
    "red_brown_fc_vertical_Bio1","red_brown_fc_vertical_Bio2", "red_brown_fc_vertical_Bio3", "red_brown_fc_vertical_Bio4", "red_brown_fc_vertical_Bio5", "red_brown_fc_vertical_mean",
    "red_brown_fc_opposite_Bio1","red_brown_fc_opposite_Bio2", "red_brown_fc_opposite_Bio3", "red_brown_fc_opposite_Bio4", "red_brown_fc_opposite_mean",
    "red_brown_fc_tension_Bio1","red_brown_fc_tension_Bio2", "red_brown_fc_tension_Bio3", "red_brown_fc_tension_Bio4", "red_brown_fc_tension_mean", 
    "b_r_p","v_t_p", "t_o_p", "v_o_p","Best.hit.Arabidopsis", "Annotated.function")]

write.table(Final_Table, "WGCNA_modules_mean_indi_blue_filter_gene_id_annotation_20240501_p_005_160genes.txt", quote = F, sep = "\t", row.names =F, col.names = T)



################################################################## boxplot for 160 genes
#================== PLOT: 160 genes mean log2FC for each bio

rm_cols <- grep("mean", colnames(Final_Table), value = TRUE)
df3 <- Final_Table[, !colnames(Final_Table) %in% rm_cols]
df3 <- df3[, -1]  # remove gene_id

# Define groupings
group_indices <- list(
  blue_brown_fc_vertical = 1:5,
  blue_brown_fc_opposite = 6:9,
  blue_brown_fc_tension  = 10:13,
  red_brown_fc_vertical  = 14:18,
  red_brown_fc_opposite  = 19:22,
  red_brown_fc_tension   = 23:26
)

# Compute mean across each bio group
df3_all <- do.call(rbind, lapply(names(group_indices), function(group) {
  df <- data.frame(Expression = colMeans(df3[, group_indices[[group]], drop = FALSE]))
  df$type <- group
  return(df)
}))

df3_all$type <- factor(df3_all$type, levels = names(group_indices))

# Pairwise t-test
df3_all_pairwise_result <- pairwise.t.test(df3_all$Expression, df3_all$type, 
                                           p.adj = "none", alternative = "greater")

# Plotting
plot_object2 <- ggplot(df3_all, aes(x = type, y = Expression, fill = type)) +
  geom_boxplot(width = 0.5) +
  geom_point(color = "black", size = 3) +
  labs(
    title = paste0("Selection in Module1: #", nrow(Final_Table)),
    x = "Conditions",
    y = "Expression changes (log2)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = rep(c("#1F77B4", "#D62728"), each = 3)) +
  theme(
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line()
  )

ggsave("box_WGCNA_blue_meanlog2FC_160genes_allbios_20240504.png",
       plot = plot_object2, width = 15, height = 20, units = "cm", dpi = 300)

write.table(df3_all, "box_WGCNA_blue_meanlog2FC_160genes_allbios_20240504.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#================== ORIGINAL expression values from 13 bio files

colnames(Final_Table)[1] <- "gene_id"

Bio_names <- c("vertical_Bio1", "vertical_Bio2", "vertical_Bio3", "vertical_Bio4", "vertical_Bio5",
               "opposite_Bio1", "opposite_Bio2", "opposite_Bio3", "opposite_Bio4",
               "tension_Bio1", "tension_Bio2", "tension_Bio3", "tension_Bio4")

names <- c("Ptr", "PtrVertical02", "PtrVertical03", "PtrVertical04", "PtrVertical01",
           "PtrOW", "PtrOpposite02", "PtrOpposite03", "PtrOpposite04", 
           "PtrTW", "PtrTension02", "PtrTension03", "PtrTension04")

Input_df_s <- list()

for (i in seq_along(names)) {
  input_file <- paste0("geneUMI_TenX_", names[i], "_brown_blue_red.txt")
  df <- read.table(input_file, sep = "\t", header = TRUE)

  selected <- df[df$gene_id %in% Final_Table$gene_id, 1:6]
  colnames(selected)[2:6] <- paste0(c("brown_exp_", "blue_exp_", "red_exp_", 
                                      "blue_brown_log2fc_", "red_brown_log2fc_"), Bio_names[i])
  Input_df_s[[i]] <- selected
}

# Check gene_id consistency
gene_ID_list <- lapply(Input_df_s, function(df) unique(df$gene_id))
if (!all(sapply(gene_ID_list[-1], function(x) identical(x, gene_ID_list[[1]])))) {
  stop("Gene ID in all files are not identical!!!")
}

# Combine all tables by gene_id
Input_df_all <- do.call(cbind, Input_df_s)

# Remove redundant gene_id columns (keep first)
Input_df_all <- Input_df_all[, -c(seq(7, ncol(Input_df_all), 6))]

# Remove log2FC columns to get only expression
Input_df_expr_only <- Input_df_all[, !grepl("log2fc", colnames(Input_df_all))]

# Merge with Final_Table to get all-in-one expression + summary
Final_Table_Final <- merge(Input_df_expr_only, Final_Table, by = "gene_id")

write.table(Final_Table_Final,
            "box_WGCNA_blue_meanlog2FC_160genes_allbios_20240504_original.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
