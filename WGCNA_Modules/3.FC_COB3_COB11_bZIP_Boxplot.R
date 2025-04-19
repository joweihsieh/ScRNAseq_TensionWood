setwd("/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/all_UMI_tables")

plot_gene_expression_boxplot <- function(
  gene_id,
  file_paths,
  file_labels,
  fc_columns = c("blue_brown_fc", "red_brown_fc"),
  output_prefix = "boxplot_output"
) {
  library(magrittr)
  library(dplyr)
  library(reshape2)
  library(ggplot2)

  umi_tables <- lapply(file_paths, function(x) read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE))

  umi_tables[[length(umi_tables)]] <- umi_tables[[length(umi_tables)]][, c("gene_id", fc_columns)]
  rownames(umi_tables[[length(umi_tables)]]) <- umi_tables[[length(umi_tables)]]$gene_id
  colnames(umi_tables[[length(umi_tables)]]) <- c("gene_id", paste0(fc_columns, "_vertical"))


  filtered_tables <- lapply(seq_along(umi_tables), function(i) {
    df <- umi_tables[[i]]
    df <- df[rownames(df) %in% gene_id , , drop = FALSE]
    if (!"gene_id" %in% colnames(df)) 
    df$gene_id <- rownames(df)
    colnames(df) <- paste0(colnames(df), "_bio", i)
    df$condition <- rownames(df)
    return(df)
  })

  merged_df <- Reduce(function(x, y) merge(x, y, by = "condition"), filtered_tables)
  rownames(merged_df) <- merged_df$condition

  condition_names <- unique(gsub("_bio\\d+$", "", grep("bio", colnames(merged_df), value = TRUE)))
  condition_names <- condition_names[!condition_names %in% "gene_id"]
  condition_data <- lapply(condition_names, function(cond) {
    subdf <- merged_df[, grepl(cond, colnames(merged_df)), drop = FALSE]
    subdf[] <- lapply(subdf, function(x) if(is.character(x)) as.numeric(x) else x)

    subdf[[paste0(cond, "_mean")]] <- rowMeans(subdf)
    return(subdf)
  })

  mean_df <- do.call(cbind, condition_data)
  mean_df$Gene.ID <- rownames(mean_df)
  long_df <- melt(mean_df, id.vars = "Gene.ID")

  long_df_nomean <- long_df[grepl("bio", long_df$variable), ]
  long_df_nomean$treatment <- gsub("_bio\\d+$", "", long_df_nomean$variable)
  long_df_nomean$treatment <- factor(long_df_nomean$treatment, levels = c(
    "blue_brown_fc_vertical", "blue_brown_fc_opposite", "blue_brown_fc_tension",
    "red_brown_fc_vertical", "red_brown_fc_opposite", "red_brown_fc_tension"
  ))

  plot_object <- ggplot(long_df_nomean, aes(x = treatment, y = value, fill = treatment)) +
    geom_boxplot(width = 0.5) +
    geom_point(color = "black", size = 3) +
    labs(
      title = paste0("Gene: ", gene_id),
      x = "Conditions",
      y = "Expression changes (log2FC)"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c(rep("#1F77B4", 3), rep("#D62728", 3))) +
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

  print(plot_object)

  ggsave(
    filename = paste0(output_prefix, "_", gene_id, ".png"),
    plot = plot_object,
    width = 15, height = 20, units = "cm", dpi = 300
  )

  # pairwise t-test
  test_result <- pairwise.t.test(long_df_nomean$value, long_df_nomean$treatment, p.adj = "none", alternative = "greater")

  write.table(long_df_nomean,
              paste0(output_prefix, "_", gene_id, ".txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  return(test_result)
}



file_paths <- c(
  "geneUMI_TenX_PtrVertical01_PtrTension02_PtrOpposite02_brown_blue_red.txt",
  "geneUMI_TenX_PtrVertical02_PtrTW_PtrOW_brown_blue_red.txt",
  "geneUMI_TenX_PtrVertical03_PtrTension03_PtrOpposite03_brown_blue_red.txt",
  "geneUMI_TenX_PtrVertical04_PtrTension04_PtrOpposite04_brown_blue_red.txt",
  "geneUMI_TenX_Ptr_brown_blue_red.txt"
)

bZIP <- plot_gene_expression_boxplot(
  gene_id = "Potri.007G006900.v4.1",
  file_paths = file_paths,
  file_labels = paste0("bio", 1:5),
  output_prefix = "bZIP_expression"
)

COB3 <- plot_gene_expression_boxplot(
  gene_id = "Potri.004G117200.v4.1",
  file_paths = file_paths,
  file_labels = paste0("bio", 1:5),
  output_prefix = "COB3_expression_0418"
)

COB11 <- plot_gene_expression_boxplot(
  gene_id = "Potri.015G060100.v4.1",
  file_paths = file_paths,
  file_labels = paste0("bio", 1:5),
  output_prefix = "COB11_expression"
)

