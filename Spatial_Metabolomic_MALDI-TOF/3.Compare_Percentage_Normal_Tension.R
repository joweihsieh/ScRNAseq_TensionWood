
###############################
# Extract 12 grids and combine each two of them into one new grid

# Load required libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# ------------------------ 1. Unified File Processing Function ------------------------
process_single_file <- function(file_path, output_dir, cutoff = 5, iqr_multiplier = 1, remove_grid_suffix = FALSE) {
  file_name <- basename(file_path)
  file_parts <- strsplit(file_name, "_")[[1]]
  treatment <- file_parts[1]
  bios <- file_parts[2]
  metabolites <- paste0(file_parts[3], "_", file_parts[4])

  data <- read.delim(file_path, sep = "\t", header = TRUE)
  if (remove_grid_suffix) colnames(data) <- gsub("_grid", "", colnames(data))

  processed <- data %>%
    mutate(group = ceiling(x / 2),
           percentage = ifelse(percentage == 0, NA, percentage),
           percentage = ifelse(percentage < cutoff, 0, percentage)) %>%
    group_by(group) %>%
    mutate(
      lower = quantile(percentage, 0.25, na.rm = TRUE) - 1.5 * iqr_multiplier * IQR(percentage, na.rm = TRUE),
      upper = quantile(percentage, 0.75, na.rm = TRUE) + 1.5 * iqr_multiplier * IQR(percentage, na.rm = TRUE),
      is_outlier = ifelse(!is.na(percentage) & (percentage < lower | percentage > upper), TRUE, FALSE)
    ) %>% ungroup()

  summary_stats <- processed %>%
    group_by(group) %>%
    summarise(mean_percentage = mean(percentage, na.rm = TRUE), .groups = "drop") %>%
    filter(group <= 6)

  minmax_ratio <- min(summary_stats$mean_percentage, na.rm = TRUE) / max(summary_stats$mean_percentage, na.rm = TRUE)
  first_last_ratio <- summary_stats$mean_percentage[summary_stats$group == 1] /
                      summary_stats$mean_percentage[summary_stats$group == 6]

  summary_stats <- summary_stats %>%
    mutate(Treatments = treatment, Bios = bios, Metabolites = metabolites,
           Min_max_ratio = minmax_ratio, First_last_ratio = first_last_ratio)

  plot_path <- file.path(output_dir, paste0(treatment, "_", bios, "_", metabolites, "_mean.png"))
  p <- ggplot(summary_stats, aes(x = factor(group), y = mean_percentage)) +
    geom_bar(stat = "identity", fill = "black") +
    labs(x = "Radial Layer", y = "Mean Percentage",
         title = paste0(treatment, "_", bios, "_", metabolites,
                        "\nMin/Max = ", round(minmax_ratio, 2),
                        ", First/Last = ", round(first_last_ratio, 2))) +
    theme_minimal()
  ggsave(plot_path, p, width = 6, height = 8, dpi = 300)

  return(list(summary = summary_stats, processed = processed))
}

process_folder <- function(input_dir, output_dir, cutoff = 5, iqr_multiplier = 1, remove_grid_suffix = FALSE) {
  files <- list.files(input_dir, pattern = "rgb_percent.txt$", full.names = TRUE)
  results <- lapply(files, function(file_path) {
    process_single_file(file_path, output_dir, cutoff, iqr_multiplier, remove_grid_suffix)
  })
  summary_all <- bind_rows(lapply(results, `[[`, "summary"))
  processed_all <- bind_rows(lapply(results, `[[`, "processed"))
  return(list(summary = summary_all, processed = processed_all))
}

# ------------------------ 2. Voting Function ------------------------
voting_summary <- function(summary_df) {
  metabolite_list <- unique(summary_df$Metabolites)
  results <- data.frame(Metabolites = metabolite_list)

  for (bio in unique(summary_df$Bios)) {
    bio_df <- summary_df %>% filter(Bios == bio)

    # Get min ratio per metabolite
    min_ratios <- bio_df %>%
      group_by(Metabolites) %>%
      summarise(min_ratio = min(Min_max_ratio, na.rm = TRUE), .groups = "drop")

    qualified_mets <- min_ratios %>% filter(min_ratio < 0.5) %>% pull(Metabolites)

    tension_vs_normal <- bio_df %>%
      filter(Metabolites %in% qualified_mets) %>%
      filter(Treatments %in% c("Normal", "Tension")) %>%
      group_by(Metabolites, Treatments) %>%
      summarise(mean_ratio = mean(First_last_ratio, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Treatments, values_from = mean_ratio) %>%
      mutate(Tension_gt_Normal = Tension > Normal)

    results <- results %>%
      left_join(min_ratios %>% filter(min_ratio < 0.5), by = "Metabolites") %>%
      rename(!!paste0("min_ratio_", bio) := min_ratio) %>%
      left_join(tension_vs_normal %>% select(Metabolites, Tension_gt_Normal), by = "Metabolites") %>%
      rename(!!paste0("Tension_gt_Normal_", bio) := Tension_gt_Normal)
  }

  return(results)
}

# ------------------------ 3. Regression Voting ------------------------
regression_voting <- function(summary_df, voting_result) {
  metabolite_list <- unique(summary_df$Metabolites)
  bios <- unique(summary_df$Bios)
  results <- data.frame(Metabolites = metabolite_list)

  for (bio in bios) {
    qualified_mets <- voting_result %>% filter(!is.na(.data[[paste0("min_ratio_", bio)]])) %>% pull(Metabolites)

    norm_slopes <- c()
    tension_slopes <- c()
    decisions <- c()
    for (met in metabolite_list) {
      if (!(met %in% qualified_mets)) {
        norm_slopes <- c(norm_slopes, NA)
        tension_slopes <- c(tension_slopes, NA)
        decisions <- c(decisions, NA)
        next
      }
      df_sub <- summary_df %>% filter(Metabolites == met, Bios == bio)
      norm <- df_sub %>% filter(Treatments == "Normal")
      tension <- df_sub %>% filter(Treatments == "Tension")
      if (nrow(norm) > 1 && nrow(tension) > 1) {
        norm_slope <- coef(lm(mean_percentage ~ group, data = norm))[2]
        tension_slope <- coef(lm(mean_percentage ~ group, data = tension))[2]
        decision <- tension_slope < norm_slope
      } else {
        norm_slope <- tension_slope <- decision <- NA
      }
      norm_slopes <- c(norm_slopes, norm_slope)
      tension_slopes <- c(tension_slopes, tension_slope)
      decisions <- c(decisions, decision)
    }
    results[[paste0("Norm_Slope_", bio)]] <- norm_slopes
    results[[paste0("Tension_Slope_", bio)]] <- tension_slopes
    results[[paste0("Regression_", bio)]] <- decisions
  }

  return(results)
}

# ------------------------ 4. Plotting ------------------------
plot_metabolite_barplot <- function(df, output_dir) {
  treatment_colors <- c("Normal_Bio1" = "grey", "Tension_Bio1" = "black",
                        "Normal_Bio2" = "grey", "Tension_Bio2" = "black",
                        "Normal_Bio3" = "grey", "Tension_Bio3" = "black")

  for (metabolite in unique(df$Metabolites)) {
    subset_data <- df %>% filter(Metabolites == metabolite)
    subset_data$Treatments_Bios <- paste0(subset_data$Treatments, "_", subset_data$Bios)

    p <- ggplot(subset_data, aes(x = factor(group), y = mean_percentage, fill = Treatments_Bios)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = metabolite, x = "Layer (Outer to Inner)", y = "Mean %") +
      facet_wrap(~Treatments_Bios, scales = "fixed") +
      scale_fill_manual(values = treatment_colors) +
      theme_minimal()

    ggsave(file.path(output_dir, paste0("mean_percentage_", metabolite, ".png")), p, width = 10, height = 6, dpi = 300)
  }
}

# ------------------------ 5. Pipeline Runner ------------------------
run_full_pipeline <- function(input_dirs, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  all_summaries <- list()
  all_processed <- list()

  for (i in seq_along(input_dirs)) {
    bio_name <- paste0("Bio", i)
    remove_suffix <- bio_name == "Bio1"  # Assume Bio1 has _grid suffix
    res <- process_folder(input_dirs[i], output_dir, remove_grid_suffix = remove_suffix)
    all_summaries[[i]] <- res$summary
    all_processed[[i]] <- res$processed
  }

  combined_summary <- bind_rows(all_summaries)
  combined_processed <- bind_rows(all_processed)

  voting_result <- voting_summary(combined_summary)
  regression_result <- regression_voting(combined_summary, voting_result)

  write.table(voting_result, file.path(output_dir, "voting_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(regression_result, file.path(output_dir, "regression_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  plot_metabolite_barplot(combined_summary, output_dir)

  return(list(summary = combined_summary, processed = combined_processed,
              voting = voting_result, regression = regression_result))
}

input_dirs <- c(
  "20241205_stem+CHCA_image/output",
  "20241205_stem+CHCA_image_2/output",
  "20241205_stem+CHCA_image_3/output"
)

output_dir <- "output_250418"

results <- run_full_pipeline(input_dirs, output_dir)

plot_summary_sums(results$summary, output_dir)
