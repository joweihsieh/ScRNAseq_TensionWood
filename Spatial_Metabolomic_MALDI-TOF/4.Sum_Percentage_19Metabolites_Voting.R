library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(reshape2)

# ------------------------ 1. Sum of percentage ------------------------

process_and_plot_by_group <- function(data, cutoff = 5, iqr_multiplier = 1, Treatments, Bios, Metabolites) {
  data %>%
    mutate(
      group = ceiling(x / 2),
      percentage = ifelse(percentage == 0, NA, percentage),
      percentage = ifelse(percentage < cutoff, 0, percentage)
    ) %>%
    group_by(group) %>%
    mutate(
      lower = quantile(percentage, 0.25, na.rm = TRUE) - iqr_multiplier * IQR(percentage, na.rm = TRUE),
      upper = quantile(percentage, 0.75, na.rm = TRUE) + iqr_multiplier * IQR(percentage, na.rm = TRUE),
      is_outlier = percentage < lower | percentage > upper
    ) %>%
    ungroup() %>%
    mutate(Treatments = Treatments, Bios = Bios, Metabolites = Metabolites) -> processed

  summary <- processed %>%
    filter(group %in% 1:6) %>%
    group_by(Name = paste(Treatments, Bios, Metabolites, sep = "_")) %>%
    summarise(sum_percentage = sum(percentage, na.rm = TRUE), count = sum(!is.na(percentage))) %>%
    mutate(
      Treatments = Treatments,
      Bios = Bios,
      Metabolites = Metabolites,
      Treatments_Bios = paste(Treatments, Bios, sep = "_")
    )
  
  list(processed = processed, summary = summary)
}

# auto-processing
process_directory <- function(input_dir, output_dir, bio_id = NULL) {
  files <- list.files(input_dir, pattern = "rgb_percent.txt$", full.names = TRUE)
  lapply(files, function(file_path) {
    file_name <- basename(file_path)
    parts <- strsplit(file_name, "_")[[1]]
    Treatments <- parts[1]
    Bios <- if (!is.null(bio_id)) bio_id else parts[2]
    Metabolites <- paste(parts[3], parts[4], sep = "_")

    data <- read.delim(file_path)
    if (bio_id == "Bio1") colnames(data) <- gsub("_grid", "", colnames(data))

    process_and_plot_by_group(data, Treatments = Treatments, Bios = Bios, Metabolites = Metabolites)
  })
}


results_Bio1 <- process_directory("20241205_stem+CHCA_image/output", "out", bio_id = "Bio1")
results_Bio2 <- process_directory("20241205_stem+CHCA_image_2/output", "out", bio_id = "Bio2")
results_Bio3 <- process_directory("20241205_stem+CHCA_image_3/output", "out", bio_id = "Bio3")

summary_all <- bind_rows(
  lapply(results_Bio1, `[[`, "summary"),
  lapply(results_Bio2, `[[`, "summary"),
  lapply(results_Bio3, `[[`, "summary")
)


norm_factor <- min(summary_all$count)
summary_all <- summary_all %>%
  mutate(sum_percentage_norm = (sum_percentage / count) * norm_factor)

write.table(summary_all, "19_metabolites_6groups_sum_percentage.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------ 2. Bar charts for sum of percentage ------------------------



plot_metabolites <- function(data, col = "sum_percentage", norm = FALSE, t_test = FALSE, output_dir = "plot") {
  if (!dir.exists(output_dir)) dir.create(output_dir)

  for (met in unique(data$Metabolites)) {
    df <- data %>% filter(Metabolites == met)
    mean_df <- df %>% group_by(Treatments) %>% summarise(mean_val = mean(.data[[col]], na.rm = TRUE), .groups = "drop")
    p <- ggplot(mean_df, aes(x = Treatments, y = mean_val)) +
      geom_bar(stat = "identity", fill = "gray", color = "black") +
      geom_point(data = df, aes(x = Treatments, y = .data[[col]]), size = 3, position = position_jitter(width = 0.2)) +
      labs(title = met, y = if (norm) "Normalized Sum %" else "Sum %") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(output_dir, paste0(met, "_barchart.png")), plot = p, width = 2, height = 6)
  }
}

plot_metabolites(summary_all, col = "sum_percentage", output_dir = "output_241212_group1-6_sum")
plot_metabolites(summary_all, col = "sum_percentage_norm", output_dir = "output_241212_group1-6_sum_norm")


# ------------------------ 3. Voting and Chi-square ------------------------

vote_df <- summary_all %>%
  group_by(Bios, Metabolites) %>%
  summarise(
    Tension = sum(sum_percentage_norm[Treatments == "Tension"], na.rm = TRUE),
    Normal = sum(sum_percentage_norm[Treatments == "Normal"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Comparison = ifelse(Tension > Normal, "Tension", "Normal"))

summary_vote <- vote_df %>%
  group_by(Metabolites) %>%
  summarise(
    Tension_win = sum(Comparison == "Tension"),
    Normal_win = sum(Comparison == "Normal"),
    .groups = "drop"
  )

vote_summary <- summary_vote %>%
  summarise(
    Tension_greater = sum(Tension_win > Normal_win),
    Normal_greater = sum(Normal_win > Tension_win)
  )

#chisq.test(as.numeric(vote_summary))
