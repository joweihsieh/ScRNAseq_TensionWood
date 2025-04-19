setwd("/Users/joweihsieh/Dropbox/YCL/tension_wood/phosphoproteome/20241221/")

library(readxl)
library(ggplot2)



sample_list_path <- "NW-TW-OW-phos_DIAhybridDDA_R1136_20241218_Peptide Quant Report_20241221_reformat.xlsx"

Genefuncs <- read.csv("/Users/joweihsieh/Dropbox/YCL/tension_wood/phosphoproteome/20240827_Phosphoproteome xylem development references/Ptrichocarpa_533_v4.1.annotation_info.txt", sep = "\t")
Genefuncs[Genefuncs == ""] <- "NA"


groups <- c("NW1", "NW2", "NW3", "NW4", 
  "D1_TW1", "D1_TW2", "D1_TW3", "D1_TW4", "D1_TW5",
  "D7_TW1", "D7_TW2", "D7_TW3", "D7_TW4",
  "D1_OW1", "D1_OW2", "D1_OW3", "D1_OW4", "D1_OW5",
  "D7_OW1", "D7_OW2", "D7_OW3", "D7_OW4")


target_groups <- c("NW1", "D1_TW1", "D7_TW1", "D1_OW1", "D7_OW1")


############# separate mixed information in PG.ProteinAccessions
df <- read_excel(sample_list_path)
df$PG.ProteinAccessions <- as.character(df$PG.ProteinAccessions)

process_accession <- function(accession) {
  if (is.na(accession) || !nzchar(accession)) {
    return(data.frame(
      PG.ProteinAccessions = NA,
      PG.Genes = NA,
      PG.ProteinDescriptions = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  parts <- strsplit(accession, ">")[[1]]
  parts <- parts[nzchar(parts)]  
  
  result <- lapply(parts, function(part) {
    PG.ProteinAccessions <- sub(" (pacid=.*)", "", part)  
    PG.Genes <- sub(".*(pacid=[^ ]+).*", "\\1", part)    
    PG.ProteinDescriptions <- sub(".*pacid=[^ ]+ ", "", part) 
    
    data.frame(protein = PG.ProteinAccessions, 
               pacid = PG.Genes, 
               PG.ProteinDescriptions = PG.ProteinDescriptions, 
               stringsAsFactors = FALSE)
  })
  
  result_df <- do.call(rbind, result)
  combined <- data.frame(
    PG.ProteinAccessions = gsub(";+", ";", paste(result_df$protein, collapse = ";")),
    PG.Genes = gsub(";+", ";", paste(result_df$pacid, collapse = ";")),
    PG.ProteinDescriptions = gsub(";+", ";", paste(result_df$PG.ProteinDescriptions, collapse = ";")),
    stringsAsFactors = FALSE
  )
  
  return(combined)
}

processed_results <- do.call(rbind, lapply(df$PG.ProteinAccessions, process_accession))

df$PG.ProteinAccessions <- processed_results$PG.ProteinAccessions
df$PG.Genes <- processed_results$PG.Genes
df$PG.ProteinDescriptions <- processed_results$PG.ProteinDescriptions



############# convert filtered into zero

matching_cols <- grepl(paste(groups, collapse = "|"), colnames(df))

# Replace "Filtered" with 0, convert to numeric, and replace NA with 0
df[, matching_cols] <- apply(df[, matching_cols], 2, function(x) {
  # Replace "Filtered" with 0
  x <- ifelse(x == "Filtered", 0, x)
  # Convert to numeric
  x <- as.numeric(x)
  # Replace NA with 0
  x[is.na(x)] <- 0
  return(x)
})

#############

process_phospho_data <- function(df_original, matching_cols, target_groups, control_indices, treatment_indices, control_prefix, treatment_prefix, label_prefix, Genefuncs) {
  df <- df_original
  df_original[, target_groups] <- df_original[, target_groups] + 0.00001

  t_test_results <- apply(df_original[, matching_cols], 1, function(row) {
    t.test(as.numeric(row[control_indices]), as.numeric(row[treatment_indices]))  # Control vs Treatment
  })

  p_values <- sapply(t_test_results, function(res) res$p.value)

  df_original[, target_groups] <- df_original[, target_groups] - 0.00001


  control_mean <- apply(df_original[, matching_cols], 1, function(row) {
    control_mean <- mean(as.numeric(row[control_indices]), na.rm = TRUE)
    return(control_mean)
  })

  treatment_mean <- apply(df_original[, matching_cols], 1, function(row) {
    treatment_mean <- mean(as.numeric(row[treatment_indices]), na.rm = TRUE)
    return(treatment_mean)
  })

  log2fold_change <- apply(df_original[, matching_cols], 1, function(row) {
    control_mean <- mean(as.numeric(row[control_indices]), na.rm = TRUE)
    treatment_mean <- mean(as.numeric(row[treatment_indices]), na.rm = TRUE)
    
 
    return(log((treatment_mean + 0.0001) / (control_mean + 0.0001), 2))
  })


  df_original[[paste0("mean_", control_prefix)]] <- control_mean
  df_original[[paste0("mean_", treatment_prefix)]] <- treatment_mean


  log2fold_change_col <- paste0("log2FC_", label_prefix)
  sig_col <- paste0("sig_", label_prefix)
  p_value_col <- paste0("pvalue_", label_prefix)

  df_original[[log2fold_change_col]] <- log2fold_change
  df_original[[p_value_col]] <- p_values

  df_original[[sig_col]] <- with(df_original, ifelse(p_values < 0.05 & abs(df_original[[log2fold_change_col]]) >= 1, "Significant", "Not Significant"))


  df_original$peptideName <- sapply(df_original$PG.ProteinAccessions, function(x) {
    parts <- strsplit(x, ";")[[1]]
    return(parts[1])
  })

  df_original_functions <- merge(df_original, Genefuncs, by = "peptideName", all.x = TRUE)
  return(df_original_functions)
}




TV_D7 <- process_phospho_data(
  df_original = df,
  matching_cols = matching_cols,
  target_groups = target_groups,
  control_indices = 1:4,      # normal
  treatment_indices = 10:13,    # tension
  control_prefix = "Normal",
  treatment_prefix = "Tension_D7",
  label_prefix = "Tension_D7_vs_Normal",
  Genefuncs = Genefuncs

)



OV_D7 <- process_phospho_data(
  df_original = df,
  matching_cols = matching_cols,
  target_groups = target_groups,
  control_indices = 1:4,      # normal
  treatment_indices = 19:22,    # opposite
  control_prefix = "Normal",
  treatment_prefix = "Opposite_D7",
  label_prefix = "Opposite_D7_vs_Normal",
  Genefuncs = Genefuncs

)


write.table(TV_D7, file = "Phospho_Peptide_Quant_ttest_TV_D7.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(OV_D7, file = "Phospho_Peptide_Quant_ttest_OV_D7.txt", sep = "\t", row.names = FALSE, quote = FALSE)


TOC132 <- c("Potri.009G131200")
MPKKK3 <- c("Potri.005G062500", "Potri.007G106800")
 

TV_D7[TV_D7$locusName%in%TOC132,]
TV_D7[TV_D7$locusName%in%MPKKK3,]


OV_D7[OV_D7$locusName%in%TOC132,]
OV_D7[OV_D7$locusName%in%MPKKK3,]












