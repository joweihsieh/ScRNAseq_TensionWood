setwd("/Users/joweihsieh/Dropbox/YCL/tension_wood/phosphoproteome/20240905")
library(ggplot2)


########################################## Get significance phosphopeptides

df <- read.csv("Phospho_Peptide Quant Report.csv")

groups <- c("VerticalBio1.P", "VerticalBio2.P", 
	"TensionBio1.P", "TensionBio2.P", 
	"OppositeBio1.P", "OppositeBio2.P")

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

# add a small value to one of the replicate of all groups
df[, matching_cols][, c(1, 3, 5)] <- df[, matching_cols][, c(1, 3, 5)] + 0.00001


# Run t-tests for each comparison
# 1. Vertical vs. Tension
t_test_TV <- apply(df[, matching_cols], 1, function(row) {
  t.test(as.numeric(row[1:2]), as.numeric(row[3:4]))  # Vertical vs Tension
})


# 2. Vertical vs. Opposite
t_test_OV <- apply(df[, matching_cols], 1, function(row) {
  t.test(as.numeric(row[1:2]), as.numeric(row[5:6]))  # Vertical vs Opposite
})


# 3. Tension vs. Opposite
t_test_TO <- apply(df[, matching_cols], 1, function(row) {
  t.test(as.numeric(row[3:4]), as.numeric(row[5:6]))  # Tension vs Opposite
})

# Extract p-values from t-test results
p_values_TV <- sapply(t_test_TV, function(res) res$p.value)
p_values_OV <- sapply(t_test_OV, function(res) res$p.value)
p_values_TO <- sapply(t_test_TO, function(res) res$p.value)

# Add the p-values back to the dataframe
df$p_value_TV <- p_values_TV 
df$p_value_OV <- p_values_OV 
df$p_value_TO <- p_values_TO 


df[, matching_cols][, c(1, 3, 5)] <- df[, matching_cols][, c(1, 3, 5)] - 0.00001


# Calculate fold change (mean ratio) for Tension  vs. Vertical
log2fold_change_TV <- apply(df[, matching_cols], 1, function(row) {
  vertical_mean <- mean(as.numeric(row[1:2]), na.rm = TRUE)
  tension_mean <- mean(as.numeric(row[3:4]), na.rm = TRUE)

  return(log((tension_mean + 0.0001) / (vertical_mean + 0.0001), 2))

})

# Calculate fold change (mean ratio) for Opposite vs. Vertical
log2fold_change_OV <- apply(df[, matching_cols], 1, function(row) {
  vertical_mean <- mean(as.numeric(row[1:2]), na.rm = TRUE)
  opposite_mean <- mean(as.numeric(row[5:6]), na.rm = TRUE)

  return(log((opposite_mean + 0.0001) / (vertical_mean + 0.0001), 2))

})

# Calculate fold change (mean ratio) for Tension vs. Opposite
log2fold_change_TO <- apply(df[, matching_cols], 1, function(row) {
  tension_mean <- mean(as.numeric(row[3:4]), na.rm = TRUE)
  opposite_mean <- mean(as.numeric(row[5:6]), na.rm = TRUE)

  return(log((tension_mean + 0.0001) / (opposite_mean + 0.0001), 2))

})

# Add fold changes back to the dataframe
df$log2fold_change_TV <- log2fold_change_TV
df$log2fold_change_OV <- log2fold_change_OV
df$log2fold_change_TO <- log2fold_change_TO



df$peptideName <- sapply(df$PG.ProteinAccessions, function(x) {
  parts <- strsplit(x, ";")[[1]]
  return(parts[1])
})


###### gene info

Genefuncs <- read.csv("/Users/joweihsieh/Dropbox/YCL/tension_wood/phosphoproteome/20240827_Phosphoproteome xylem development references/Ptrichocarpa_533_v4.1.annotation_info.txt", sep = "\t")
Genefuncs[Genefuncs == ""] <- "NA"



df_functions <- merge(df, Genefuncs, by = "peptideName", all.x = T)

df_functions$sig_TV <- with(df_functions, ifelse(p_value_TV < 0.05 & abs(log2fold_change_TV) >= 1, "Significant", "Not Significant"))
df_functions$sig_OV <- with(df_functions, ifelse(p_value_OV < 0.05 & abs(log2fold_change_OV) >= 1, "Significant", "Not Significant"))
df_functions$sig_TO <- with(df_functions, ifelse(p_value_TO < 0.05 & abs(log2fold_change_TO) >= 1, "Significant", "Not Significant"))


write.table(df_functions, file = "Phospho_Peptide_Quant_ttest_functions.txt", sep = "\t", row.names = FALSE, quote = FALSE)

