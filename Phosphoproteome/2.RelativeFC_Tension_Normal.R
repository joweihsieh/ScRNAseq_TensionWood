###### Phosphoproteomics data were analyzed using relative fold change. After performing statistical testing, a table of significantly phosphorylated peptides was retained.
###### Tension D7 vd Normal
library(dplyr)

setwd("/Users/joweihsieh/Dropbox/YCL/tension_wood/phosphoproteome/20241221")

TV_D7_df1 <- read.csv("/Users/joweihsieh/Dropbox/YCL/tension_wood/phosphoproteome/20240905/Phospho_Peptide_Quant_ttest_functions.txt", sep = "\t", header = T)
TV_D7_df2 <- read.csv("Phospho_Peptide_Quant_ttest_TV_D7.txt", sep = "\t", header = T)

TV_D7_df1_s <- TV_D7_df1[, c(1:11)]
TV_D7_df2_s <- TV_D7_df2[, c(1:11, 17:20, 35:45)]



get_matching_accession <- function(precursor_id1, accession1, precursor_id2, accession2) {
  if (is.na(precursor_id1) || is.na(precursor_id2)) {
    return(NA)
  }
  
  if (precursor_id1 != precursor_id2) {
    return(NA)
  }
  
  accessions1 <- unlist(strsplit(accession1, ";"))
  accessions2 <- unlist(strsplit(accession2, ";"))
  
  matching <- intersect(accessions1, accessions2)
  
  if (length(matching) > 0) {
    return(matching[1])
  } else {
    return(NA)
  }
}

TV_D7_df1_s$matching_accession <- mapply(
  get_matching_accession,
  TV_D7_df1_s$EG.PrecursorId,
  TV_D7_df1_s$PG.ProteinAccessions,
  TV_D7_df2_s$EG.PrecursorId[match(TV_D7_df1_s$EG.PrecursorId, TV_D7_df2_s$EG.PrecursorId)],
  TV_D7_df2_s$PG.ProteinAccessions[match(TV_D7_df1_s$EG.PrecursorId, TV_D7_df2_s$EG.PrecursorId)]
)

TV_D7_df2_s$matching_accession <- mapply(
  get_matching_accession,
  TV_D7_df2_s$EG.PrecursorId,
  TV_D7_df2_s$PG.ProteinAccessions,
  TV_D7_df1_s$EG.PrecursorId[match(TV_D7_df2_s$EG.PrecursorId, TV_D7_df1_s$EG.PrecursorId)],
  TV_D7_df1_s$PG.ProteinAccessions[match(TV_D7_df2_s$EG.PrecursorId, TV_D7_df1_s$EG.PrecursorId)]
)

TV_D7_df1_s_matched <- TV_D7_df1_s[!is.na(TV_D7_df1_s$matching_accession), ]
TV_D7_df2_s_matched <- TV_D7_df2_s[!is.na(TV_D7_df2_s$matching_accession), ]



TV_D7_df1_s_matched$full_name <- paste0(TV_D7_df1_s_matched$matching_accession, "_", TV_D7_df1_s_matched$EG.PrecursorId)
TV_D7_df2_s_matched$full_name <- paste0(TV_D7_df2_s_matched$matching_accession, "_", TV_D7_df2_s_matched$EG.PrecursorId)


combine_df <- merge(TV_D7_df1_s_matched, TV_D7_df2_s_matched[,8:ncol(TV_D7_df2_s_matched)], by = "full_name")


combine_df <- combine_df %>%
  mutate(across(matches("^(VerticalBio|TensionBio|NW|D7_TW)"), ~ .x + 0.00001)) %>%
  mutate(VerticalBio_Avg = rowMeans(select(., starts_with("VerticalBio")))) %>%
  mutate(across(starts_with("TensionBio"), ~ .x / VerticalBio_Avg, .names = "Relative_FC_{.col}")) %>%
  mutate(NW_Avg = rowMeans(select(., starts_with("NW")))) %>%
  mutate(across(starts_with("D7_TW"), ~ .x / NW_Avg, .names = "Relative_FC_{.col}")) %>%
  mutate(VerticalBio_Avg = 1, NW_Avg = 1)



##### MWU_test_pvalue

avg_columns <- grep("Avg", names(combine_df), value = TRUE)
fc_columns <- grep("Relative_FC_", names(combine_df), value = TRUE)


for (i in 1:nrow(combine_df)) {
    fc_data <- as.numeric(combine_df[i, fc_columns])
    avg_data <- as.numeric(combine_df[i, avg_columns])
    test_result <- wilcox.test(fc_data, avg_data, paired = FALSE, exact = TRUE)  
    combine_df[i, "MWU_test_pvalue"] <- test_result$p.value / 2
}


combine_df$Sig <- ifelse(combine_df$MWU_test_pvalue < 0.05, "Significant", "Not_Significant")
combine_df_sig <- combine_df[combine_df$Sig == c("Significant"),]

combine_df_sig_clean <- combine_df_sig[!apply(combine_df_sig, 1, function(row) all(is.na(row))), ]



combine_df_sig_clean_p <- combine_df_sig_clean[grepl("Phospho", combine_df_sig_clean$EG.PrecursorId),]



columns_to_check <- colnames(combine_df_sig_clean_p)[grepl("Relative_FC", colnames(combine_df_sig_clean_p))]
combine_df_sig_clean_p$Trend <- apply(combine_df_sig_clean_p[, columns_to_check], 1, function(row) {
  if (all(row > 1, na.rm = TRUE)) {
    return("Up")
  } 
  else if (all(row < 1, na.rm = TRUE)) {
    return("Down")
  } 
  else {
    return("Mix")
  }
})


write.table(combine_df_sig_clean_p, "Phospho_Peptide_Quant_MWU_Test_TV_D7_twobatches_sig_only.txt", sep = "\t", quote = F, row.names = F, col.names = T)

combine_df_sig_clean_p[combine_df_sig_clean_p$matching_accession.x %in% "Potri.015G142600.1.p",]
combine_df_sig_clean_p[combine_df_sig_clean_p$matching_accession.x %in% "Potri.005G062500.1.p",]
combine_df_sig_clean_p[combine_df_sig_clean_p$matching_accession.x %in% "Potri.007G106800.1.p",]


combine_df_sig_clean_p$PG.GeneAccessions <- gsub("\\.\\d+\\.p", ".v4.1", combine_df_sig_clean_p$PG.ProteinAccessions)


combine_df_sig_expanded <- combine_df_sig_clean_p %>%
  mutate(PG.GeneAccessions_split = strsplit(PG.GeneAccessions, ";")) %>%
  unnest(PG.GeneAccessions_split)


