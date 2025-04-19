################# We directly use significance IDs based on their ABA 1hr stastistic results
library(dplyr)
library(tidyr)
library(readxl)


setwd("/Users/joweihsieh/Dropbox/YCL/tension_wood/phosphoproteome/20241223_ABA_JAIle/")

#################################################### Our data
############# if PG.ProteinAccessions is mixed, we should separate mixed information in PG.ProteinAccessions

combine_df_sig_clean_p <- read.csv("/Users/joweihsieh/Dropbox/YCL/tension_wood/phosphoproteome/20241221/Phospho_Peptide_Quant_MWU_Test_TV_D7_twobatches_sig_only.txt", sep = "\t", header = T)
columns_to_check <- colnames(combine_df_sig_clean_p)[grepl("Relative_FC", colnames(combine_df_sig_clean_p))]

####################################################
#################################################### Step2. matching

Orth_Table <- read.table("/Users/joweihsieh/Dropbox/YCL/tension_wood/phosphoproteome/20241223_ABA_JAIle/Orthogroups_reformat.tsv", sep = "\t", header = T)


#################################################### ABA Ref 1

Orth_Table <- read.table("/Users/joweihsieh/Dropbox/YCL/tension_wood/phosphoproteome/20241223_ABA_JAIle/Orthogroups_reformat.tsv", sep = "\t", header = T)
ABA_At_dedup_ID <- unique(ABA_At_dedup$best_arabi_gene)
ABA_At_orth <- unique(Orth_Table[Orth_Table$Genes %in% ABA_At_dedup_ID, "Orthogroup"])

Orth_Table_ABA_At_orth_Ptr <- Orth_Table[Orth_Table$Orthogroup %in% ABA_At_orth & Orth_Table$Species %in% "Ptrichocarpa",]
Orth_Table_ABA_At_orth_Ptr_ID <- unique(Orth_Table_ABA_At_orth_Ptr$Genes)


write.table(Orth_Table_ABA_At_orth_Ptr, "20241224_ABA_JA_phosphopeptides_AthDataOnly_Ptr_OG.txt", sep = "\t", quote = F, col.names = T, row.names = F)


#################################################### ABA Ref 2
ABA2 <- read_excel("20241229_ABA responsive phosphopeptides_Dr. Hsu.xlsx")
ABA2_At <- ABA2[grepl("AT", ABA2$`Gene ID`) ,]
ABA2_At$best_arabi_gene <- sub("\\.1$", "", ABA2_At$`Gene ID`)
ABA2_At_dedup_ID <- unique(ABA2_At$best_arabi_gene)
ABA2_At_orth <- unique(Orth_Table[Orth_Table$Genes %in% ABA2_At_dedup_ID, "Orthogroup"])


Orth_Table_ABA2_At_orth_Ptr <- Orth_Table[Orth_Table$Orthogroup %in% ABA2_At_orth & Orth_Table$Species %in% "Ptrichocarpa",]
Orth_Table_ABA2_At_orth_Ptr_ID <- unique(Orth_Table_ABA2_At_orth_Ptr$Genes)

write.table(Orth_Table_ABA2_At_orth_Ptr, "20241229_ABA responsive phosphopeptides_Dr. Hsu_Ptr_OG.txt", sep = "\t", quote = F, col.names = T, row.names = F)


#################################################### JAIle

JA <- read_excel("20250111_JA_phosphopeptides.xlsx")
JA_At <- JA[grepl("AT", JA$`Accession nember`) ,]
JA_At <- JA_At %>%
  mutate(`best_arabi_gene` = gsub("\\.\\d+$", "", `Accession nember`))
JA_At_dedup_ID <- unique(JA_At$best_arabi_gene)
JA_At_orth <- unique(Orth_Table[Orth_Table$Genes %in% JA_At_dedup_ID, "Orthogroup"])


Orth_Table_JA_At_orth_Ptr <- Orth_Table[Orth_Table$Orthogroup %in% JA_At_orth & Orth_Table$Species %in% "Ptrichocarpa",]
Orth_Table_JA_At_orth_Ptr_ID <- unique(Orth_Table_JA_At_orth_Ptr$Genes)


write.table(Orth_Table_JA_At_orth_Ptr, "20250111_JA_phosphopeptides_OG.txt", sep = "\t", quote = F, col.names = T, row.names = F)
