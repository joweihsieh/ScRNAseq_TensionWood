
#######https://data.jgi.doe.gov/refine-download/phytozome?q=populus&expanded=Phytozome-533
df <- read.csv("/Users/joweihsieh/Dropbox/YCL/tension_wood/Module1_GO/Ptrichocarpa_533_v4.1.synonym.txt", sep = "\t", header = F)

df$POPTR_Combined <- apply(df[, -1], 1, function(row) {
  paste(row[grep("^POPTR", row)], collapse = ",")
})


result <- df[, c("V1", "POPTR_Combined")]
result$V1_base <- sub("\\.[0-9]+$", "", result$V1)
colnames(result) <- c("protein_v4_1", "v3", "v4_1")

write.table(result, "/Users/joweihsieh/Dropbox/YCL/tension_wood/Module1_GO/Ptrichocarpa_533_v4.1.synonym_20241227.txt", sep = "\t", quote = T, row.names = F, col.names = T)