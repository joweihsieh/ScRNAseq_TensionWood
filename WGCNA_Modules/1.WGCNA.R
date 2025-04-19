

library(magrittr)
library(colorRamp2)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(tidyr)
library(GO.db)
library(class)
library(cluster)
library(impute)
library(Hmisc)
library(WGCNA)
library(dplyr)
library(hrbrthemes)

######################### Step1. Import data

Table <- readRDS("integration_PtrTensionWoodALL_13samples_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds")
Integ_Clu <- read.csv("integration_PtrTensionWoodALL_13samples_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base_PC30_assingedColors.csv")


PtrVert1 <- read.csv("geneUMI_TenX_PtrVertical01.csv")
PtrVert1$Barcode <- paste0("geneUMI_TenX_PtrVertical01_", PtrVert1$Barcode)
PtrTens2 <- read.csv("geneUMI_TenX_PtrTension02.csv")
PtrTens2$Barcode <- paste0("geneUMI_TenX_PtrTension02_", PtrTens2$Barcode)
PtrOp2 <- read.csv("geneUMI_TenX_PtrOpposite02.csv")
PtrOp2$Barcode <- paste0("geneUMI_TenX_PtrOpposite02_", PtrOp2$Barcode)

PtrVert1_Clu <- merge(PtrVert1, Integ_Clu, by = "Barcode")
PtrTens2_Clu <- merge(PtrTens2, Integ_Clu, by = "Barcode")
PtrOp2_Clu <- merge(PtrOp2, Integ_Clu, by = "Barcode")


PtrVert1_Clu_s <- PtrVert1_Clu[PtrVert1_Clu$Cluster %in% c(1,2,7),]
PtrTens2_Clu_s <- PtrTens2_Clu[PtrTens2_Clu$Cluster %in% c(1,2,7),]
PtrOp2_Clu_s <- PtrOp2_Clu[PtrOp2_Clu$Cluster %in% c(1,2,7),]


#####
#mycol <- colorpanel(n=40,low="white",high="darkred")

#data <- as.matrix(t(PtrVert1_Clu_s[2:34700]))

#out <- heatmap.2(data,key.title=NA,col=mycol,trace="none",
#    distfun = function(x)dist(x,method="euclidean"), hclustfun = function(x){hclust(x, method = 'complete')}, 
#    mar=c(8,20),dendrogram='row', Rowv=TRUE, cexCol = 0.8, cexRow = 0.65, Colv=FALSE,colsep=c(1:8))

#jpeg(filename="Heatmap_brown_blue_red.jpg", width=4500, height=3000, res=400)
#out
#dev.off()

######################### Step2. Fold change between brown to red or blue

df <- data.frame(matrix(0, 34699, 8))
colnames(df) <- c("gene_id", "brown_mean", "blue_mean", "red_mean", "blue_brown_fc", "red_brown_fc", "blue_brown_p", "red_brown_p")
df$gene_id <- colnames(PtrVert1_Clu_s) [2:34700]

for (i in 2:34700){
  print(i)
  sets <- PtrVert1_Clu_s[,c(1, i, 34704)]
  colnames(sets)[2] <- "expression"
  df[(i-1), "brown_mean"] <- mean(sets[sets$Cluster %in% 2, "expression"])
  df[(i-1), "blue_mean"] <- mean(sets[sets$Cluster %in% 1, "expression"])
  df[(i-1), "red_mean"] <- mean(sets[sets$Cluster %in% 7, "expression"])
  df[(i-1), "blue_brown_fc"] <- log((df[(i-1), "blue_mean"]+1)/(df[(i-1), "brown_mean"]+1),2)
  df[(i-1), "red_brown_fc"] <- log((df[(i-1), "red_mean"]+1)/(df[(i-1), "brown_mean"]+1),2)


  pairwise_result <- pairwise.t.test(sets$expression, sets$Cluster, p.adj = "none")
  pairwise_raw_p_values <- pairwise_result$p.value
  print("brown -> blue")
  df[(i-1), "blue_brown_p"] <- pairwise_raw_p_values[1,1]
  print("brown -> red")
  df[(i-1), "red_brown_p"] <- pairwise_raw_p_values[2,2]

}

write.table(df, "geneUMI_TenX_PtrVertical01_brown_blue_red.txt", quote = F, sep="\t", row.names = F, col.names = T)


#### tension
df2 <- data.frame(matrix(0, 34699, 8))
colnames(df2) <- c("gene_id", "brown_mean", "blue_mean", "red_mean", "blue_brown_fc", "red_brown_fc", "blue_brown_p", "red_brown_p")
df2$gene_id <- colnames(PtrTens2_Clu_s) [2:34700]

for (i in 2:34700){
  print(i)
  sets <- PtrTens2_Clu_s[,c(1, i, 34704)]
  colnames(sets)[2] <- "expression"
  df2[(i-1), "brown_mean"] <- mean(sets[sets$Cluster %in% 2, "expression"])
  df2[(i-1), "blue_mean"] <- mean(sets[sets$Cluster %in% 1, "expression"])
  df2[(i-1), "red_mean"] <- mean(sets[sets$Cluster %in% 7, "expression"])
  df2[(i-1), "blue_brown_fc"] <- log((df2[(i-1), "blue_mean"]+1)/(df2[(i-1), "brown_mean"]+1),2)
  df2[(i-1), "red_brown_fc"] <- log((df2[(i-1), "red_mean"]+1)/(df2[(i-1), "brown_mean"]+1),2)


  pairwise_result <- pairwise.t.test(sets$expression, sets$Cluster, p.adj = "none")
  pairwise_raw_p_values <- pairwise_result$p.value
  print("brown -> blue")
  df2[(i-1), "blue_brown_p"] <- pairwise_raw_p_values[1,1]
  print("brown -> red")
  df2[(i-1), "red_brown_p"] <- pairwise_raw_p_values[2,2]

}


write.table(df2, "geneUMI_TenX_PtrTension02_brown_blue_red.txt", quote = F, sep="\t", row.names = F, col.names = T)

#### op
df3 <- data.frame(matrix(0, 34699, 8))
colnames(df3) <- c("gene_id", "brown_mean", "blue_mean", "red_mean", "blue_brown_fc", "red_brown_fc", "blue_brown_p", "red_brown_p")
df3$gene_id <- colnames(PtrOp2_Clu_s) [2:34700]

for (i in 2:34700){
  print(i)
  sets <- PtrOp2_Clu_s[,c(1, i, 34704)]
  colnames(sets)[2] <- "expression"
  df3[(i-1), "brown_mean"] <- mean(sets[sets$Cluster %in% 2, "expression"])
  df3[(i-1), "blue_mean"] <- mean(sets[sets$Cluster %in% 1, "expression"])
  df3[(i-1), "red_mean"] <- mean(sets[sets$Cluster %in% 7, "expression"])
  df3[(i-1), "blue_brown_fc"] <- log((df3[(i-1), "blue_mean"]+1)/(df3[(i-1), "brown_mean"]+1), 2)
  df3[(i-1), "red_brown_fc"] <- log((df3[(i-1), "red_mean"]+1)/(df3[(i-1), "brown_mean"]+1), 2)


  pairwise_result <- pairwise.t.test(sets$expression, sets$Cluster, p.adj = "none")
  pairwise_raw_p_values <- pairwise_result$p.value
  print("brown -> blue")
  df3[(i-1), "blue_brown_p"] <- pairwise_raw_p_values[1, 1]
  print("brown -> red")
  df3[(i-1), "red_brown_p"] <- pairwise_raw_p_values[2, 2]

}


write.table(df3, "geneUMI_TenX_PtrOpposite02_brown_blue_red.txt", quote = F, sep="\t", row.names = F, col.names = T)


######################### Step3. WGCNA

### pre
df <- read.table("geneUMI_TenX_PtrVertical01_brown_blue_red.txt",header=T)
df2 <- read.table("geneUMI_TenX_PtrTension02_brown_blue_red.txt",header=T)
df3 <- read.table("geneUMI_TenX_PtrOpposite02_brown_blue_red.txt",header=T)


df_combine <- cbind(df[,c("blue_brown_fc", "red_brown_fc")], df2[,c("blue_brown_fc", "red_brown_fc")], df3[,c("blue_brown_fc", "red_brown_fc")])
rownames(df_combine) <- df$gene_id
df_combine$delta <- apply(df_combine, 1 , function(x) max(x) - min(x))
df_combine2 <- df_combine[(df_combine$delta) >= 0.1,]
data <- as.matrix(df_combine2[, 1:6])



datExpr <- t(data)
A <- adjacency(datExpr)
TOM_sim <- TOMsimilarity(A, TOMType = "signed")
TOM_dist <- TOMdist(A)
clu <- hclust(as.dist(TOM_dist), method="average")   # get hierarchical clustering on TOM distance with average linking
color_clu <- cutreeStaticColor(clu, cutHeight = 0.75, minSize = 50) # cut static and get colors for modules


## plot module clustring
png(filename = "WGCNA_clustering_mean_module_log_expression.png", width = 1600, height = 1600, res = 200);

par(mfrow = c(2,1),mar = c(2,4,1,1))
plot(clu, main = "", labels = F, xlab = "", sub = "")
plotColorUnderTree(clu,colors=data.frame(module = color_clu), rowText = data.frame(module = color_clu))

dev.off()

## get data
table <- as.data.frame(cbind(data, color_clu))
table$row_names <- rownames(table)

#module: gene_id, value and color for GO
colnames(df_combine2)[ncol(df_combine2)] <- c("row_names")
df_combine2$gene_id <- rownames(df_combine2)
ID_number <- df_combine2[, c("gene_id", "row_names")]
table2 <- merge(table, ID_number, by = "row_names")
write.table(table2, "WGCNA_modules_mean_log_expression_value_colors_gene_id.txt", row.names = F, quote = F, col.names = T, sep="\t")




######################### Step 4. Plotting heatmap for each module
 


table2 <- read.table("WGCNA_modules_mean_log_expression_value_colors_gene_id.txt",sep="\t",header=T)
colnames(table2) <- c("row_names", "blue_brown_fc_vertical", "red_brown_fc_vertical", 
  "blue_brown_fc_tension", "red_brown_fc_tension", "blue_brown_fc_opposite", "red_brown_fc_opposite", "color_clu", "gene_id")


table2 <- table2[,c("row_names", "blue_brown_fc_vertical", "blue_brown_fc_opposite", 
  "blue_brown_fc_tension", "red_brown_fc_vertical", "red_brown_fc_opposite", "red_brown_fc_tension", "color_clu", "gene_id")]
colors2 <- unique(table2$color_clu)




## color setting

start_color <- "#AC9300"  
end_color <- "#21FF58"    

color_interpolator <- colorRampPalette(c(start_color, "white", end_color))
hmcol2 <- color_interpolator(100)
print(interpolated_colors)

## heatmap

for (i in colors2){
  module = i
  data2 = as.matrix(table2[table2$color_clu == module, c(2:7)])
  names = paste0("heatmap_WGCNA_mean_module_log_expression_",module, "20240507.png")
  png(filename=names, width = 2000, height = 1600, res=200);
  heatmap.2(data2, hclustfun=function(d) hclust(d, method="average"), 
    main=module, distfun = function(x) dist(x, method = "canberra"), 
    dendrogram = "row", labRow = NA, col = hmcol2, trace = "none", 
    sepwidth = c(0.02,0.02,0.02), cex.lab = 0.5, mar = c(15,2), Colv = F, las = 2)
  dev.off()
}

