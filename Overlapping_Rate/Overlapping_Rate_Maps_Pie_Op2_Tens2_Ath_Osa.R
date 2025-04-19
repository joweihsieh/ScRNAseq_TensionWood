library(Seurat)
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(igraph)
library(Matrix)
library(MASS)
library(ggplot2)
oriPar = par(no.readonly=T)

setwd("/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Multi_species_analysis/all_data_rds/")


############################################################ Step 0. Convert geneid into clusterID due to interspecies

### Opposite Bio2
PtrOp2 <- read.table("PtrOp2_filitered_feature_bc_matrix/features.tsv.gz",sep="\t")
colnames(PtrOp2) <- c("geneName","geneName_2","gene_expression")

Ortho_table <- readRDS("RDS_ortho_cluster_table.rds")
Ortho_table_Ptr <- Ortho_table[Ortho_table$speciesName=="PoT",]

PtrOp2_ortho <- merge(PtrOp2, Ortho_table_Ptr, by = "geneName")
stopifnot(nrow(PtrOp2_ortho) == nrow(PtrOp2))


PtrOp2_ortho$clusterID_2 <- paste0("Cluster_",PtrOp2_ortho$clusterID)
PtrOp2_ortho$clusterID <- PtrOp2_ortho$clusterID_2
write.table(PtrOp2_ortho[,c("clusterID","clusterID_2","gene_expression")], file = gzfile("features.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


PtrOp2_b <- read.table("PtrOp2_filitered_feature_bc_matrix/barcodes.tsv.gz",sep="\t")
PtrOp2_b$V1 <- paste0("PtrOp2_",PtrOp2_b$V1)
write.table(PtrOp2_b, file = gzfile("barcodes.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

### Tension Bio2
PtrTens2 <- read.table("PtrTens2_filtered_feature_bc_matrix/features.tsv.gz",sep="\t")
colnames(PtrTens2) <- c("geneName","geneName_2","gene_expression")

Ortho_table <- readRDS("RDS_ortho_cluster_table.rds")
Ortho_table_Ptr <- Ortho_table[Ortho_table$speciesName=="PoT",]

PtrTens2_ortho <- merge(PtrTens2, Ortho_table_Ptr, by = "geneName")
stopifnot(nrow(PtrTens2_ortho) == nrow(PtrTens2))


PtrTens2_ortho$clusterID_2 <- paste0("Cluster_",PtrTens2_ortho$clusterID)
PtrTens2_ortho$clusterID <- PtrTens2_ortho$clusterID_2
write.table(PtrTens2_ortho[,c("clusterID","clusterID_2","gene_expression")], file = gzfile("features.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

PtrTens2_b <- read.table("PtrTens2_filtered_feature_bc_matrix/barcodes.tsv.gz",sep="\t")
PtrTens2_b$V1 <- paste0("PtrTens2_",PtrTens2_b$V1)
write.table(PtrTens2_b, file = gzfile("barcodes.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

############################################################ Step 1. integration - Opposite Bio2 vs Ath

Ptr.data = Read10X(data.dir = 'PtrOp2_filtered_feature_bc_matrix')
Ptr = CreateSeuratObject(counts = Ptr.data, project = 'PtrOp2', min.cells = 3, min.features = 200)
Ptr

Ath2021L.data = Read10X(data.dir = 'Ath_2021L_filtered_feature_bc_matrix')
Ath2021L = CreateSeuratObject(counts = Ath2021L.data, project = 'Ath2021L', min.cells = 3, min.features = 200)
Ath2021L


Ath2019R1.data = Read10X(data.dir = 'Ath_2019R1_filtered_feature_bc_matrix')
Ath2019R2.data = Read10X(data.dir = 'Ath_2019R2_filtered_feature_bc_matrix')
Ath2019R3.data = Read10X(data.dir = 'Ath_2019R3_filtered_feature_bc_matrix')
Ath2019R1 = CreateSeuratObject(counts = Ath2019R1.data, project = 'Ath2019R1', min.cells = 3, min.features = 200)
Ath2019R2 = CreateSeuratObject(counts = Ath2019R2.data, project = 'Ath2019R2', min.cells = 3, min.features = 200)
Ath2019R3 = CreateSeuratObject(counts = Ath2019R3.data, project = 'Ath2019R3', min.cells = 3, min.features = 200)
Ath2019R1
Ath2019R2
Ath2019R3

#Normalize the data
Ptr = NormalizeData(Ptr, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2021L = NormalizeData(Ath2021L, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2019R1 = NormalizeData(Ath2019R1, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2019R2 = NormalizeData(Ath2019R2, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2019R3 = NormalizeData(Ath2019R3, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
Ptr = FindVariableFeatures(Ptr, selection.method = 'vst', nfeatures = 2000)
Ath2021L = FindVariableFeatures(Ath2021L, selection.method = 'vst', nfeatures = 2000)
Ath2019R1 = FindVariableFeatures(Ath2019R1, selection.method = 'vst', nfeatures = 2000)
Ath2019R2 = FindVariableFeatures(Ath2019R2, selection.method = 'vst', nfeatures = 2000)
Ath2019R3 = FindVariableFeatures(Ath2019R3, selection.method = 'vst', nfeatures = 2000)

#Identify the integration anchors
integration_anchors_Ath2021L = FindIntegrationAnchors(object.list = list(Ptr,Ath2021L),
                                                      anchor.features = 2000,
                                                      scale = TRUE,
                                                      reduction = 'cca',
                                                      l2.norm = TRUE,
                                                      k.anchor = 5)
integration_anchors_Ath2019R = FindIntegrationAnchors(object.list = list(Ptr,Ath2019R1,Ath2019R2,Ath2019R3),
                                                      anchor.features = 2000,
                                                      scale = TRUE,
                                                      reduction = 'cca',
                                                      l2.norm = TRUE,
                                                      k.anchor = 5)

Combined_object_Ath2021L = IntegrateData(anchorset = integration_anchors_Ath2021L)
Combined_object_Ath2019R = IntegrateData(anchorset = integration_anchors_Ath2019R,
                                         preserve.order = T,
                                         sample.tree = rbind(c(-1,-2),
                                                             c(1,-3),
                                                             c(2,-4)))


#Run the standard workflow for visualization and clustering
Combined_object_Ath2021L %<>% ScaleData %>% RunPCA(npcs=30) %>% RunUMAP(reduction='pca', dims=1:30)
Combined_object_Ath2021L %<>% FindNeighbors(reduction='pca', dims=1:30, k.param=3) %>% FindClusters(resolution=0.5)
Combined_object_Ath2021L@reductions$umap@cell.embeddings[,2] %<>% multiply_by(-1)
rotation_angle = 60 * (pi/180)
rotation_matrix = matrix(c(cos(rotation_angle),-sin(rotation_angle),
                           sin(rotation_angle),cos(rotation_angle)),nrow=2,byrow=T)
Combined_object_Ath2021L@reductions$umap@cell.embeddings %<>% 
    multiply_by_matrix(rotation_matrix) %>% 
    set_colnames(c('UMAP_1','UMAP_2'))

saveRDS(Combined_object_Ath2021L,'RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Ath2021L.rds')
#Combined_object_Ath2021L = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Ath2021L.rds')


Combined_object_Ath2019R %<>% ScaleData %>% RunPCA(npcs=30) %>% RunUMAP(reduction='pca', dims=1:30)
Combined_object_Ath2019R %<>% FindNeighbors(reduction='pca', dims=1:30, k.param=3) %>% FindClusters(resolution=0.5)
Combined_object_Ath2019R@reductions$umap@cell.embeddings[,1] %<>% multiply_by(-1)

saveRDS(Combined_object_Ath2019R,'RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Ath2019R.rds')
#Combined_object_Ath2019R = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Ath2019R.rds')

############################################################ Step 2. Generating RDS 



getMSTsubtreeCenter = function(projection, sp1, sp2){
    print('Calculate the distance between each pair of cells')
    distMatrix = as.matrix(dist(projection)) %>% Matrix(sparse=T)
    stopifnot(sum(distMatrix==0) == nrow(projection))

    print('Create the graph from adjacent matrix')
    graphFull = graph_from_adjacency_matrix(distMatrix,
                                            mode='undirected',weighted=T)

    print('Construct the MST')
    graphMST = mst(graphFull)

    print('Remove inter-species edges')
    edgeVname = attr(E(graphMST),'vnames')
    delEdge = edgeVname %>% grep(sp1,.,value=T) %>% grep(sp2,.,value=T)
    graphCutMST = delete_edges(graphMST,delEdge)

    print('Extract the subgraph centers')
    subgraphCenter = c()
    candidateVertices = attr(V(graphCutMST),'name')
    while(length(candidateVertices)>0){
        pickedVertex = candidateVertices[1]
        pickedVertices = attr(subcomponent(graphCutMST,pickedVertex),'name')
        pickedGraph = induced_subgraph(graphCutMST,pickedVertices)
        pickedCloseness = closeness(pickedGraph)
        if(length(pickedCloseness)==1){
            subgraphCenter %<>% c(names(pickedCloseness))
        }else{
            subgraphCenter %<>% c(names(which.max(pickedCloseness)))
        }
        candidateVertices %<>% setdiff(pickedVertices)
    }
    return(subgraphCenter)
}



runUMAPandSaveSubtreeCenter <- function(rdsFilePath, prefix, sp1, sp2) {
    # Read the RDS file
    ClaPAIR_combinedObject = readRDS(rdsFilePath)

    # Run UMAP
    ClaPAIR_combinedObject <- RunUMAP(
        object = ClaPAIR_combinedObject,
        reduction = "pca",
        dims = 1:30,
        seed.use = 42,
        min.dist = 0.3, # 0.3 # c(0.001, 0.5)
        n.neighbors = 30, # 30L # c(5, 50)
        umap.method = "uwot", metric = "cosine"
    )

    # Get UMAP projections
    ClaPAIR_projectionUMAP = ClaPAIR_combinedObject@reductions$umap@cell.embeddings

    print('Calculating MST and saving subtree center')
    ClaPAIR_subtreeCenter = getMSTsubtreeCenter(ClaPAIR_projectionUMAP, sp1, sp2)
    saveRDS(ClaPAIR_subtreeCenter, paste0('RDS_', prefix, '_subtreeCenter.rds'))
}


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Ath2021L.rds',
    prefix = 'PtrOp2Ath2021L',
    sp1 = "PtrOp2_",
    sp2 = "Ath_2021L"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Ath2019R.rds',
    prefix = 'PtrOp2Ath2019R',
    sp1 = "PtrOp2_",
    sp2 = "Ath_2019R"
)
############################################################ Step 3.1 overlapping: pie - Opposite Bio2 vs Ath
########## PtrOp2Ath2021L
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Ath2021L.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp2Ath2021L_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp2Ath2021L',
    plot_name = 'PtrOp2Ath2021L',
    sp1 = "PtrOp2_",
    df_path = paste0('Overlap_heatpie_', "PtrOp2Ath2021L", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp2Ath2021L_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrOp2Ath2021L_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)



########## PtrOp2Ath2019R
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Ath2019R.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp2Ath2019R_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp2Ath2019R',
    plot_name = 'PtrOp2Ath2019R',
    sp1 = "PtrOp2_",
    df_path = paste0('Overlap_heatpie_', "PtrOp2Ath2019R", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp2Ath2019R_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrOp2Ath2019R_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)

############################################################ Step 3.2 overlapping: map - Opposite Bio2 vs Ath


########## PtrOp2Ath2021L
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Ath2021L.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp2Ath2021L_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp2Ath2021L',
    plot_name = 'PtrOp2Ath2021L',
    sp1 = "PtrOp2_"
    )



########## PtrOp2Ath2019R
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Ath2019R.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp2Ath2019R_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp2Ath2019R',
    plot_name = 'PtrOp2Ath2019R',
    sp1 = "PtrOp2_"
)


############################################################
############################################################
############################################################
############################################################

############################################################ Step 1. integration - Tension Bio2 vs Ath

Ptr.data = Read10X(data.dir = 'PtrTens2_filtered_feature_bc_matrix')
Ptr = CreateSeuratObject(counts = Ptr.data, project = 'PtrTens2', min.cells = 3, min.features = 200)
Ptr

Ath2021L.data = Read10X(data.dir = 'Ath_2021L_filtered_feature_bc_matrix')
Ath2021L = CreateSeuratObject(counts = Ath2021L.data, project = 'Ath2021L', min.cells = 3, min.features = 200)
Ath2021L


Ath2019R1.data = Read10X(data.dir = 'Ath_2019R1_filtered_feature_bc_matrix')
Ath2019R2.data = Read10X(data.dir = 'Ath_2019R2_filtered_feature_bc_matrix')
Ath2019R3.data = Read10X(data.dir = 'Ath_2019R3_filtered_feature_bc_matrix')
Ath2019R1 = CreateSeuratObject(counts = Ath2019R1.data, project = 'Ath2019R1', min.cells = 3, min.features = 200)
Ath2019R2 = CreateSeuratObject(counts = Ath2019R2.data, project = 'Ath2019R2', min.cells = 3, min.features = 200)
Ath2019R3 = CreateSeuratObject(counts = Ath2019R3.data, project = 'Ath2019R3', min.cells = 3, min.features = 200)
Ath2019R1
Ath2019R2
Ath2019R3

#Normalize the data
Ptr = NormalizeData(Ptr, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2021L = NormalizeData(Ath2021L, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2019R1 = NormalizeData(Ath2019R1, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2019R2 = NormalizeData(Ath2019R2, normalization.method = 'LogNormalize', scale.factor = 10000)
Ath2019R3 = NormalizeData(Ath2019R3, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
Ptr = FindVariableFeatures(Ptr, selection.method = 'vst', nfeatures = 2000)
Ath2021L = FindVariableFeatures(Ath2021L, selection.method = 'vst', nfeatures = 2000)
Ath2019R1 = FindVariableFeatures(Ath2019R1, selection.method = 'vst', nfeatures = 2000)
Ath2019R2 = FindVariableFeatures(Ath2019R2, selection.method = 'vst', nfeatures = 2000)
Ath2019R3 = FindVariableFeatures(Ath2019R3, selection.method = 'vst', nfeatures = 2000)

#Identify the integration anchors
integration_anchors_Ath2021L = FindIntegrationAnchors(object.list = list(Ptr,Ath2021L),
                                                      anchor.features = 2000,
                                                      scale = TRUE,
                                                      reduction = 'cca',
                                                      l2.norm = TRUE,
                                                      k.anchor = 5)
integration_anchors_Ath2019R = FindIntegrationAnchors(object.list = list(Ptr,Ath2019R1,Ath2019R2,Ath2019R3),
                                                      anchor.features = 2000,
                                                      scale = TRUE,
                                                      reduction = 'cca',
                                                      l2.norm = TRUE,
                                                      k.anchor = 5)

Combined_object_Ath2021L = IntegrateData(anchorset = integration_anchors_Ath2021L)
Combined_object_Ath2019R = IntegrateData(anchorset = integration_anchors_Ath2019R,
                                         preserve.order = T,
                                         sample.tree = rbind(c(-1,-2),
                                                             c(1,-3),
                                                             c(2,-4)))


#Run the standard workflow for visualization and clustering
Combined_object_Ath2021L %<>% ScaleData %>% RunPCA(npcs=30) %>% RunUMAP(reduction='pca', dims=1:30)
Combined_object_Ath2021L %<>% FindNeighbors(reduction='pca', dims=1:30, k.param=3) %>% FindClusters(resolution=0.5)
Combined_object_Ath2021L@reductions$umap@cell.embeddings[,2] %<>% multiply_by(-1)
rotation_angle = 60 * (pi/180)
rotation_matrix = matrix(c(cos(rotation_angle),-sin(rotation_angle),
                           sin(rotation_angle),cos(rotation_angle)),nrow=2,byrow=T)
Combined_object_Ath2021L@reductions$umap@cell.embeddings %<>% 
    multiply_by_matrix(rotation_matrix) %>% 
    set_colnames(c('UMAP_1','UMAP_2'))

saveRDS(Combined_object_Ath2021L,'RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Ath2021L.rds')
#Combined_object_Ath2021L = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Ath2021L.rds')


Combined_object_Ath2019R %<>% ScaleData %>% RunPCA(npcs=30) %>% RunUMAP(reduction='pca', dims=1:30)
Combined_object_Ath2019R %<>% FindNeighbors(reduction='pca', dims=1:30, k.param=3) %>% FindClusters(resolution=0.5)
Combined_object_Ath2019R@reductions$umap@cell.embeddings[,1] %<>% multiply_by(-1)

saveRDS(Combined_object_Ath2019R,'RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Ath2019R.rds')
#Combined_object_Ath2019R = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Ath2019R.rds')

############################################################ Step 2. Generating RDS - Tension Bio2 vs Ath 



getMSTsubtreeCenter = function(projection, sp1, sp2){
    print('Calculate the distance between each pair of cells')
    distMatrix = as.matrix(dist(projection)) %>% Matrix(sparse=T)
    stopifnot(sum(distMatrix==0) == nrow(projection))

    print('Create the graph from adjacent matrix')
    graphFull = graph_from_adjacency_matrix(distMatrix,
                                            mode='undirected',weighted=T)

    print('Construct the MST')
    graphMST = mst(graphFull)

    print('Remove inter-species edges')
    edgeVname = attr(E(graphMST),'vnames')
    delEdge = edgeVname %>% grep(sp1,.,value=T) %>% grep(sp2,.,value=T)
    graphCutMST = delete_edges(graphMST,delEdge)

    print('Extract the subgraph centers')
    subgraphCenter = c()
    candidateVertices = attr(V(graphCutMST),'name')
    while(length(candidateVertices)>0){
        pickedVertex = candidateVertices[1]
        pickedVertices = attr(subcomponent(graphCutMST,pickedVertex),'name')
        pickedGraph = induced_subgraph(graphCutMST,pickedVertices)
        pickedCloseness = closeness(pickedGraph)
        if(length(pickedCloseness)==1){
            subgraphCenter %<>% c(names(pickedCloseness))
        }else{
            subgraphCenter %<>% c(names(which.max(pickedCloseness)))
        }
        candidateVertices %<>% setdiff(pickedVertices)
    }
    return(subgraphCenter)
}



runUMAPandSaveSubtreeCenter <- function(rdsFilePath, prefix, sp1, sp2) {
    # Read the RDS file
    ClaPAIR_combinedObject = readRDS(rdsFilePath)

    # Run UMAP
    ClaPAIR_combinedObject <- RunUMAP(
        object = ClaPAIR_combinedObject,
        reduction = "pca",
        dims = 1:30,
        seed.use = 42,
        min.dist = 0.3, # 0.3 # c(0.001, 0.5)
        n.neighbors = 30, # 30L # c(5, 50)
        umap.method = "uwot", metric = "cosine"
    )

    # Get UMAP projections
    ClaPAIR_projectionUMAP = ClaPAIR_combinedObject@reductions$umap@cell.embeddings

    print('Calculating MST and saving subtree center')
    ClaPAIR_subtreeCenter = getMSTsubtreeCenter(ClaPAIR_projectionUMAP, sp1, sp2)
    saveRDS(ClaPAIR_subtreeCenter, paste0('RDS_', prefix, '_subtreeCenter.rds'))
}


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Ath2021L.rds',
    prefix = 'PtrTens2Ath2021L',
    sp1 = "PtrTens2_",
    sp2 = "Ath_2021L"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Ath2019R.rds',
    prefix = 'PtrTens2Ath2019R',
    sp1 = "PtrTens2_",
    sp2 = "Ath_2019R"
)



############################################################ Step 3.1 overlapping: pie - Tension Bio2 vs Ath
########## PtrTens2Ath2021L
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Ath2021L.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens2Ath2021L_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens2Ath2021L',
    plot_name = 'PtrTens2Ath2021L',
    sp1 = "PtrTens2_",
    df_path = paste0('Overlap_heatpie_', "PtrTens2Ath2021L", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrTens2Ath2021L_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrTens2Ath2021L_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)



########## PtrTens2Ath2019R
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Ath2019R.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens2Ath2019R_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens2Ath2019R',
    plot_name = 'PtrTens2Ath2019R',
    sp1 = "PtrTens2_",
    df_path = paste0('Overlap_heatpie_', "PtrTens2Ath2019R", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrTens2Ath2019R_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrTens2Ath2019R_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)

############################################################ Step 3.2 overlapping: map - Tension Bio2 vs Ath

########## PtrTens2Ath2021L
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Ath2021L.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens2Ath2021L_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens2Ath2021L',
    plot_name = 'PtrTens2Ath2021L',
    sp1 = "PtrTens2_"
    )



########## PtrTens2Ath2019R
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Ath2019R.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens2Ath2019R_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens2Ath2019R',
    plot_name = 'PtrTens2Ath2019R',
    sp1 = "PtrTens2_"
)


############################################################
############################################################
############################################################
############################################################
############################################################ Step 1. integration - Opposite Bio2 vs Osa

Ptr.data = Read10X(data.dir = 'PtrOp2_filtered_feature_bc_matrix')
Ptr = CreateSeuratObject(counts = Ptr.data, project = 'PtrOp2', min.cells = 3, min.features = 200)
Ptr

Osa2021R1.data = Read10X(data.dir = 'Osa_2021R1_filtered_feature_bc_matrix')
Osa2021R2.data = Read10X(data.dir = 'Osa_2021R2_filtered_feature_bc_matrix')
Osa2021R1 = CreateSeuratObject(counts = Osa2021R1.data, project = 'Osa2021R1', min.cells = 3, min.features = 200)
Osa2021R2 = CreateSeuratObject(counts = Osa2021R2.data, project = 'Osa2021R2', min.cells = 3, min.features = 200)
Osa2021R1
Osa2021R2

#Normalize the data
Ptr = NormalizeData(Ptr, normalization.method = 'LogNormalize', scale.factor = 10000)
Osa2021R1 = NormalizeData(Osa2021R1, normalization.method = 'LogNormalize', scale.factor = 10000)
Osa2021R2 = NormalizeData(Osa2021R2, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
Ptr = FindVariableFeatures(Ptr, selection.method = 'vst', nfeatures = 2000)
Osa2021R1 = FindVariableFeatures(Osa2021R1, selection.method = 'vst', nfeatures = 2000)
Osa2021R2 = FindVariableFeatures(Osa2021R2, selection.method = 'vst', nfeatures = 2000)

#Identify the integration anchors
integration_anchors_Osa2021R = FindIntegrationAnchors(object.list = list(Ptr,Osa2021R1,Osa2021R2),
                                                      anchor.features = 2000,
                                                      scale = TRUE,
                                                      reduction = 'cca',
                                                      l2.norm = TRUE,
                                                      k.anchor = 5)

Combined_object_Osa2021R = IntegrateData(anchorset = integration_anchors_Osa2021R,
                                         preserve.order = T,
                                         sample.tree = rbind(c(-3,-2),
                                                             c(-1,1)))


#Run the standard workflow for visualization and clustering
Combined_object_Osa2021R %<>% ScaleData %>% RunPCA(npcs=30) %>% RunUMAP(reduction='pca', dims=1:30)
Combined_object_Osa2021R %<>% FindNeighbors(reduction='pca', dims=1:30, k.param=3) %>% FindClusters(resolution=0.5)

# Combined_object_Osa2021R@reductions$umap@cell.embeddings[,2] %<>% multiply_by(-1)
rotation_angle = -45 * (pi/180)
rotation_matrix = matrix(c(cos(rotation_angle),-sin(rotation_angle),
                           sin(rotation_angle),cos(rotation_angle)),nrow=2,byrow=T)
Combined_object_Osa2021R@reductions$umap@cell.embeddings %<>%
    multiply_by_matrix(rotation_matrix) %>%
    set_colnames(c('UMAP_1','UMAP_2'))

saveRDS(Combined_object_Osa2021R,'RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Osa2021R.rds')
# Combined_object_Osa2021R = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Osa2021R.rds')


############################################################ Step 1. integration - Tension Bio2 vs Osa

Ptr.data = Read10X(data.dir = 'PtrTens2_filtered_feature_bc_matrix')
Ptr = CreateSeuratObject(counts = Ptr.data, project = 'PtrTens2', min.cells = 3, min.features = 200)
Ptr

Osa2021R1.data = Read10X(data.dir = 'Osa_2021R1_filtered_feature_bc_matrix')
Osa2021R2.data = Read10X(data.dir = 'Osa_2021R2_filtered_feature_bc_matrix')
Osa2021R1 = CreateSeuratObject(counts = Osa2021R1.data, project = 'Osa2021R1', min.cells = 3, min.features = 200)
Osa2021R2 = CreateSeuratObject(counts = Osa2021R2.data, project = 'Osa2021R2', min.cells = 3, min.features = 200)
Osa2021R1
Osa2021R2

#Normalize the data
Ptr = NormalizeData(Ptr, normalization.method = 'LogNormalize', scale.factor = 10000)
Osa2021R1 = NormalizeData(Osa2021R1, normalization.method = 'LogNormalize', scale.factor = 10000)
Osa2021R2 = NormalizeData(Osa2021R2, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
Ptr = FindVariableFeatures(Ptr, selection.method = 'vst', nfeatures = 2000)
Osa2021R1 = FindVariableFeatures(Osa2021R1, selection.method = 'vst', nfeatures = 2000)
Osa2021R2 = FindVariableFeatures(Osa2021R2, selection.method = 'vst', nfeatures = 2000)

#Identify the integration anchors
integration_anchors_Osa2021R = FindIntegrationAnchors(object.list = list(Ptr,Osa2021R1,Osa2021R2),
                                                      anchor.features = 2000,
                                                      scale = TRUE,
                                                      reduction = 'cca',
                                                      l2.norm = TRUE,
                                                      k.anchor = 5)

Combined_object_Osa2021R = IntegrateData(anchorset = integration_anchors_Osa2021R,
                                         preserve.order = T,
                                         sample.tree = rbind(c(-3,-2),
                                                             c(-1,1)))


#Run the standard workflow for visualization and clustering
Combined_object_Osa2021R %<>% ScaleData %>% RunPCA(npcs=30) %>% RunUMAP(reduction='pca', dims=1:30)
Combined_object_Osa2021R %<>% FindNeighbors(reduction='pca', dims=1:30, k.param=3) %>% FindClusters(resolution=0.5)

# Combined_object_Osa2021R@reductions$umap@cell.embeddings[,2] %<>% multiply_by(-1)
rotation_angle = -45 * (pi/180)
rotation_matrix = matrix(c(cos(rotation_angle),-sin(rotation_angle),
                           sin(rotation_angle),cos(rotation_angle)),nrow=2,byrow=T)
Combined_object_Osa2021R@reductions$umap@cell.embeddings %<>%
    multiply_by_matrix(rotation_matrix) %>%
    set_colnames(c('UMAP_1','UMAP_2'))

saveRDS(Combined_object_Osa2021R,'RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Osa2021R.rds')
# Combined_object_Osa2021R = readRDS('RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Osa2021R.rds')


############################################################ Step 2. Generating RDS 

############## subtree (due to limited memory, use tapir)
library(igraph)
library(magrittr)
library(Matrix)
library(MASS)
library(Seurat)
library(ggplot2)
#/work4/home/joweihsieh/miniconda3
setwd("/work4/home/joweihsieh/test/20240215_tens")

oriPar = par(no.readonly=T)


getMSTsubtreeCenter = function(projection, sp1, sp2){
    print('Calculate the distance between each pair of cells')
    distMatrix = as.matrix(dist(projection)) %>% Matrix(sparse=T)
    stopifnot(sum(distMatrix==0) == nrow(projection))

    print('Create the graph from adjacent matrix')
    graphFull = graph_from_adjacency_matrix(distMatrix,
                                            mode='undirected',weighted=T)

    print('Construct the MST')
    graphMST = mst(graphFull)

    print('Remove inter-species edges')
    edgeVname = attr(E(graphMST),'vnames')
    delEdge = edgeVname %>% grep(sp1,.,value=T) %>% grep(sp2,.,value=T)
    graphCutMST = delete_edges(graphMST,delEdge)

    print('Extract the subgraph centers')
    subgraphCenter = c()
    candidateVertices = attr(V(graphCutMST),'name')
    while(length(candidateVertices)>0){
        pickedVertex = candidateVertices[1]
        pickedVertices = attr(subcomponent(graphCutMST,pickedVertex),'name')
        pickedGraph = induced_subgraph(graphCutMST,pickedVertices)
        pickedCloseness = closeness(pickedGraph)
        if(length(pickedCloseness)==1){
            subgraphCenter %<>% c(names(pickedCloseness))
        }else{
            subgraphCenter %<>% c(names(which.max(pickedCloseness)))
        }
        candidateVertices %<>% setdiff(pickedVertices)
    }
    return(subgraphCenter)
}



runUMAPandSaveSubtreeCenter <- function(rdsFilePath, prefix, sp1, sp2) {
    # Read the RDS file
    ClaPAIR_combinedObject = readRDS(rdsFilePath)

    # Run UMAP
    ClaPAIR_combinedObject <- RunUMAP(
        object = ClaPAIR_combinedObject,
        reduction = "pca",
        dims = 1:30,
        seed.use = 42,
        min.dist = 0.3, # 0.3 # c(0.001, 0.5)
        n.neighbors = 30, # 30L # c(5, 50)
        umap.method = "uwot", metric = "cosine"
    )

    # Get UMAP projections
    ClaPAIR_projectionUMAP = ClaPAIR_combinedObject@reductions$umap@cell.embeddings

    print('Calculating MST and saving subtree center')
    ClaPAIR_subtreeCenter = getMSTsubtreeCenter(ClaPAIR_projectionUMAP, sp1, sp2)
    saveRDS(ClaPAIR_subtreeCenter, paste0('RDS_', prefix, '_subtreeCenter.rds'))
}


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrOp2_Osa2021R.rds',
    prefix = 'PtrOp2Osa2021R',
    sp1 = "PtrOp2_",
    sp2 = "Osa_2021R"
)


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_PtrTens2_Osa2021R.rds',
    prefix = 'PtrTens2Osa2021R',
    sp1 = "PtrTens2_",
    sp2 = "Osa_2021R"
)

