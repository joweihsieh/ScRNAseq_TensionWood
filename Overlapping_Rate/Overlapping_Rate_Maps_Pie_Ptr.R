library(igraph)
library(magrittr)
library(Matrix)
library(MASS)
library(Seurat)
library(ggplot2)

setwd("/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Multi_species_analysis/all_data_rds")

oriPar = par(no.readonly=T)

############################################################################## Step 1. Generating RDS 


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
    rdsFilePath = 'integration_PtrOp12_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrOp12_rep',
    sp1 = "OW_",
    sp2 = "Opposite02_"
)



runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp13_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrOp13_rep',
    sp1 = "OW_",
    sp2 = "Opposite03_"
)



runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp14_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrOp14_rep',
    sp1 = "OW_",
    sp2 = "Opposite04_"
)


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp23_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrOp23_rep',
    sp1 = "Opposite02_",
    sp2 = "Opposite03_"
)


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp24_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrOp24_rep',
    sp1 = "Opposite02_",
    sp2 = "Opposite04_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp34_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrOp34_rep',
    sp1 = "Opposite03_",
    sp2 = "Opposite04_"
)


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens12_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrTens12_rep',
    sp1 = "TW_",
    sp2 = "Tension02_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens13_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrTens13_rep',
    sp1 = "TW_",
    sp2 = "Tension03_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens14_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrTens14_rep',
    sp1 = "TW_",
    sp2 = "Tension04_"
)


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens23_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrTens23_rep',
    sp1 = "Tension02_",
    sp2 = "Tension03_"
)


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens24_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrTens24_rep',
    sp1 = "Tension02_",
    sp2 = "Tension04_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens34_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrTens34_rep',
    sp1 = "Tension03_",
    sp2 = "Tension04_"
)



runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert12_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrVert12_rep',
    sp1 = "Vertical01_",
    sp2 = "Vertical02_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert13_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrVert13_rep',
    sp1 = "Vertical01_",
    sp2 = "Vertical03_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert14_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrVert14_rep',
    sp1 = "Vertical01_",
    sp2 = "Vertical04_"
)


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert23_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrVert23_rep',
    sp1 = "Vertical02_",
    sp2 = "Vertical03_"
)


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert24_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrVert24_rep',
    sp1 = "Vertical02_",
    sp2 = "Vertical04_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert34_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrVert34_rep',
    sp1 = "Vertical03_",
    sp2 = "Vertical04_"
)



runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert1_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrVert1',
    sp1 = "Ptr_",
    sp2 = "Vertical01_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrVert2',
    sp1 = "Ptr_",
    sp2 = "Vertical02_"
)


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert3_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrVert3',
    sp1 = "Ptr_",
    sp2 = "Vertical03_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert4_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrVert4',
    sp1 = "Ptr_",
    sp2 = "Vertical04_"
)



runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens1_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrTens1',
    sp1 = "Ptr_",
    sp2 = "TW_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrTens2',
    sp1 = "Ptr_",
    sp2 = "Tension02_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens3_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrTens3',
    sp1 = "Ptr_",
    sp2 = "Tension03_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens4_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrTens4',
    sp1 = "Ptr_",
    sp2 = "Tension04_"
)


runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp1_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrOp1',
    sp1 = "Ptr_",
    sp2 = "OW_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrOp2',
    sp1 = "Ptr_",
    sp2 = "Opposite02_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp3_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrOp3',
    sp1 = "Ptr_",
    sp2 = "Opposite03_"
)

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp4_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds',
    prefix = 'PtrOp4',
    sp1 = "Ptr_",
    sp2 = "Opposite04_"
)

############################################################################## Step 2. maps


runUMAPandSaveSubtreeCenter <- function(rdsFilePath) {
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
    return(ClaPAIR_projectionUMAP)
}



centerContourPlot = function(subgraphCenter = ClaPAIR_subtreeCenter,
                             projectionUMAP = ClaPAIR_projectionUMAP,
                             lowerPercentile = 5,
                             higherPercentile = 40,
                             main = 'PtrTens12_rep',
                             plot_name = 'PtrTens12_rep',
                             sp1 = "TW_",
                             df_path = NULL){
    
    Ptr_projectionUMAP = projectionUMAP %>% extract(grepl(sp1,rownames(.)),)
    Ptr_density_map = kde2d(Ptr_projectionUMAP[,1],
                            Ptr_projectionUMAP[,2],n=500,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))
    get_territory_density = function(Coor){
        out = Ptr_density_map$z[max(which(Ptr_density_map$x < Coor[1])),
                                max(which(Ptr_density_map$y < Coor[2]))]
        return(out)
    }
    Ptr_territory_density = Ptr_projectionUMAP %>% apply(1,get_territory_density)
    norm_factor = 1/sum(Ptr_territory_density)
    # print(norm_factor)
    
    Ptr_SubgraphCenter = subgraphCenter %>% extract(grepl(sp1,.))
    Sp2_SubgraphCenter = subgraphCenter %>% extract(!grepl(sp1,.))
    
    num_Ptr =  sum(grepl(sp1,rownames(projectionUMAP)))
    num_Sp2 =  sum(!grepl(sp1,rownames(projectionUMAP)))
    num_total = nrow(projectionUMAP)
    num_center_Ptr = length(Ptr_SubgraphCenter)
    num_center_Sp2 = length(Sp2_SubgraphCenter)

    get_subgraphcenter_density = function(partSubgraphCenter){
        allCenterCoordinate = projectionUMAP[subgraphCenter,]
        partCenterCoordinate = projectionUMAP[partSubgraphCenter,]
        density_map = kde2d(partCenterCoordinate[,1],
                            partCenterCoordinate[,2],n=500,
                            h=apply(allCenterCoordinate,2,bandwidth.nrd)/2,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))
        return(density_map)
    }
    Ptr_subgraphcenter_density = get_subgraphcenter_density(Ptr_SubgraphCenter)
    Sp2_subgraphcenter_density = get_subgraphcenter_density(Sp2_SubgraphCenter)
    stopifnot(Ptr_subgraphcenter_density$x==Sp2_subgraphcenter_density$x)
    stopifnot(Ptr_subgraphcenter_density$y==Sp2_subgraphcenter_density$y)
    
    merge_subgraphcenter_density = Ptr_subgraphcenter_density
    merge_subgraphcenter_density$z =
        (num_Ptr/num_total)*num_center_Ptr*Ptr_subgraphcenter_density$z +
        (num_Sp2/num_total)*num_center_Sp2*Sp2_subgraphcenter_density$z
    
    merge_subgraphcenter_density$z %<>% multiply_by(norm_factor)
    
    message(
      "Plot density max:",
      round(max(merge_subgraphcenter_density$z), 5)
    )
    message(
      "Plot density Q75:",
      round(quantile(merge_subgraphcenter_density$z, 0.75), 5)
    )
    message(
      "Plot density Q50:",
      round(quantile(merge_subgraphcenter_density$z, 0.50), 5)
    )
    message(
      "Plot density Q25:",
      round(quantile(merge_subgraphcenter_density$z, 0.25), 5)
    )
    message(
      "Plot density min:",
      round(min(merge_subgraphcenter_density$z), 5)
    )
    
     territory_map = kde2d(projectionUMAP[,1],
                          projectionUMAP[,2],n=500,h=0.02,
                          lims=c(min(projectionUMAP[,1])-0.5,
                                 max(projectionUMAP[,1])+0.5,
                                 min(projectionUMAP[,2])-0.5,
                                 max(projectionUMAP[,2])+0.5))
    
    png(paste0(plot_name,'.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    {
        plot(NA,
             xlim=range(projectionUMAP[,'umap_1']),
             ylim=range(projectionUMAP[,'umap_2']),
             xlab='',ylab='',axes=F,main=main)
        Levels = seq(0,1,length.out=500)
        lowerRank = 500 * lowerPercentile/100
        higherRank = 500 * higherPercentile/100
        .filled.contour(merge_subgraphcenter_density$x,
                        merge_subgraphcenter_density$y,
                        merge_subgraphcenter_density$z,
                        levels=Levels,
                        col=c(colorRampPalette(c('#AED2F5','#FECB71'))(lowerRank),
                              colorRampPalette(c('#FECB71','#F8696B'))(higherRank-lowerRank),
                              colorRampPalette(c('#F8696B','#713020'))(500-higherRank)))
         .filled.contour(territory_map$x,
                        territory_map$y,
                        ifelse(territory_map$z>0,1,0),
                        levels=c(0,0.5),col=c('white',NA))
        contour(territory_map$x,territory_map$y,ifelse(territory_map$z>0,1,0),
                levels=0.5,lwd=0.8,drawlabels=F,add=T)
    }
    dev.off()
}



##################### PtrTens12_rep
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens12_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens12_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens12_rep',
    plot_name = 'PtrTens12_rep',
    sp1 = "TW_"
    )


##################### PtrTens13_rep


ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens13_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens13_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens13_rep',
    plot_name = 'PtrTens13_rep',
    sp1 = "TW_"
    )

##################### PtrTens14_rep

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens14_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens14_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens14_rep',
    plot_name = 'PtrTens14_rep',
    sp1 = "TW_"
    )


##################### PtrTens23_rep
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens23_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens23_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens23_rep',
    plot_name = 'PtrTens23_rep',
    sp1 = "Tension02_"
    )



##################### PtrTens24_rep


ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens24_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens24_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens24_rep',
    plot_name = 'PtrTens24_rep',
    sp1 = "Tension02_"
    )




##################### PtrTens34_rep

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens34_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens34_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens34_rep',
    plot_name = 'PtrTens34_rep',
    sp1 = "Tension03_"
    )



##################### PtrOp12_rep
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp12_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp12_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp12_rep',
    plot_name = 'PtrOp12_rep',
    sp1 = "OW_"
    )


##################### PtrOp13_rep


ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp13_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp13_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp13_rep',
    plot_name = 'PtrOp13_rep',
    sp1 = "OW_"
    )



##################### PtrOp14_rep

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp14_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp14_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp14_rep',
    plot_name = 'PtrOp14_rep',
    sp1 = "OW_"
    )




##################### PtrOp23_rep
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp23_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp23_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp23_rep',
    plot_name = 'PtrOp23_rep',
    sp1 = "Opposite02_"
    )



##################### PtrOp24_rep


ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp24_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp24_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp24_rep',
    plot_name = 'PtrOp24_rep',
    sp1 = "Opposite02_"
    )



##################### PtrOp34_rep

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp34_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp34_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp34_rep',
    plot_name = 'PtrOp34_rep',
    sp1 = "Opposite03_"
    )



##################### PtrVert12_rep

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert12_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrVert12_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert12_rep',
    plot_name = 'PtrVert12_rep',
    sp1 = "Vertical01_"
    )




##################### PtrVert13_rep


ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert13_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrVert13_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert13_rep',
    plot_name = 'PtrVert13_rep',
    sp1 = "Vertical01_"
    )



##################### PtrVert14_rep

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert14_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrVert14_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert14_rep',
    plot_name = 'PtrVert14_rep',
    sp1 = "Vertical01_"
    )


##################### PtrVert23_rep
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert23_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrVert23_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert23_rep',
    plot_name = 'PtrVert23_rep',
    sp1 = "Vertical02_"
    )


##################### PtrVert24_rep


ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert24_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrVert24_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert24_rep',
    plot_name = 'PtrVert24_rep',
    sp1 = "Vertical02_"
    )



##################### PtrVert34_rep

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert34_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrVert34_rep_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert34_rep',
    plot_name = 'PtrVert34_rep',
    sp1 = "Vertical03_"
    )



##################### PtrVert1
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert1_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrVert1_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert1',
    plot_name = 'PtrVert1',
    sp1 = "Ptr_"
    )


##################### PtrVert2
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrVert2_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert2',
    plot_name = 'PtrVert2',
    sp1 = "Ptr_"
    )


##################### PtrVert3
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert3_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrVert3_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert3',
    plot_name = 'PtrVert3',
    sp1 = "Ptr_",
    )


##################### PtrVert4
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert4_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrVert4_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert4',
    plot_name = 'PtrVert4',
    sp1 = "Ptr_",
    )

#

##################### PtrTens1

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens1_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens1_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens1',
    plot_name = 'PtrTens1',
    sp1 = "Ptr_"
    )

##################### PtrTens2

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens2_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens2',
    plot_name = 'PtrTens2',
    sp1 = "Ptr_"
    )


##################### PtrTens3

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens3_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens3_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens3',
    plot_name = 'PtrTens3',
    sp1 = "Ptr_"
    )

##################### PtrTens4

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens4_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens4_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens4',
    plot_name = 'PtrTens4',
    sp1 = "Ptr_"
    )

##################### PtrOp1

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp1_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp1_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp1',
    plot_name = 'PtrOp1',
    sp1 = "Ptr_"
    )


##################### PtrOp2

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp2_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp2',
    plot_name = 'PtrOp2',
    sp1 = "Ptr_"
    )


##################### PtrOp3

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp3_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp3_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp3',
    plot_name = 'PtrOp3',
    sp1 = "Ptr_"
    )

##################### PtrOp4

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp4_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp4_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp4',
    plot_name = 'PtrOp4',
    sp1 = "Ptr_"
    )





############################################################################## Step 3. pie



centerContourPlot_pie = function(subgraphCenter = ClaPAIR_subtreeCenter,
                             projectionUMAP = ClaPAIR_projectionUMAP,
                             lowerPercentile = 5,
                             higherPercentile = 40,
                             main = 'PtrVert2',
                             plot_name = 'PtrVert2',
                             sp1 = "Ptr_",
                             df_path = paste0('Overlap_heatpie_', "PtrVert2", '_on_Ptr.csv')){
    
    Ptr_projectionUMAP = projectionUMAP %>% extract(grepl(sp1,rownames(.)),)
    Ptr_density_map = kde2d(Ptr_projectionUMAP[,1],
                            Ptr_projectionUMAP[,2],n=500,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))
    get_territory_density = function(Coor){
        out = Ptr_density_map$z[max(which(Ptr_density_map$x < Coor[1])),
                                max(which(Ptr_density_map$y < Coor[2]))]
        return(out)
    }
    Ptr_territory_density = Ptr_projectionUMAP %>% apply(1,get_territory_density)
    norm_factor = 1/sum(Ptr_territory_density)
    # print(norm_factor)
    
    Ptr_SubgraphCenter = subgraphCenter %>% extract(grepl(sp1,.))
    Sp2_SubgraphCenter = subgraphCenter %>% extract(!grepl(sp1,.))
    
    num_Ptr =  sum(grepl(sp1,rownames(projectionUMAP)))
    num_Sp2 =  sum(!grepl(sp1,rownames(projectionUMAP)))
    num_total = nrow(projectionUMAP)
    num_center_Ptr = length(Ptr_SubgraphCenter)
    num_center_Sp2 = length(Sp2_SubgraphCenter)

    get_subgraphcenter_density = function(partSubgraphCenter){
        allCenterCoordinate = projectionUMAP[subgraphCenter,]
        partCenterCoordinate = projectionUMAP[partSubgraphCenter,]
        density_map = kde2d(partCenterCoordinate[,1],
                            partCenterCoordinate[,2],n=500,
                            h=apply(allCenterCoordinate,2,bandwidth.nrd)/2,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))
        return(density_map)
    }
    Ptr_subgraphcenter_density = get_subgraphcenter_density(Ptr_SubgraphCenter)
    Sp2_subgraphcenter_density = get_subgraphcenter_density(Sp2_SubgraphCenter)
    stopifnot(Ptr_subgraphcenter_density$x==Sp2_subgraphcenter_density$x)
    stopifnot(Ptr_subgraphcenter_density$y==Sp2_subgraphcenter_density$y)
    
    merge_subgraphcenter_density = Ptr_subgraphcenter_density
    merge_subgraphcenter_density$z =
        (num_Ptr/num_total)*num_center_Ptr*Ptr_subgraphcenter_density$z +
        (num_Sp2/num_total)*num_center_Sp2*Sp2_subgraphcenter_density$z
    
    merge_subgraphcenter_density$z %<>% multiply_by(norm_factor)
    
    
    territory_map = kde2d(projectionUMAP[,1],
                          projectionUMAP[,2],n=500,h=0.02,
                          lims=c(min(projectionUMAP[,1])-0.5,
                                 max(projectionUMAP[,1])+0.5,
                                 min(projectionUMAP[,2])-0.5,
                                 max(projectionUMAP[,2])+0.5))
    
    message(
      "Cell region density max:",
      round(max(merge_subgraphcenter_density$z[territory_map$z > 0]), 5)
    )
    message(
      "Cell region density Q75:",
      round(quantile(
        merge_subgraphcenter_density$z[territory_map$z > 0], 0.75), 5)
    )
    message(
      "Cell region density Q50:",
      round(quantile(
        merge_subgraphcenter_density$z[territory_map$z > 0], 0.50), 5)
    )
    message(
      "Cell region density Q25:",
      round(quantile(
        merge_subgraphcenter_density$z[territory_map$z > 0], 0.25), 5)
    )
    message(
      "Cell region density min:",
      round(min(merge_subgraphcenter_density$z[territory_map$z > 0]), 5)
    )

    {
        n_col_levels <- 500
        lower_rank <- n_col_levels * lowerPercentile / 100
        higher_rank <- n_col_levels * higherPercentile / 100
        filled_col <-
            c(colorRampPalette(c('#AED2F5','#FECB71'))(lower_rank),
              colorRampPalette(c('#FECB71','#F8696B'))(higher_rank-lower_rank),
              colorRampPalette(c('#F8696B','#713020'))(500-higher_rank))
        col_counts <-
            cut(
                merge_subgraphcenter_density$z[territory_map$z > 0],
                breaks = seq(0, 1, length.out = n_col_levels + 1)
            ) %>%
            table() %>%
            as.vector()
        col_counts_df <- data.frame(filled_col, col_counts)
#Densities from 0 to 1 are divided into 500 bins with different color shading, with proportions of different densities shown in a pie chart in each panel.  
#remove #B4D1EA  
#total is 60   
#sum(col_counts_df$col_counts)
#75898
#(75898-60)/75898

        agg_filled_col <- sapply(
            seq(n_col_levels / 5),
            function(i) {
                out <- col_counts_df$filled_col[(i - 1) * 5 + 3]
                return(out)
            }
        )
        agg_col_counts <- sapply(
            seq(n_col_levels / 5),
            function(i) {
                out <- sum(col_counts_df$col_counts[1:5 + (i - 1) * 5])
                return(out)
            }
        )
        col_counts_df <- data.frame(
            filled_col = agg_filled_col,
            col_counts = agg_col_counts
        )
        
        col_counts_df$col_props <-
            col_counts_df$col_counts / sum(col_counts_df$col_counts)
        
        if (!is.null(df_path)) write.csv(col_counts_df, df_path, row.names = FALSE)
        
        col_counts_df$factor_filled_col <-
            factor(col_counts_df$filled_col, levels = col_counts_df$filled_col)
        ggplot(
            col_counts_df,
            aes(x = "", y = col_props, fill = factor_filled_col)
        ) +
            geom_bar(stat = "identity", colour = "white", size = 0.05) +
            scale_fill_manual("legend", values = setNames(filled_col, filled_col)) +
            coord_polar("y", start = 0) +
            theme_void() +
            theme(legend.position="none")
        ggsave(
            paste0("Heathist_", plot_name, ".png"),
            width = 20, height = 15, units = "cm", dpi = 300
        )

    }
}



##################### PtrVert1
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = '../integration_PtrVert1_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('../RDS_PtrVert1_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert1',
    plot_name = 'PtrVert1',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrVert1", '_on_Ptr.csv')
    )

ovelaps <- read.csv("Overlap_heatpie_PtrVert1_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100

#94.25304
#94.4923
##################### PtrVert2
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = '../integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('../RDS_PtrVert2_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert2',
    plot_name = 'PtrVert2',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrVert2", '_on_Ptr.csv')
    )

#

#a = readRDS("integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds")
#b = readRDS("RDS_PtrVert2_subtreeCenter.rds")


ovelaps <- read.csv("Overlap_heatpie_PtrVert2_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100
#99.52444575
percent
##################### PtrVert3
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = '../integration_PtrVert3_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('../RDS_PtrVert3_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert3',
    plot_name = 'PtrVert3',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrVert3", '_on_Ptr.csv')
    )

#

#a = readRDS("integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds")
#b = readRDS("RDS_PtrVert2_subtreeCenter.rds")


ovelaps <- read.csv("Overlap_heatpie_PtrVert3_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100
#99.03562

percent
##################### PtrVert4
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = '../integration_PtrVert4_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('../RDS_PtrVert4_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert4',
    plot_name = 'PtrVert4',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrVert4", '_on_Ptr.csv')
    )

#

#a = readRDS("integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds")
#b = readRDS("RDS_PtrVert2_subtreeCenter.rds")


ovelaps <- read.csv("Overlap_heatpie_PtrVert4_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100
#99.88212
percent
##################### PtrVert12
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = '../integration_PtrVert12_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('../RDS_PtrVert12_rep_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert12_rep',
    plot_name = 'PtrVert12_rep',
    sp1 = "PtrVertical01_",
    df_path = paste0('Overlap_heatpie_', "PtrVert12_rep", '_on_Ptr.csv')
    )

#

#a = readRDS("integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds")
#b = readRDS("RDS_PtrVert2_subtreeCenter.rds")


ovelaps <- read.csv("Overlap_heatpie_PtrVert12_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


#99.69081
#99.5536
##################### PtrVert13
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = '../integration_PtrVert13_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('../RDS_PtrVert13_rep_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert13_rep',
    plot_name = 'PtrVert13_rep',
    sp1 = "PtrVertical01_",
    df_path = paste0('Overlap_heatpie_', "PtrVert13_rep", '_on_Ptr.csv')
    )

#

#a = readRDS("integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds")
#b = readRDS("RDS_PtrVert2_subtreeCenter.rds")


ovelaps <- read.csv("Overlap_heatpie_PtrVert13_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100
#99.18869


##################### PtrVert14
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = '../integration_PtrVert14_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('../RDS_PtrVert14_rep_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert14_rep',
    plot_name = 'PtrVert14_rep',
    sp1 = "PtrVertical01_",
    df_path = paste0('Overlap_heatpie_', "PtrVert14_rep", '_on_Ptr.csv')
    )

#

#a = readRDS("integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds")
#b = readRDS("RDS_PtrVert2_subtreeCenter.rds")


ovelaps <- read.csv("Overlap_heatpie_PtrVert14_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100
percent
#99.41958
##################### PtrVert23
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = '../integration_PtrVert23_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('../RDS_PtrVert23_rep_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert23_rep',
    plot_name = 'PtrVert23_rep',
    sp1 = "PtrVertical02_",
    df_path = paste0('Overlap_heatpie_', "PtrVert23_rep", '_on_Ptr.csv')
    )

#

#a = readRDS("integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds")
#b = readRDS("RDS_PtrVert2_subtreeCenter.rds")


ovelaps <- read.csv("Overlap_heatpie_PtrVert23_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100

percent
#99.90961
##################### PtrVert24
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = '../integration_PtrVert24_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('../RDS_PtrVert24_rep_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert24_rep',
    plot_name = 'PtrVert24_rep',
    sp1 = "PtrVertical02_",
    df_path = paste0('Overlap_heatpie_', "PtrVert24_rep", '_on_Ptr.csv')
    )

#

#a = readRDS("integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds")
#b = readRDS("RDS_PtrVert2_subtreeCenter.rds")


ovelaps <- read.csv("Overlap_heatpie_PtrVert24_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100



##################### PtrVert34
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrVert34_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrVert34_rep_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrVert34_rep',
    plot_name = 'PtrVert34_rep',
    sp1 = "PtrVertical03_",
    df_path = paste0('Overlap_heatpie_', "PtrVert34_rep", '_on_Ptr.csv')
    )

#

#a = readRDS("integration_PtrVert2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds")
#b = readRDS("RDS_PtrVert2_subtreeCenter.rds")


ovelaps <- read.csv("Overlap_heatpie_PtrVert34_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100

percent

##################### PtrTens12_rep
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens12_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens12_rep_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens12_rep',
    plot_name = 'PtrTens12_rep',
    sp1 = "TW_",
    df_path = paste0('Overlap_heatpie_', "PtrTens12_rep", '_on_Ptr.csv')

    )

ovelaps <- read.csv("Overlap_heatpie_PtrTens12_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrTens12_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)

##################### PtrTens13_rep


ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens13_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens13_rep_subtreeCenter.rds')
if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens13_rep',
    plot_name = 'PtrTens13_rep',
    sp1 = "TW_",
    df_path = paste0('Overlap_heatpie_', "PtrTens13_rep", '_on_Ptr.csv')

    )


ovelaps <- read.csv("Overlap_heatpie_PtrTens13_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100

input <- "Overlap_heatpie_PtrTens13_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)

##################### PtrTens14_rep

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens14_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens14_rep_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens14_rep',
    plot_name = 'PtrTens14_rep',
    sp1 = "TW_",
    df_path = paste0('Overlap_heatpie_', "PtrTens14_rep", '_on_Ptr.csv')

    )

ovelaps <- read.csv("Overlap_heatpie_PtrTens14_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100

input <- "Overlap_heatpie_PtrTens14_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)


##################### PtrTens23_rep
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens23_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens23_rep_subtreeCenter.rds')
if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens23_rep',
    plot_name = 'PtrTens23_rep',
    sp1 = "Tension02_",
    df_path = paste0('Overlap_heatpie_', "PtrTens23_rep", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrTens23_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100

input <- "Overlap_heatpie_PtrTens23_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)


##################### PtrTens24_rep


ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens24_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens24_rep_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens24_rep',
    plot_name = 'PtrTens24_rep',
    sp1 = "Tension02_",
    df_path = paste0('Overlap_heatpie_', "PtrTens24_rep", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrTens24_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100



input <- "Overlap_heatpie_PtrTens24_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)



##################### PtrTens34_rep

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens34_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens34_rep_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens34_rep',
    plot_name = 'PtrTens34_rep',
    sp1 = "Tension03_",
    df_path = paste0('Overlap_heatpie_', "PtrTens34_rep", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrTens34_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100



input <- "Overlap_heatpie_PtrTens34_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)


##################### PtrOp12_rep
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp12_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp12_rep_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp12_rep',
    plot_name = 'PtrOp12_rep',
    sp1 = "OW_",
    df_path = paste0('Overlap_heatpie_', "PtrOp12_rep", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp12_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100



input <- "Overlap_heatpie_PtrOp12_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)

##################### PtrOp13_rep


ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp13_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp13_rep_subtreeCenter.rds')
if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}

centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp13_rep',
    plot_name = 'PtrOp13_rep',
    sp1 = "OW_",
    df_path = paste0('Overlap_heatpie_', "PtrOp13_rep", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp13_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrOp13_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)

##################### PtrOp14_rep

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp14_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp14_rep_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp14_rep',
    plot_name = 'PtrOp14_rep',
    sp1 = "OW_",
    df_path = paste0('Overlap_heatpie_', "PtrOp14_rep", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp14_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100



input <- "Overlap_heatpie_PtrOp14_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)


##################### PtrOp23_rep
ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp23_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp23_rep_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp23_rep',
    plot_name = 'PtrOp23_rep',
    sp1 = "Opposite02_",
    df_path = paste0('Overlap_heatpie_', "PtrOp23_rep", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp23_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100

input <- "Overlap_heatpie_PtrOp23_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)


##################### PtrOp24_rep


ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp24_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp24_rep_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp24_rep',
    plot_name = 'PtrOp24_rep',
    sp1 = "Opposite02_",
    df_path = paste0('Overlap_heatpie_', "PtrOp24_rep", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp24_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrOp24_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)

##################### PtrOp34_rep

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp34_rep_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp34_rep_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp34_rep',
    plot_name = 'PtrOp34_rep',
    sp1 = "Opposite03_",
    df_path = paste0('Overlap_heatpie_', "PtrOp34_rep", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp34_rep_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrOp34_rep_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)


##################### PtrTens1

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens1_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens1_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens1',
    plot_name = 'PtrTens1',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrTens1", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrTens1_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrTens1_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)

##################### PtrTens2

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens2_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens2',
    plot_name = 'PtrTens2',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrTens2", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrTens2_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrTens2_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)


##################### PtrTens3

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens3_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens3_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens3',
    plot_name = 'PtrTens3',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrTens3", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrTens3_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrTens3_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)

##################### PtrTens4

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrTens4_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrTens4_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrTens4',
    plot_name = 'PtrTens4',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrTens4", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrTens4_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrTens4_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)


##################### PtrOp1

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp1_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp1_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp1',
    plot_name = 'PtrOp1',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrOp1", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp1_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrOp1_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)


##################### PtrOp2

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp2_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp2',
    plot_name = 'PtrOp2',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrOp2", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp2_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrOp2_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)


##################### PtrOp3

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp3_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp3_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp3',
    plot_name = 'PtrOp3',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrOp3", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp3_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrOp3_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)

##################### PtrOp4

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'integration_PtrOp4_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOp4_subtreeCenter.rds')

if (is.null(ClaPAIR_subtreeCenter)) {
  stop("ClaPAIR_subtreeCenter not found. Stopping execution.")
}
centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'PtrOp4',
    plot_name = 'PtrOp4',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "PtrOp4", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_PtrOp4_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_PtrOp4_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)




################## Ath

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Ath2019R.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrAth2019R_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'Ptr_Ath2019R',
    plot_name = 'Ptr_Ath2019R',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "Ptr_Ath2019R", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_Ptr_Ath2019R_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_Ptr_Ath2019R_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)




################## AthL (todo) map Rstudio

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Ath2021L.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrAth2021L_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'Ptr_Ath2021L',
    plot_name = 'Ptr_Ath2021L',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "Ptr_Ath2021L", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_Ptr_Ath2021L_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_Ptr_Ath2021L_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)

################## rice (todo) map Rstudio

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_SF10000_AF2000_KA5_Ptr_Osa2021R.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_PtrOsa2021R_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'Ptr_Osa2021R',
    plot_name = 'Ptr_Osa2021R',
    sp1 = "Ptr_",
    df_path = paste0('Overlap_heatpie_', "Ptr_Osa2021R", '_on_Ptr.csv')

    )
ovelaps <- read.csv("Overlap_heatpie_Ptr_Osa2021R_on_Ptr.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100


input <- "Overlap_heatpie_Ptr_Osa2021R_on_Ptr.csv"
write.table(input, "percent_values.txt", col.names = FALSE, append = TRUE)
write.table(percent, "percent_values.txt", col.names = FALSE, append = TRUE)







