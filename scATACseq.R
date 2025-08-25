
###Note：All the directories in this chunk of codes are my local directory
###You need to customize your directory to run the whole process
##Prepare all the packages for the following running
library(devtools)
library(ArchR)
library(pheatmap)
library(tidyr)
library(Rsamtools)
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(ggplot2)
library(stringr)
library(SingleR)
library(scran)
library(dittoSeq)
library(mclust)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Human genome annotation bank
library(org.Hs.eg.db)
library(igraph)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(UpSetR)
set.seed(44)
addArchRGenome("hg38")
addArchRThreads(threads = 4)

###First part started
#Download the 6 patients from GSE173682: Patient 6 to Patient 11

input.file.list = c('~/data_my_scATAC/Data/GSM5276949_37EACL_ATAC_fragments.tsv.gz',
                    '~/data_my_scATAC/Data/GSM5276950_38FE7L_ATAC_fragments.tsv.gz',
                    '~/data_my_scATAC/Data/GSM5276951_3BAE2L_ATAC_fragments.tsv.gz',
                    '~/data_my_scATAC/Data/GSM5276952_3E5CFL_ATAC_fragments.tsv.gz',
                    '~/data_my_scATAC/Data/GSM5276953_3CCF1L_ATAC_fragments.tsv.gz',
                    '~/data_my_scATAC/Data/GSM5276954_3E4D1L_ATAC_fragments.tsv.gz')

sampleNames = c("P1","P2","P3","P4","P5","P6")


#Create ArrowFiles 
ArrowFiles = createArrowFiles(inputFiles = input.file.list,
                              sampleNames = sampleNames,
                              minTSS = 0,
                              minFrags = 0,
                              addTileMat = TRUE,
                              addGeneScoreMat = TRUE,
                              excludeChr = c("chrM", "chrY", "chrX"))

#Create ArchRProject
projncov = ArchRProject(ArrowFiles = ArrowFiles,
                        outputDirectory = "ATAC_out", # output directory
                        copyArrows = TRUE) # save the backup for ArrowFiles

#Determine the doublets

doubScores <- addDoubletScores(input = ArrowFiles,
                               k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                               knnMethod = "UMAP",
                               useMatrix = "TileMatrix",
                               nTrials=5,
                               LSIMethod = 1,
                               scaleDims = F,
                               corCutOff = 0.75,
                               UMAPParams = list(n_neighbors =30,
                                                 min_dist = 0.3,
                                                 metric = "cosine",
                                                 verbose =T),dimsToUse = 1:50)


# Filter out outlier low quality cells and doublets
# GMM for fragments per cell

for (i in sampleNames) {
  projncov.i <- projncov[projncov$Sample == i, ]
  projncov.i <- projncov.i[!is.na(projncov.i$nFrags) & projncov.i$nFrags > 0, ]
  projncov.i <- projncov.i[!is.na(projncov.i$TSSEnrichment) & projncov.i$TSSEnrichment >= 0, ]
  
  # GMM for fragments per cell
  depth.clust <- Mclust(log10(projncov.i$nFrags), G = 2)
  projncov.i$depth.cluster <- depth.clust$classification
  projncov.i$depth.cluster.uncertainty <- depth.clust$uncertainty
  
  plot1 <- ggPoint(
    x = log10(projncov.i$nFrags),
    y = log10(projncov.i$TSSEnrichment + 1),
    color = as.character(projncov.i$depth.cluster),
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n", i, " log10(fragments)"))
  
  ggsave(paste0(i, "_depth.pdf"), plot = plot1, width = 4, height = 4)
  
  # GMM for TSS per cell
  TSS.clust <- Mclust(log10(projncov.i$TSSEnrichment + 1), G = 2)
  projncov.i$TSS.cluster <- TSS.clust$classification
  projncov.i$TSS.cluster.uncertainty <- TSS.clust$uncertainty
  
  plot2 <- ggPoint(
    x = log10(projncov.i$nFrags),
    y = log10(projncov.i$TSSEnrichment + 1),
    color = as.character(projncov.i$TSS.cluster),
    discrete = TRUE,
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n", i, " TSS Enrichment"))
  
  ggsave(paste0(i, "_TSS.pdf"), plot = plot2, width = 4, height = 4)
  
  # Save TSS-clustered data
  df.TSS <- data.frame(projncov.i$cellNames, projncov.i$TSS.cluster, projncov.i$TSS.cluster.uncertainty, projncov.i$TSSEnrichment)
  colnames(df.TSS) <- c("cellNames", "TSS.cluster", "TSS.cluster.uncertainty", "TSSEnrichment")
  df.TSS <- dplyr::filter(df.TSS, TSS.cluster == 2, TSS.cluster.uncertainty <= 0.05)
  saveRDS(df.TSS, paste0("df_TSS_", i, ".rds"))
  
  # Save depth-clustered data
  df.depth <- data.frame(projncov.i$cellNames, projncov.i$depth.cluster, projncov.i$depth.cluster.uncertainty, projncov.i$nFrags)
  colnames(df.depth) <- c("cellNames", "depth.cluster", "depth.cluster.uncertainty", "nFrags")
  df.depth <- dplyr::filter(df.depth, depth.cluster == 2, depth.cluster.uncertainty <= 0.05)
  saveRDS(df.depth, paste0("df_depth_", i, ".rds"))
  
  # QC thresholds plot
  if (nrow(df.TSS) > 0 & nrow(df.depth) > 0) {
    plot3 <- ggPoint(
      x = log10(projncov.i$nFrags),
      y = log10(projncov.i$TSSEnrichment + 1),
      colorDensity = TRUE,
      continuousSet = "sambaNight",
      xlabel = "log10(unique fragments)",
      ylabel = "log10(TSS Enrichment+1)"
    ) +
      geom_hline(yintercept = log10(min(df.TSS$TSSEnrichment) + 1), linetype = "dashed") +
      geom_vline(xintercept = log10(min(df.depth$nFrags)), linetype = "dashed") +
      ggtitle(paste0("QC thresholds:\n", i))
    
    ggsave(paste0(i, "_QC.pdf"), plot = plot3, width = 4, height = 4)
  }
  
  # Doublet Enrichment plot
  if (nrow(df.TSS) > 0 & nrow(df.depth) > 0) {
    plot4 <- ggPoint(
      x = log10(projncov.i$nFrags),
      y = log10(projncov.i$TSSEnrichment + 1),
      color = projncov.i$DoubletEnrichment,
      discrete = FALSE,
      continuousSet = "sambaNight",
      xlabel = "log10(unique fragments)",
      ylabel = "log10(TSS Enrichment+1)"
    ) +
      geom_hline(yintercept = log10(min(df.TSS$TSSEnrichment) + 1), linetype = "dashed") +
      geom_vline(xintercept = log10(min(df.depth$nFrags)), linetype = "dashed") +
      ggtitle(paste0("Doublet Enrichment:\n", i))
    
    ggsave(paste0(i, "_doublets.pdf"), plot = plot4, width = 4, height = 4)
  }
}

# Set the QC files directory
qc_path <- "~/data_my_scATAC/Data/scATACseq_QC_files"

# Read depth QC controlled cells

list.depth <- list.files(path = qc_path, pattern = "^df_depth", full.names = TRUE)
df.depth.list <- lapply(list.depth, readRDS)
df.depth <- do.call(rbind, df.depth.list)
colnames(df.depth) <- c("cellNames", "depth.cluster", "depth.cluster.uncertainty", "nFrags")

# Read the TSS QC controlled cells
list.TSS <- list.files(path = qc_path, pattern = "^df_TSS", full.names = TRUE)
df.TSS.list <- lapply(list.TSS, readRDS)
df.TSS <- do.call(rbind, df.TSS.list)
colnames(df.TSS) <- c("cellNames", "TSS.cluster", "TSS.cluster.uncertainty", "TSSEnrichment")

# Take the intersect，select the cells meeting the standard for both QC
cellsPass <- intersect(df.TSS$cellNames, df.depth$cellNames)
cellsFail <- projncov$cellNames[!(projncov$cellNames %in% cellsPass)]

# Filter out the lowQC cells
projncov.filter <- projncov[projncov$cellNames %in% cellsPass, ]

# Remove the Doublets 
projncov <- filterDoublets(projncov.filter, filterRatio = 1, cutEnrich = 1, cutScore = -Inf)

# save the clean `projncov`
saveRDS(projncov, "~/data_my_scATAC/Data/projncov_filtered.rds")


# Fragment Size Histogram
plot1 <- plotFragmentSizes(projncov) + ggtitle("Fragment Size Histogram")
ggsave(file.path(qc_path, "Frags_hist.pdf"), plot = plot1, width = 6, height = 4)

# TSS Enrichment distribution
plot2 <- plotTSSEnrichment(projncov) + ggtitle("TSS Enrichment")
ggsave(file.path(qc_path, "TSS.pdf"), plot = plot2, width = 6, height = 4)


#####scATAC formal analysis started 
# Perform LSI reduction and clustering with ATAC data only
# Add LSI dimreduc

proj.filter <- addIterativeLSI(
  ArchRProj = projncov,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 4,
  LSIMethod = 2,
  scaleDims = T,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  UMAPParams = list(n_neighbors =30,
                    min_dist = 0.3,
                    metric = "cosine",
                    verbose =FALSE),
  varFeatures = 25000,
  dimsToUse = 1:50,
  binarize = T,
  corCutOff = 0.75,
  force = T,
  seed=44
)

#Addclusters 

proj.filter <- addClusters(
  input = proj.filter,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "ATAC_clusters",
  resolution = 0.7,
  dimsToUse = 1:50,force = T
)

#AddUMAP

proj.filter <- addUMAP(proj.filter,nNeighbors = 30,minDist = 0.3,dimsToUse = 1:50,metric = "cosine",force = T,reducedDims = "IterativeLSI")

# We keep the cluster with  >= 400 counts
valid_clusters <- names(cluster_counts[cluster_counts >= 400])
proj.filter <- proj.filter[proj.filter$ATAC_clusters %in% valid_clusters, ]

#Save the file regularly
proj.filter=readRDS("~/data_my_scATAC/Data/projncov_filtered.rds")

#Visualize the UMAP map 
p7 <- plotEmbedding(ArchRProj = proj.filter,
                    colorBy = "cellColData",
                    name = "ATAC_clusters",
                    embedding = "UMAP",
                    size = 0.5,plotWidth = 10, 
                    plotHeight = 8 )# 


p7

#Polish the UMAP map 

p7 <- p7 + theme(
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 14),
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold")  # Bold the title
) +
  labs(x = "UMAP_1", y = "UMAP_2", color = "Clusters") +  
  ggtitle("scATAC clusters")              


p7

#use getMarkerFeatures to get the marker genes for each cluster
markersGS <- getMarkerFeatures(ArchRProj = proj.filter,
                               useMatrix = "GeneScoreMatrix",
                               groupBy = "ATAC_clusters",
                               bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon",
                               maxCells = 1000)

# Use getMarkers to select the marker genes for each cluster，we use "FDR <= 0.05 & Log2FC >= 0.5" here

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")


# Use addImputeWeights to enhance signal visualization in downstream plots
proj.filter <- addImputeWeights(ArchRProj = proj.filter)

# set the NK common marker genes
NK_markers=c("NCAM1","KLRC1","KLRC2","B3GAT1","FCGR3A","CD69","KLRD1","NCR1","NCR2","NCR3","KLRK1","PTPRC")

# Visualize the NK common marker genes by heatmap 
heatmapGS <- plotMarkerHeatmap(seMarker = markersGS, #from getMarkerFeatures
                               cutOff = "FDR <= 0.01 & Log2FC >= 1", # Threshold
                               labelMarkers = NK_markers,# NK common marker genes 
                               binaryClusterRows = TRUE,
                               clusterCols = TRUE,
                               transpose = FALSE)

hm <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "right", annotation_legend_side = "bot")

hm

#Save the file regularly
saveRDS(proj.filter,"~/data_my_scATAC/Data/projncov_filtered_3.rds")



###Stop the scATAC processing temporarily and process the matched scRNA data from the same patients

###Run the process steps one by one including NK cell annotations and then merge them together, Script in another R file

###Integrate scRNA与scATAC data to annotate the NK cells in scATAC clusters

#Read the scATAC and merged scRNA rds files 

proj.filter=readRDS("~/data_my_scATAC/Data/projncov_filtered_3.rds")
rna=readRDS("~/data_my_scATAC/Single_RNA_data/rna_merged_processed.rds")

#This code constructs a groupList that maps matched RNA and ATAC cells for each sample, 
#enabling sample-aware integration in addGeneIntegrationMatrix()

groupList <- SimpleList()
for (i in levels(factor(proj.filter$Sample))){
  
  rna.sub <- rna[,rna$Simple_ID_new == i]
  RNA.cells <- colnames(rna.sub)
  
  idxSample <- BiocGenerics::which(proj.filter$Sample == i)
  cellsSample <- proj.filter$cellNames[idxSample]
  proj.filter_2 <- proj.filter[cellsSample, ]
  ATAC.cells <- proj.filter_2$cellNames
  
  groupList[[i]] <- SimpleList(
    ATAC = ATAC.cells,
    RNA = RNA.cells
  )
}

# Run addGeneIntegrationMatrix
set.seed(10)
proj.filter <- addGeneIntegrationMatrix(ArchRProj = proj.filter, #ATAC proj object
                                        useMatrix = "GeneScoreMatrix",
                                        matrixName = "GeneIntegrationMatrix",
                                        reducedDims = "IterativeLSI",
                                        seRNA = rna, #scRNA object
                                        groupList = groupList,
                                        addToArrow = TRUE,
                                        groupRNA = "SingleR.NK.merged",# NK Annotation in scRNA object 
                                        nameCell = "predictedCell",
                                        nameGroup = "predictedGroup",
                                        nameScore = "predictedScore",
                                        transferParams = list(dims = 1:50),
                                        force=TRUE)

#Visualize the integration result
plotEmbedding(proj.filter, colorBy = "cellColData", name = "predictedGroup",size = 0.05)

#Determine which ATAC cluster belongs to NK cluster according to the prediction result

cM <- confusionMatrix(proj.filter$ATAC_clusters, proj.filter$predictedGroup)
cM <- cM / Matrix::rowSums(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
mapLabs <- cbind(rownames(cM), labelNew)
pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), border_color = "black")

proj.filter$Clusters_sub <- mapLabs[match(proj.filter$ATAC_clusters, mapLabs[,1]),2]



#Extract the data，beautify the scATAC UMAP wiith the NK cell annotation

cell_types <- unique(proj.filter$predictedGroup)
df <- proj.filter@embeddings@listData$UMAP$df
df$Sample <- proj.filter$Sample[match(rownames(df), proj.filter$cellNames)]
df$Cluster <- proj.filter$ATAC_clusters[match(rownames(df), proj.filter$cellNames)]
df$Group <- proj.filter$predictedGroup[match(rownames(df), proj.filter$cellNames)]
colnames(df)[1:2] <- c("UMAP1", "UMAP2")

#Visualize the UMAP map with NK cells annotated
ggplot(df, aes(x = UMAP1, y = UMAP2, col = Group)) +
  geom_point(size = 0.2) +
  theme_classic(base_size = 14) +
  scale_color_manual(
    values = c("#fc8d62", "#a6cee3")
  ) +
  labs(
    x = "UMAP_1", y = "UMAP_2",
    color = "Group"
  ) +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3))
  ) +
  annotate(
    "text", x = 1, y = 6, label = "NK cells",
    size = 12,  color = "black"
  )

#save the scATAC object with the integration and NK annotation 

saveRDS(proj.filter,"~/data_my_scATAC/Data/proj_LSI_GeneScores_Annotations_Int.rds")

###First part finished





###Seond part started: scATAC NK data 

#read the scATAC object 
proj.filter <- readRDS("~/data_my_scATAC/Data/proj_LSI_GeneScores_Annotations_Int.rds")

#Extract the scATAC NK part 

proj.NK <- subsetArchRProject(
  ArchRProj = proj.filter,
  cells = proj.filter$cellNames[proj.filter$predictedGroup == "NK cells"],
  outputDirectory = "NK_cells_ATAC_proj", # 输出路径，可改
  dropCells = TRUE,force=T
)

#save the original NK ATAC object 

saveRDS(proj.NK,"~/data_my_scATAC/Data/proj.NK.rds")
proj.NK=readRDS("~/data_my_scATAC/Data/proj.NK.rds")

#Add IterativeLSI 
proj.NK.filter <- addIterativeLSI(
  ArchRProj = proj.NK,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 4,
  LSIMethod = 2,
  scaleDims = T,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 1000,
    n.start = 10
  ),
  UMAPParams = list(n_neighbors =30,
                    min_dist = 0.3,
                    metric = "cosine",
                    verbose =FALSE),
  varFeatures =  50000,
  dimsToUse = 1:50,
  binarize = T,
  corCutOff = 0.5,
  force = T,
  seed=40
)

#Add cluster
proj.NK.filter <- addClusters(
  input = proj.NK.filter,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "ATAC_NK_clusters",
  resolution = 1.5,
  dimsToUse = 1:20,force = T
)

#Add UMAP 
proj.NK.filter <- addUMAP(proj.NK.filter,nNeighbors = 10,minDist = 0.3,dimsToUse = 1:20,metric = "cosine",force = T,reducedDims = "IterativeLSI")


#Viusalize the scATAC NK UMAP 

p7 <- plotEmbedding(ArchRProj = proj.NK.filter,
                    colorBy = "cellColData",
                    name = "ATAC_NK_clusters",
                    embedding = "UMAP",
                    size = 1,plotWidth = 10,  # 
                    plotHeight = 8 )#  


p7
p7 <- p7 + theme(
  legend.text = element_text(size = 14),   # 
  legend.title = element_text(size = 14)   # 
)

p7


#Cleansing the NK ATAC clusters, remove the C1, which contains very few cells  

proj.NK.filter$ATAC_NK_clusters_original <- proj.NK.filter$ATAC_NK_clusters
cells_to_keep <- proj.NK.filter$cellNames[proj.NK.filter$ATAC_NK_clusters != "C1"]
proj.NK.filter <- subsetArchRProject(
  ArchRProj = proj.NK.filter,
  cells = cells_to_keep,
  outputDirectory = "NK_cells_ATAC_proj_no_C1",
  dropCells = TRUE,
  force = TRUE
)


clusters <- unique(proj.NK.filter$ATAC_NK_clusters)
#order the cluster (optional)
clusters <- sort(clusters)
#rename the clusters 
new_names <- paste0("C", seq_along(clusters))

# replace the old names with new names
proj.NK.filter$ATAC_NK_clusters <- plyr::mapvalues(
  proj.NK.filter$ATAC_NK_clusters,
  from = clusters,
  to = new_names
)

#re-visualize the UMAP without C1 

p7 <- plotEmbedding(ArchRProj = proj.NK.filter,
                    colorBy = "cellColData",
                    name = "ATAC_NK_clusters",
                    embedding = "UMAP",
                    size = 1,plotWidth = 10, labelSize = 5, # 
                    plotHeight = 8)# 


p7
p7 <- p7 + theme(
  legend.text = element_text(size = 14),  
  legend.title = element_text(size = 14) 
)
p7 <- p7 + theme(
  legend.text = element_text(size = 15),
  legend.title = element_text(size = 15),
  axis.title = element_text(size = 15),
  axis.text = element_text(size = 15),
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold") 
) +
  labs(x = "UMAP_1", y = "UMAP_2", color = "Clusters") +        
  ggtitle("NK ATAC clusters before QC")                #

p7



#addImputeWeights
proj.NK.filter <- addImputeWeights(ArchRProj =  proj.NK.filter)

#save the cleansed NK scATAC data
saveRDS(proj.NK.filter,"~/data_my_scATAC/Data/proj.NK_filtered.rds")
proj.NK.filter=readRDS("~/data_my_scATAC/Data/proj.NK_filtered.rds")


#Integration of NK scATAC with NK scRNA 
#Read the NK scRNA with the NK subsets annotations
rna_NK=readRDS("~/data_my_scATAC/Single_RNA_data/rna_NK.rds")
proj.NK.filter=readRDS("~/data_my_scATAC/Data/proj.NK_filtered.rds")
#addGeneIntegrationMatrix 
set.seed(10)
proj.NK.filter<- addGeneIntegrationMatrix(ArchRProj = proj.NK.filter, #ATAC proj
                                          useMatrix = "GeneScoreMatrix",
                                          matrixName = "GeneIntegrationMatrix",
                                          reducedDims = "IterativeLSI",
                                          seRNA = rna_NK, #NK scRNA 
                                          #groupList = groupList,
                                          addToArrow = TRUE,
                                          groupRNA = "NK_cell_type",#scRNA NK annotations
                                          nameCell = "predictedCell",
                                          nameGroup = "predictedGroup",
                                          nameScore = "predictedScore",
                                          transferParams = list(dims = 1:20),
                                          force=TRUE)




cM <- confusionMatrix(proj.NK.filter$ATAC_NK_clusters, proj.NK.filter$predictedGroup)
cM <- cM / Matrix::rowSums(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
mapLabs <- cbind(rownames(cM), labelNew)

#Visualize the distribution of the 
pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black",
  fontsize_row = 14,  
  fontsize_col = 14  
)
proj.NK.filter$Clusters_sub <- mapLabs[match(proj.NK.filter$ATAC_NK_clusters, mapLabs[,1]),2]

#Annotate the NK clusters. One note: Although C1 and C8 are annotated as adaptive NK cells, considering 
#C1 and C8 are separated in UMAP and C8 has highest annotation cell number, we defined the C8 as adaptive NK cells finally
#For C5, it doesn't have dominant annotated NK clusters, therefore, we defined it as" low quality" cluster

cluster_annotation <- data.frame(
  Cluster = paste0("C", 1:8),
  Clusters_sub = c(
    "low_q",                 # C1
    "NK_cluster_7",          # C2
    "NK_cluster_1",          # C3
    "Inflamed_NK_cells",     # C4
    "low_q",                 # C5
    "NK_cluster_0",          # C6
    "Conventional_NK_cells", # C7
    "Adaptive_NK_cells"      # C8
  )
)

#Apply the finalized annotation
proj.NK.filter$Clusters_sub <- cluster_annotation$Clusters_sub[
  match(proj.NK.filter$ATAC_NK_clusters, cluster_annotation$Cluster)
]
proj.NK.filter <- proj.NK.filter[proj.NK.filter$Clusters_sub != "low_q", ]

# Visualized the UMAP with the final annotation 
p1=plotEmbedding(proj.NK.filter, colorBy = "cellColData", name = "Clusters_sub",embedding = "UMAP",size = 2,labelSize = 5)
p1
p1 <- p1 + theme(
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 14),
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  plot.title = element_text(hjust = 0.5, size = 18, face = "bold")  
) +
  labs(x = "UMAP_1", y = "UMAP_2", color = "Clusters") +        
  ggtitle("NK ATAC clusters")              

p1
saveRDS(proj.NK.filter,"~/data_my_scATAC/Data/proj.NK.filter_LSI_GeneScores_Annotations_Int_final.rds")


########################Second step finished 

###Third step starated 
proj.NK.filter=readRDS("~/data_my_scATAC/Data/proj.NK.filter_LSI_GeneScores_Annotations_Int_final.rds")
#Create the Pseudo-bulk Replicates
proj.NK.filter <- addGroupCoverages(
  ArchRProj=proj.NK.filter,
  groupBy="Clusters_sub",
  minCells = 50,
  force=TRUE
)

#Directory for MACS2
pathToMacs2 <-"/disk2/user/yizhsu/biosoft/miniconda3/envs/ATAC/bin/macs2"

#Identify the peaks
proj.NK.filter <- addReproduciblePeakSet(
  ArchRProj = proj.NK.filter,
  groupBy = "Clusters_sub",
  peaksPerCell = 500,
  pathToMacs2 = pathToMacs2,
  force = TRUE)

# Add Peak Matrix
proj.NK.filter <- addPeakMatrix(ArchRProj = proj.NK.filter, force = TRUE)
saveRDS(proj.NK.filter,"~/data_my_scATAC/Data/proj.NK.filter__Annotations_final_with_peaks.rds")

proj.NK.filter=readRDS("~/data_my_scATAC/Data/proj.NK.filter__Annotations_final_with_peaks.rds")

#Identify Marker Peaks
all_ClusterCell <- unique(proj.NK.filter$Clusters_sub)
markerPeaks <- getMarkerFeatures(ArchRProj = proj.NK.filter,
                                 useMatrix = "PeakMatrix",
                                 groupBy = "Clusters_sub",
                                 useGroups = all_ClusterCell,
                                 bias = c("TSSEnrichment", "log10(nFrags)"),
                                 testMethod = "wilcoxon")

#Select the peaks 
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.25")

#Integration of peaks information for the following operation
marker_Peaks = markerList@listData
for(x in seq_along(marker_Peaks)){
  marker_Peaks[[x]] <- as.data.frame(marker_Peaks[[x]])
}
for(x in seq_along(marker_Peaks)){
  if(nrow(marker_Peaks[[x]])==0){
    marker_Peaks[[x]][1,] <- names(marker_Peaks[x])
  }else{
    marker_Peaks[[x]]$Clusters_cell <- names(marker_Peaks[x])
  }
}
#Exclude the empty data frame
valid_marker_Peaks <- marker_Peaks[sapply(marker_Peaks, nrow) > 1]
# Integrate and add the cluster source 
marker_Peaks_df <- do.call(rbind, lapply(names(valid_marker_Peaks), function(name) {
  df <- valid_marker_Peaks[[name]]
  df$Cluster <- name
  return(df)
}))

#Save the files for the following process 
saveRDS(marker_Peaks_df, "~/data_my_scATAC/Data/marker_Peaks_clean.rds")
write.csv(marker_Peaks_df, file = '~/data_my_scATAC/Data/marker_Peaks.csv')
marker_Peaks_df=readRDS("~/data_my_scATAC/Data/marker_Peaks_clean.rds")
marker_Peaks_df <- marker_Peaks_df %>%
  mutate(peak_id = paste0(seqnames, ":", start, "-", end))

#Heatmap for the DAR peaks
sig_peaks <- unique(marker_Peaks_df$peak_id)

peak_ids <- paste0(
  rowData(markerPeaks)$seqnames, ":",
  rowData(markerPeaks)$start, "-",
  rowData(markerPeaks)$end
)

mat <- assay(markerPeaks, "Mean") 
mat_sig <- mat[rownames(mat) %in% sig_peaks, ]
mat_z <- t(scale(t(mat_sig)))

clusters <- unique(marker_Peaks_df$Cluster)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#313695", "white", "red")
)

# heatmap visualization
Heatmap(
  mat_z,
  name = "Row Z-scores",
  col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "red")),  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 15),        
  column_title = "NK cell subtypes",
  #column_names_rot = 45, 
  column_title_gp = gpar(fontsize = 15),        
  row_title = "Differentially accessible peaks (DAPs)",
  row_title_gp = gpar(fontsize = 16),           
  border = TRUE,                                                    
  heatmap_legend_param = list(
    title = "Row Z-scores",
    title_gp = gpar(fontsize = 15),
    labels_gp = gpar(fontsize = 15),
    legend_height = unit(4, "cm")
  )
)

###cis elements distribution
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


#Annotate the peaks using Chipseek
peak_unique <- marker_Peaks_df %>%
  distinct(seqnames, start, end, Cluster)

peaks_gr <- GRanges(
  seqnames = peak_unique$seqnames,
  ranges = IRanges(start = peak_unique$start, end = peak_unique$end),
  Cluster = peak_unique$Cluster
)

peakAnno <- annotatePeak(
  peaks_gr, 
  tssRegion = c(-3000, 3000),  
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, 
  annoDb = "org.Hs.eg.db"
)

peak_annot_df <- as.data.frame(peakAnno)
peak_annot_df$peak <- paste0(
  peak_annot_df$seqnames, ":", peak_annot_df$start, "-", peak_annot_df$end
)


peak_annot_df <- peak_annot_df %>%
  mutate(cis_category = case_when(
    grepl("Promoter", annotation, ignore.case = TRUE)   ~ "Promoter (±3kb)",
    grepl("5' UTR", annotation, ignore.case = TRUE)     ~ "5' UTR",
    grepl("3' UTR", annotation, ignore.case = TRUE)     ~ "3' UTR",
    grepl("Exon", annotation, ignore.case = TRUE)       ~ "Exon",
    grepl("Intron", annotation, ignore.case = TRUE)     ~ "Intron",
    grepl("Downstream", annotation, ignore.case = TRUE) ~ "Downstream",
    grepl("Intergenic", annotation, ignore.case = TRUE) ~ "Distal intergenic",
    TRUE                                                ~ "Other"
  ))

pal <- c(
  "Promoter (±3kb)"    = "#6A8DCF",  #
  "Intron"             = "#B6A5CF",  #
  "Distal intergenic"  = "#F4B96B",  # 
  "Exon"               = "#7FBF7F",  # 
  "3' UTR"              = "#B38C7A",  # 
  "5' UTR"              = "#B0B0B0"   # 
)
lev <- names(pal)
threshold <- 0.5  # 

# 
anno_sum <- tibble::as_tibble(peak_annot_df) |>
  dplyr::transmute(cis_category = as.character(cis_category)) |>
  dplyr::filter(!is.na(cis_category)) |>
  dplyr::count(cis_category, name = "n") |>
  dplyr::mutate(
    percent = 100 * n / sum(n),
    label   = paste0(round(percent, 1), "%"),
    cis_category = factor(cis_category, levels = lev)
  ) |>
  dplyr::arrange(cis_category) |>
  dplyr::mutate(ypos = cumsum(percent) - percent/2) 

# 
ggplot(anno_sum, aes(x = 1, y = percent, fill = cis_category)) +
  geom_col(color = "white", width = 1) +

  geom_text(
    data = subset(anno_sum, percent > threshold),
    aes(label = label, group = cis_category),
    position = position_stack(vjust = 0.5),
    size = 7, show.legend = FALSE
  ) +
  geom_segment(
    data = subset(anno_sum, percent <= threshold),
    aes(x = 1.00, xend = 1.10, y = ypos, yend = ypos),
    inherit.aes = FALSE, color = "grey40", linewidth = 0.3
  ) +
  geom_text(
    data = subset(anno_sum, percent <= threshold),
    aes(x = 1.5, y = ypos, label = label),
    inherit.aes = FALSE, size = 7, hjust = 0
  ) +
  coord_polar(theta = "y", start = pi/2, direction = -1) + 
  scale_fill_manual(values = pal, breaks = lev) +
  xlim(0.5, 1.6) + 
  theme_void() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        plot.margin = margin(t = 2, r = 10, b = 2, l = 2),
        plot.title = element_text( hjust = 0.5, size = 20,
                                   margin = margin(b = 0))) +
  labs(title = "Distribution of cis-regulatory elements(CREs) in DAPs",
       fill = "Regions of genome")


write.csv(peak_annot_df, "~/data_my_scATAC/Data/annotated_peaks_with_clusters.csv", row.names = FALSE)

##

pal <- c(
  "Promoter (±3kb)"    = "#6A8DCF",
  "Intron"             = "#B6A5CF",
  "Distal intergenic"  = "#F4B96B",
  "Exon"               = "#7FBF7F",
  "3' UTR"             = "#B38C7A",
  "5' UTR"             = "#B0B0B0"
)
lev <- names(pal)


cluster_cis <- tibble::as_tibble(peak_annot_df) %>%
  dplyr::mutate(
    cis_category = factor(as.character(cis_category), levels = lev)
  ) %>%
  dplyr::count(Cluster, cis_category, name = "n") %>%
  dplyr::group_by(Cluster) %>%
  dplyr::mutate(
    percent = 100 * n / sum(n),
    label = paste0(round(percent, 1), "%")  
  ) %>%
  dplyr::ungroup()


ggplot(cluster_cis, aes(x = Cluster, y = percent, fill = cis_category)) +
  geom_col(color = "white") +
  geom_text(
    data = subset(cluster_cis, percent >= 0.3),  
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 6
  ) +
  scale_fill_manual(values = pal, breaks = lev) +
  labs(
    y = "Percentage of DAPs locations in genome",
    fill = "Regions of genome",
    title = "DAPs distribution per NK cell subtypes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.text.x = element_text(size = 18),        
    axis.title = element_text(size = 16), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=18),
    axis.text.y = element_text(size=17),
    plot.title = element_text(hjust = 0.5, size = 18
   )
  )






#Get the top13 peaks for each cluster
topPeaks_list <- lapply(names(markerList), function(cluster){
  df <- as.data.frame(markerList[[cluster]])
  if(nrow(df) > 0){
    df <- df[order(df$Log2FC, decreasing = TRUE), ]
    df_top <- head(df, 13)
    df_top$Cluster <- cluster
    return(df_top)
  } else {
    return(NULL)
  }
})

topPeaks_df <- do.call(rbind, topPeaks_list)
topPeaks <- unique(rownames(topPeaks_df))

#subset markerPeaks object
markerPeaks_subset <- markerPeaks[rownames(markerPeaks) %in% topPeaks, ]

#Visualize the differentially enriched(DE) peaks for each cluster by heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks_subset,
  cutOff = "FDR <=0.1 & Log2FC >= 0.25",  
  labelMarkers = NULL,          
  binaryClusterRows = T,
  clusterCols = T,
  transpose = F
)

#Visualize the DE peaks for adaptive NK cells by volcano map
MAplot= markerPlot(seMarker = markerPeaks, name = "Adaptive_NK_cells", cutOff = "FDR <= 0.1 & Log2FC >= 0.25", plotAs = "MA")
MAplot <- MAplot +
  labs(
    x = "Log2 Mean",            
    y = "Log2 Fold Change" ,color="Groups"    
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16 ) 
  ) +
  geom_point(size = 2) 

#Identify the motif information 

proj.NK.filter <- addMotifAnnotations(ArchRProj = proj.NK.filter, motifSet = "cisbp", name = "Motif", force = TRUE)


#Peak to gene analysis 
proj.NK.filter <- addPeak2GeneLinks(
  ArchRProj = proj.NK.filter,
  reducedDims = "IterativeLSI",  
  useMatrix = "GeneIntegrationMatrix", # or "GeneIntegrationMatrix" if you integrated RNA
  corCutOff = 0.4,              # Threshold for correlation,0.45 as default
  maxDist = 250000               # Maximum distance
)


proj.NK.filter=readRDS("~/data_my_scATAC/Data/proj.NK.filter_final_with_all_analysis.rds")
#Extract the correlation result from the peak-to-gene analysis 

p2g_df <- getPeak2GeneLinks(
  ArchRProj = proj.NK.filter,
  corCutOff = 0.4,
  returnLoops = FALSE
)


p2g_df_df <- as.data.frame(p2g_df)
peaks <- getPeakSet(proj.NK.filter)
genes <- rowData(getMatrixFromProject(proj.NK.filter, useMatrix = "GeneIntegrationMatrix"))

p2g_df_df$peak <- paste0(
  seqnames(peaks)[p2g_df_df$idxATAC], ":",
  start(peaks)[p2g_df_df$idxATAC], "-",
  end(peaks)[p2g_df_df$idxATAC]
)

p2g_df_df$gene <- genes$name[p2g_df_df$idxRNA]

#Read the annotated peak csv file
peak_annot_df=read.csv("~/data_my_scATAC/Data/annotated_peaks_with_clusters.csv")

#Extract the p2g data with the peaks in adaptive and conventional NK cells

p2g_filtered <- p2g_df_df[p2g_df_df$peak %in% peak_annot_df$peak, ]
p2g_filtered <- left_join(
  p2g_filtered,
  peak_annot_df[, c("peak", "Cluster", "SYMBOL")],
  by = "peak"
)
p2g_filtered <- p2g_filtered[order(p2g_filtered$Cluster),]

#save the p2g final csv file
write.csv(p2g_filtered,"~/data_my_scATAC/Data/sig_peak2gene_with_cluster.csv",row.names = F)


#Read the Motif-Matches-In-Peaks.rds files containing the matched TF motifs
motif_matches <- readRDS("/disk2/user/yizhsu/data_my_scATAC/Single_RNA_data/NK_cells_ATAC_proj_no_C1/Annotations/Motif-Matches-In-Peaks.rds")

#######Extract the matrix containing the information whether TF motifs apppear in certain peaks
motif_matrix <- assay(motif_matches, "matches") 

#Make the coordinates of location of peaks
peak_coords <- paste0(
  seqnames(rowRanges(motif_matches)), ":",
  start(rowRanges(motif_matches)), "-",
  end(rowRanges(motif_matches))
)

# Replace the rownames with the peak coordinates
rownames(motif_matrix) <- peak_coords
# Update the colnames
motif_names <- gsub("_[0-9]+$", "", colnames(motif_matrix))
colnames(motif_matrix) <- motif_names

# Transder to data.frame format
motif_df_full <- as.data.frame(as.matrix(motif_matrix))
motif_df_full$peak <- rownames(motif_matrix)

#Read the p2g csv file
p2g_filtered=read.csv("~/data_my_scATAC/Data/sig_peak2gene_with_cluster.csv")


motif_subset <- motif_df_full[motif_df_full$peak %in% p2g_filtered$peak, ]
motif_subset_annot <- merge(
  motif_subset,                    
  p2g_filtered[, c("peak", "Cluster")],  
  by = "peak"
)
motif_cols <- sapply(motif_subset_annot, is.logical)
motif_matrix_only <- motif_subset_annot[, motif_cols]

# Find the name of TF motif for each"TRUE" 
motif_list_per_peak <- apply(motif_matrix_only, 1, function(x) {
  tf_names <- names(x)[x == TRUE]
  return(tf_names)
})

names(motif_list_per_peak) <- motif_subset_annot$peak

# Transfer the list to data.frame
motif_peak_df <- data.frame(
  peak = names(motif_list_per_peak),
  matched_TFs = sapply(motif_list_per_peak, function(x) paste(x, collapse = ", "))
)

motif_peak_df <- merge(
  motif_peak_df,
  p2g_filtered[, c("peak", "Cluster","gene","SYMBOL")],
  by = "peak"
)

motif_peak_df <- motif_peak_df[order(motif_peak_df$Cluster), ]
motif_peak_df <- unique(motif_peak_df)

#Save the csv file
write.csv(motif_peak_df,"~/data_my_scATAC/Data/Motif_TF_gene_with_clusters.csv",row.names = FALSE)

motif_peak_df=read.csv("~/data_my_scATAC/Data/Motif_TF_gene_with_clusters.csv")
peak_annot_df=read.csv("~/data_my_scATAC/Data/annotated_peaks_with_clusters.csv")




#Choose the TOP50 TF motifs for adaptive NK cluster and conventional NK cluster

motif_peak_df=read.csv("~/data_my_scATAC/Data/Motif_TF_gene_with_clusters.csv")

#For adaptive NK cells
adaptive_df <- motif_peak_df[motif_peak_df$Cluster == "Adaptive_NK_cells", ]

adaptive_long <- adaptive_df %>%
  dplyr::select(peak, gene, matched_TFs) %>%
  separate_rows(matched_TFs, sep = ", ")

adaptive_long$matched_TFs <- trimws(adaptive_long$matched_TFs)

#Summarize the frequency(count) of the TFs in adaptive NK cells
tf2gene_adaptive <- adaptive_long %>%
  dplyr::filter(gene != "") %>%  
  dplyr::group_by(matched_TFs) %>%
  dplyr::summarise(
    genes = paste(unique(gene), collapse = ", "),
    gene_count = n_distinct(gene),
    .groups = "drop"
  ) %>%
  arrange(dplyr::desc(gene_count))



TF_aNK <- tf2gene_adaptive$matched_TFs

#For conventional NK cells
motif_peak_df=read.csv("~/data_my_scATAC/Data/Motif_TF_gene_with_clusters.csv")
conventional_df <- motif_peak_df[motif_peak_df$Cluster == "Conventional_NK_cells", ]

conventional_long <- conventional_df %>%
  dplyr::select(peak, gene, matched_TFs) %>%
  separate_rows(matched_TFs, sep = ", ")

conventional_long$matched_TFs <- trimws(conventional_long$matched_TFs)

#Summarize the frequency(count) of the TFs in conventional NK cells
tf2gene_conventional <- conventional_long %>%
  dplyr::filter(gene != "") %>%  
  dplyr::group_by(matched_TFs) %>%
  dplyr::summarise(
    genes = paste(unique(gene), collapse = ", "),
    gene_count = n_distinct(gene),
    .groups = "drop"
  ) %>%
  arrange(dplyr::desc(gene_count))

TF_cNK <- tf2gene_conventional$matched_TFs

################
##DAR heatmap with TFs

motif_peak_df=read.csv("~/data_my_scATAC/Data/Motif_TF_gene_with_clusters.csv")
peak_annot_df=read.csv("~/data_my_scATAC/Data/annotated_peaks_with_clusters.csv")

#motif_matches <- readRDS("/disk2/user/yizhsu/data_my_scATAC/Single_RNA_data/NK_cells_ATAC_proj_no_C1/Annotations/Motif-Matches-In-Peaks.rds")
#mm <- assay(motif_matches, "matches") 

#peak_coords <- paste0(seqnames(rowRanges(motif_matches)), ":", 
#                      start(rowRanges(motif_matches)), "-", 
#                      end(rowRanges(motif_matches)))
#rownames(mm) <- peak_coords
#colnames(mm) <- gsub("_[0-9]+$", "", colnames(mm))

#marker_Peaks_df=readRDS("~/data_my_scATAC/Data/marker_Peaks_clean.rds")
#marker_Peaks_df <- marker_Peaks_df %>%
#  mutate(peak_id = paste0(seqnames, ":", start, "-", end))
#sig_peaks <- unique(marker_Peaks_df$peak_id)
#peak_ids <- paste0(
#  rowData(markerPeaks)$seqnames, ":",
#  rowData(markerPeaks)$start, "-",
#  rowData(markerPeaks)$end
#)
#mat <- assay(markerPeaks, "Mean") 
#mat_sig <- mat[rownames(mat) %in% sig_peaks, ]
#mat_z <- t(scale(t(mat_sig)))

ht0  <- Heatmap(mat_z, cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE)
ht0d <- draw(ht0)
row_order_names <- rownames(mat_z)[unlist(row_order(ht0d))]   

topN <- 25
topTF <- head(tf2gene_adaptive$matched_TFs, topN) %>% 
  trimws() %>%
  union("STAT2")

tf_peak_long <- motif_peak_df %>%
  filter(Cluster == "Adaptive_NK_cells") %>%
  transmute(peak, Cluster, gene, SYMBOL, TFs = matched_TFs) %>%
  separate_rows(TFs, sep = ",\\s*") %>%
  mutate(TF = trimws(TFs)) %>%
  filter(TF %in% topTF) %>%
  arrange(TF) %>%               
  distinct(TF, .keep_all = TRUE) 

mat_ord <- mat_z[row_order_names, , drop = FALSE]
pos <- match(tf_peak_long$peak, rownames(mat_ord))
ok  <- !is.na(pos)
pos <- pos[ok]
labs <- tf_peak_long$TF[ok]

lab_by_pos <- tapply(labs, pos, function(v) paste(sort(unique(v)), collapse = "\n"))
pos_uniq  <- as.integer(names(lab_by_pos))
labs_uniq <- as.character(lab_by_pos)

left_ann <- rowAnnotation(
  `Top TFs (aNK)` = anno_mark(
    at = pos_uniq,
    labels = labs_uniq,
    side = "left",
    link_width = unit(2, "mm"),
    labels_gp = gpar(fontsize = 12, lineheight = 0.9)
  ),
  width = unit(6, "cm")
)

ht_main <- Heatmap(
  mat_ord,
  name = "Row Z-scores",
  col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "red")),
  cluster_rows = FALSE,              
  cluster_columns = TRUE,
  show_row_names = FALSE,            
  left_annotation = left_ann,
  border = TRUE
)
draw(ht_main)




#Visualize the adaptive NK TF counts and TF-gene regulatory map and 
#Top50 TF
#The count number for TOP TFs in adaptive NK cells

top_tf_adaptive <- motif_peak_df %>%
  filter(Cluster == "Adaptive_NK_cells") %>%
  pull(matched_TFs) %>%
  strsplit(", ") %>%
  unlist() %>%
  table() %>%
  sort(decreasing = TRUE)

#The count number for TOP TFs in conventional NK cells
top_tf_conventional <- motif_peak_df %>%
  filter(Cluster == "Conventional_NK_cells") %>%
  pull(matched_TFs) %>%
  strsplit(", ") %>%
  unlist() %>%
  table() %>%
  sort(decreasing = TRUE)

top50_tfs <- names(head(top_tf_adaptive, 50)) 
adaptive_top50 <- adaptive_long %>%
  filter(matched_TFs %in% top50_tfs)

tf_gene_edges_top50 <- adaptive_top50 %>%
  dplyr::select(from = matched_TFs, to = gene) %>%
  filter(!is.na(from), !is.na(to), from != "", to != "")

net <- tbl_graph(edges = tf_gene_edges_top50, directed = TRUE) %>%
  mutate(
    type   = ifelse(name %in% tf_gene_edges_top50$from, "TF", "Gene"),
    degree = centrality_degree(mode = "all")
  )

nodes <- as_tibble(net, active = "nodes")
order_tf   <- nodes %>% filter(type == "TF")   %>% arrange(desc(degree)) %>% mutate(y = row_number())
order_gene <- nodes %>% filter(type == "Gene") %>% arrange(desc(degree)) %>% mutate(y = row_number())
order_gene$y <- rescale(order_gene$y, to = range(order_tf$y))

nodes_pos <- bind_rows(order_tf, order_gene) %>%
  mutate(x = ifelse(type == "TF", 0, 0.30)) %>%    
  arrange(match(name, nodes$name))

node_names <- nodes$name
net <- net %>% 
  activate(edges) %>% mutate(to_gene = node_names[to]) %>% 
  activate(nodes)


lay <- create_layout(net, layout = "manual", x = nodes_pos$x, y = nodes_pos$y)

g <- ggraph(lay) +
  geom_edge_link(
    aes(color = to_gene), 
    arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
    end_cap = circle(2.2, "mm"),
    edge_width = 0.4, edge_alpha = 0.55, lineend = "round"
  ) +
  geom_node_point(aes(size = degree, color = type), alpha = 0.9) +
  geom_node_text(
    data = as.data.frame(lay) %>% dplyr::filter(type == "TF"),
    aes(x = x, y = y, label = name, color = type),
    hjust = 1, nudge_x = -0.025, vjust = 0.5, size = 6, show.legend = FALSE
  ) +
  geom_node_text(
    data = as.data.frame(lay) %>% dplyr::filter(type == "Gene"),
    aes(x = x, y = y, label = name, color = type),
    hjust = 0, nudge_x =  0.025, vjust = 0.5, size = 6, show.legend = FALSE
  ) +
  scale_color_manual(values = c(TF = "#005AB5", Gene = "#E66100")) +
  scale_edge_colour_discrete(guide = "none") +
  scale_size(range = c(3, 8)) +
  coord_cartesian(xlim = c(-0.35, 1.05), clip = "off") +
  theme_graph(base_family = "Arial") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 30, 10, 30)) +
  theme(
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    plot.title   = element_text(size = 16, hjust = 0.5)
  )+
  labs(
       color = "Node Type", size = "Degree", edge_color = "Target Gene")

print(g)





#Visualize the counts for TOP50 TFs in adaptive NK cells 

top_tf_adaptive_df <- as.data.frame(top_tf_adaptive)
colnames(top_tf_adaptive_df) <- c("TF", "count")

ggplot(head(top_tf_adaptive_df, 50), aes(x = reorder(TF, count), y = count, fill = count)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = count),
    hjust = -0.2, size = 5, color = "black"
  ) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top 50 TF Motifs Enriched in Adaptive NK Cells",
    x = "Enriched TF Motif",
    y = "Number of Peaks with Motif",
    fill = "Motif Count"
  ) +
  theme(
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray80", linetype = "dotted"),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradientn(
    colours = c("#3B528B", "#21908C", "#5DC863", "#FDE725")  # 紫→蓝→绿→黄
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    limits = c(0, max(top_tf_adaptive_df$count[1:50]) * 1.15)
  )

##GO enrichment for TOP50 TFs 
top50_TF <- head(top_tf_adaptive_df$TF, 50)

gene_df <- bitr(top50_TF, 
                fromType = "SYMBOL",
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene         = gene_df$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",       #or "MF", "CC", "ALL"
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)      
ego_df <- as.data.frame(ego)
ego_df <- ego_df %>%
  arrange(desc(Count)) %>%
  head(30)  # 

ggplot(ego_df, aes(x = reorder(Description, Count), 
                   y = Count, 
                   fill = p.adjust)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), 
            hjust = -0.2, 
            size = 4, 
            color = "black") +
  coord_flip() +
  scale_fill_gradient(low = "#d73027", high = "#fee090", name = "p.adjust") +
  labs(
    title = "GO Enrichment of Top 50 TFs in Adaptive NK Cells",
    x = "GO Term",
    y = "Gene Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

dotplot(ego, showCategory = 20) + 
  ggtitle("GO Enrichment of Top 50 TFs in Adaptive NK Cells") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),  
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16),  
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16), 
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 16)   
  )


ekegg <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = "hsa",       
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH"
)

ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
ekegg_df <- as.data.frame(ekegg@result) %>%
  arrange(desc(Count)) %>%
  head(30)
ggplot(ekegg_df, aes(x = reorder(Description, Count), 
                     y = Count, 
                     fill = p.adjust)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), 
            hjust = -0.2, 
            size = 4, 
            color = "black") +
  coord_flip() +
  scale_fill_gradient(low = "#2171b5", high = "#c6dbef",
                      name = "adj. p-value") +
  labs(
    title = "KEGG Enrichment of Top 50 TFs in Adaptive NK Cells",
    x = "KEGG Pathway",
    y = "Gene Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))







#Narrow down the TFs by overalapping with the differentially enriched TF regulons from PySCENIC 

IPLA_NK_rss <- read.csv("~/IPLA_tumor_NK_Pyscenic/IPLA_NK_rss_clusters_SD.csv",row.names = 1)
TF_from_pySCENIC <- colnames(IPLA_NK_rss)

#Take the TOP50 TFs for aNK and cNK 
TF_aNK <- unique(na.omit(tf2gene_adaptive$matched_TFs[1:50]))
TF_cNK <- unique(na.omit(tf2gene_conventional$matched_TFs[1:50]))

# Remove "...",only retain the gene names for TF from ATAC data
TF_from_pySCENIC <- sub("\\.\\.\\..*$", "", colnames(IPLA_NK_rss))

#Get the intersected genes between aNK cells in ATAC data and pySCENIC DE TFs, cNK cells in ATAC data and pySCENIC DE TFs   
shared_aNK_PySCENIC <- intersect(TF_aNK, TF_from_pySCENIC)
shared_cNK_PySCENIC <- intersect(TF_cNK, TF_from_pySCENIC)

#Take the union between shared_aNK_PySCENIC and shared_cNK_PySCENIC 
shared_all_union <- Reduce(union, list(shared_aNK_PySCENIC, shared_cNK_PySCENIC))

#Use Upset plot to show the overlapped TFs in adaptive NK cluster; conventional NK cluster and PySCENIC adaptive NK cells 
Upset_list_TF = list(pySCENIC_DA_TFs=TF_from_pySCENIC,
                    scATAC_aNK_enriched_TFs_TOP50=TF_aNK,
                    scATAC_cNK_enriched_TFs_TOP50=TF_cNK)

upset(fromList(Upset_list_TF), 
      sets = c("pySCENIC_DA_TFs", "scATAC_aNK_enriched_TFs_TOP50", "scATAC_cNK_enriched_TFs_TOP50"),
      order.by = "freq", 
      sets.bar.color = c("#FF8C00", "#4DAF4A", "#B64E89"),  
      text.scale = c(3, 3,2.5,2.5,2.5,3)) 


#Similar process but only visualize the 8 TFs in 21 TFs overlapped in aNK in ATAC data and pySCENIC  

adaptive_filtered <- adaptive_long %>%
  filter(matched_TFs %in% shared_aNK_PySCENIC)

tf_gene_edges_filtered <- adaptive_filtered %>%
  select(from = matched_TFs, to = gene) %>%
  filter(!is.na(from), !is.na(to), from != "", to != "")

net <- tbl_graph(edges = tf_gene_edges_filtered, directed = TRUE)
net <- net %>%
  mutate(
    type = ifelse(name %in% tf_gene_edges_filtered$from, "TF", "Gene"),
    degree = centrality_degree(mode = "all")
  )
set.seed(123)
layout <- create_layout(net, layout = "fr")

# 
ggraph(layout) +
  geom_edge_link(
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    end_cap = circle(3, "mm"),
    edge_width = 0.6,
    edge_alpha = 0.4,
    edge_colour = "gray50"
  ) +
  geom_node_point(
    aes(color = type, size = degree),
    alpha = 0.9
  ) +
  geom_node_text(
    aes(label = name, color = type),
    repel = TRUE,
    force = 1.2,                           
    box.padding = unit(0.6, "lines"),      
    size = 6,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("TF" = "#3182bd", "Gene" = "#e6550d")) +
  scale_size(range = c(3, 8)) +
  theme_graph(base_family = "Arial") +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "TF → Gene Regulatory Network in Adaptive NK Cells (Overlapped)",
    color = "Node Type",
    size = "Degree"
  )

#Choose the target genes which are CRCP and MTFP1 to draw the peak-to-gene correlation browser track map

plotPDF(
  plotBrowserTrack(
    ArchRProj = proj.NK.filter,
    groupBy = "Clusters_sub",
    geneSymbol = "CRCP",
    loops = getPeak2GeneLinks(proj.NK.filter)
  ),
  name = "CRCP_peak2gene_browserTrack.pdf",
  ArchRProj = proj.NK.filter,
  addDOC = FALSE,      
  width = 6,
  height = 5
)

plotPDF(
  plotBrowserTrack(
    ArchRProj = proj.NK.filter,
    groupBy = "Clusters_sub",
    geneSymbol = "MTFP1",
    loops = getPeak2GeneLinks(proj.NK.filter)
  ),
  name = "MTFP1_peak2gene_browserTrack.pdf",
  ArchRProj = proj.NK.filter,
  addDOC = FALSE,      
  width = 6,
  height = 5
)


