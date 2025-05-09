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
library(readxl)
library(ggpubr)
set.seed(44)
addArchRGenome("hg38")
addArchRThreads(threads = 4)



################Start to process the scRNA data matched with scATAC data#############

###Download the scRNA data from GSE173682
#P6 patient
counts <- Read10X(data.dir = "~/data_my_scATAC/Single_RNA_data/GSM5276938_P6")
rna <- CreateSeuratObject(counts = counts, min.cells = 3)

# Mito percentage 
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")

# Quality control
rna@meta.data$nCount_RNA_outlier_2mad <- isOutlier(log(rna$nCount_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$nFeature_RNA_outlier_2mad <- isOutlier(log(rna$nFeature_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$percent_mt_outlier_2mad <- isOutlier(log1p(rna$percent.mt), log = FALSE, type = "higher", nmads = 2)

rna <- subset(rna, subset = nCount_RNA_outlier_2mad == FALSE &
                nFeature_RNA_outlier_2mad == FALSE &
                percent_mt_outlier_2mad == FALSE)

#Standard workflow
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna <- ScaleData(rna)
rna <- RunPCA(rna)
rna <- FindNeighbors(rna, dims = 1:20)
rna <- FindClusters(rna, resolution = 0.7)
rna <- RunUMAP(rna, dims = 1:20)

# Cell annotations using SingleR
rna.sce <- as.SingleCellExperiment(rna)

#HPCA
ref.data.HPCA <- celldex::HumanPrimaryCellAtlasData()
saveRDS(ref.data.HPCA, file = "~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
#Blueprint
ref.data.BED <- celldex::BlueprintEncodeData()
saveRDS(ref.data.BED, file = "~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")

ref.data.HPCA <- readRDS("~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
ref.data.BED <- readRDS("~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")



predictions.HPCA <- SingleR(test = rna.sce, ref = ref.data.HPCA,
                            labels = ref.data.HPCA$label.main, assay.type.test = "logcounts")
predictions.BED <- SingleR(test = rna.sce, ref = ref.data.BED,
                           labels = ref.data.BED$label.main, assay.type.test = "logcounts")

# Save the annotation result
rna$SingleR.HPCA <- predictions.HPCA$pruned.labels
rna$SingleR.BED <- predictions.BED$pruned.labels

table(rna$SingleR.HPCA)
table(rna$SingleR.BED)

Idents(rna)="SingleR.BED"
DimPlot(rna,reduction = "umap",label = T)

rna$SingleR.NK.merged <- ifelse(
  grepl("NK", rna$SingleR.HPCA, ignore.case = TRUE) | 
    grepl("NK", rna$SingleR.BED, ignore.case = TRUE),
  "NK cells",
  "Other"
)
#Save the P6
saveRDS(rna, "~/data_my_scATAC/Single_RNA_data/processed_rds/p6_processed.rds")




#P7 patient
counts <- Read10X(data.dir = "~/data_my_scATAC/Single_RNA_data/GSM5276939_P7")
rna <- CreateSeuratObject(counts = counts, min.cells = 3)

# Mito percentage 
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")

# Quality control
rna@meta.data$nCount_RNA_outlier_2mad <- isOutlier(log(rna$nCount_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$nFeature_RNA_outlier_2mad <- isOutlier(log(rna$nFeature_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$percent_mt_outlier_2mad <- isOutlier(log1p(rna$percent.mt), log = FALSE, type = "higher", nmads = 2)

rna <- subset(rna, subset = nCount_RNA_outlier_2mad == FALSE &
                nFeature_RNA_outlier_2mad == FALSE &
                percent_mt_outlier_2mad == FALSE)

#Standard workflow
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna <- ScaleData(rna)
rna <- RunPCA(rna)
rna <- FindNeighbors(rna, dims = 1:20)
rna <- FindClusters(rna, resolution = 0.7)
rna <- RunUMAP(rna, dims = 1:20)

# Cell annotations using SingleR
rna.sce <- as.SingleCellExperiment(rna)

#HPCA
ref.data.HPCA <- celldex::HumanPrimaryCellAtlasData()
saveRDS(ref.data.HPCA, file = "~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
#Blueprint
ref.data.BED <- celldex::BlueprintEncodeData()
saveRDS(ref.data.BED, file = "~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")

ref.data.HPCA <- readRDS("~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
ref.data.BED <- readRDS("~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")



predictions.HPCA <- SingleR(test = rna.sce, ref = ref.data.HPCA,
                            labels = ref.data.HPCA$label.main, assay.type.test = "logcounts")
predictions.BED <- SingleR(test = rna.sce, ref = ref.data.BED,
                           labels = ref.data.BED$label.main, assay.type.test = "logcounts")

# Save the annotation result
rna$SingleR.HPCA <- predictions.HPCA$pruned.labels
rna$SingleR.BED <- predictions.BED$pruned.labels

table(rna$SingleR.HPCA)
table(rna$SingleR.BED)

Idents(rna)="SingleR.BED"
DimPlot(rna,reduction = "umap",label = T)

rna$SingleR.NK.merged <- ifelse(
  grepl("NK", rna$SingleR.HPCA, ignore.case = TRUE) | 
    grepl("NK", rna$SingleR.BED, ignore.case = TRUE),
  "NK cells",
  "Other"
)
#Save the P7
saveRDS(rna, "~/data_my_scATAC/Single_RNA_data/processed_rds/p7_processed.rds")




#P8 patient
counts <- Read10X(data.dir = "~/data_my_scATAC/Single_RNA_data/GSM5276940_P8")
rna <- CreateSeuratObject(counts = counts, min.cells = 3)

# Mito percentage 
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")

# Quality control
rna@meta.data$nCount_RNA_outlier_2mad <- isOutlier(log(rna$nCount_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$nFeature_RNA_outlier_2mad <- isOutlier(log(rna$nFeature_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$percent_mt_outlier_2mad <- isOutlier(log1p(rna$percent.mt), log = FALSE, type = "higher", nmads = 2)

rna <- subset(rna, subset = nCount_RNA_outlier_2mad == FALSE &
                nFeature_RNA_outlier_2mad == FALSE &
                percent_mt_outlier_2mad == FALSE)

#Standard workflow
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna <- ScaleData(rna)
rna <- RunPCA(rna)
rna <- FindNeighbors(rna, dims = 1:20)
rna <- FindClusters(rna, resolution = 0.7)
rna <- RunUMAP(rna, dims = 1:20)

# Cell annotations using SingleR
rna.sce <- as.SingleCellExperiment(rna)

#HPCA
ref.data.HPCA <- celldex::HumanPrimaryCellAtlasData()
saveRDS(ref.data.HPCA, file = "~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
#Blueprint
ref.data.BED <- celldex::BlueprintEncodeData()
saveRDS(ref.data.BED, file = "~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")

ref.data.HPCA <- readRDS("~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
ref.data.BED <- readRDS("~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")



predictions.HPCA <- SingleR(test = rna.sce, ref = ref.data.HPCA,
                            labels = ref.data.HPCA$label.main, assay.type.test = "logcounts")
predictions.BED <- SingleR(test = rna.sce, ref = ref.data.BED,
                           labels = ref.data.BED$label.main, assay.type.test = "logcounts")

# Save the annotation result
rna$SingleR.HPCA <- predictions.HPCA$pruned.labels
rna$SingleR.BED <- predictions.BED$pruned.labels

table(rna$SingleR.HPCA)
table(rna$SingleR.BED)

Idents(rna)="SingleR.BED"
DimPlot(rna,reduction = "umap",label = T)

rna$SingleR.NK.merged <- ifelse(
  grepl("NK", rna$SingleR.HPCA, ignore.case = TRUE) | 
    grepl("NK", rna$SingleR.BED, ignore.case = TRUE),
  "NK cells",
  "Other"
)
#Save the P8
saveRDS(rna, "~/data_my_scATAC/Single_RNA_data/processed_rds/p8_processed.rds")





#P9 patient
counts <- Read10X(data.dir = "~/data_my_scATAC/Single_RNA_data/GSM5276941_P9")
rna <- CreateSeuratObject(counts = counts, min.cells = 3)

# Mito percentage 
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")

# Quality control
rna@meta.data$nCount_RNA_outlier_2mad <- isOutlier(log(rna$nCount_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$nFeature_RNA_outlier_2mad <- isOutlier(log(rna$nFeature_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$percent_mt_outlier_2mad <- isOutlier(log1p(rna$percent.mt), log = FALSE, type = "higher", nmads = 2)

rna <- subset(rna, subset = nCount_RNA_outlier_2mad == FALSE &
                nFeature_RNA_outlier_2mad == FALSE &
                percent_mt_outlier_2mad == FALSE)

#Standard workflow
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna <- ScaleData(rna)
rna <- RunPCA(rna)
rna <- FindNeighbors(rna, dims = 1:20)
rna <- FindClusters(rna, resolution = 0.7)
rna <- RunUMAP(rna, dims = 1:20)

# Cell annotations using SingleR
rna.sce <- as.SingleCellExperiment(rna)

#HPCA
ref.data.HPCA <- celldex::HumanPrimaryCellAtlasData()
saveRDS(ref.data.HPCA, file = "~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
#Blueprint
ref.data.BED <- celldex::BlueprintEncodeData()
saveRDS(ref.data.BED, file = "~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")

ref.data.HPCA <- readRDS("~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
ref.data.BED <- readRDS("~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")



predictions.HPCA <- SingleR(test = rna.sce, ref = ref.data.HPCA,
                            labels = ref.data.HPCA$label.main, assay.type.test = "logcounts")
predictions.BED <- SingleR(test = rna.sce, ref = ref.data.BED,
                           labels = ref.data.BED$label.main, assay.type.test = "logcounts")

# Save the annotation result
rna$SingleR.HPCA <- predictions.HPCA$pruned.labels
rna$SingleR.BED <- predictions.BED$pruned.labels

table(rna$SingleR.HPCA)
table(rna$SingleR.BED)

Idents(rna)="SingleR.BED"
DimPlot(rna,reduction = "umap",label = T)

rna$SingleR.NK.merged <- ifelse(
  grepl("NK", rna$SingleR.HPCA, ignore.case = TRUE) | 
    grepl("NK", rna$SingleR.BED, ignore.case = TRUE),
  "NK cells",
  "Other"
)
#Save the P9
saveRDS(rna, "~/data_my_scATAC/Single_RNA_data/processed_rds/p9_processed.rds")




#P10 patient
counts <- Read10X(data.dir = "~/data_my_scATAC/Single_RNA_data/GSM5276942_P10")
rna <- CreateSeuratObject(counts = counts, min.cells = 3)

# Mito percentage 
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")

# Quality control
rna@meta.data$nCount_RNA_outlier_2mad <- isOutlier(log(rna$nCount_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$nFeature_RNA_outlier_2mad <- isOutlier(log(rna$nFeature_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$percent_mt_outlier_2mad <- isOutlier(log1p(rna$percent.mt), log = FALSE, type = "higher", nmads = 2)

rna <- subset(rna, subset = nCount_RNA_outlier_2mad == FALSE &
                nFeature_RNA_outlier_2mad == FALSE &
                percent_mt_outlier_2mad == FALSE)

#Standard workflow
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna <- ScaleData(rna)
rna <- RunPCA(rna)
rna <- FindNeighbors(rna, dims = 1:20)
rna <- FindClusters(rna, resolution = 0.7)
rna <- RunUMAP(rna, dims = 1:20)

# Cell annotations using SingleR
rna.sce <- as.SingleCellExperiment(rna)

#HPCA
ref.data.HPCA <- celldex::HumanPrimaryCellAtlasData()
saveRDS(ref.data.HPCA, file = "~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
#Blueprint
ref.data.BED <- celldex::BlueprintEncodeData()
saveRDS(ref.data.BED, file = "~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")

ref.data.HPCA <- readRDS("~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
ref.data.BED <- readRDS("~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")



predictions.HPCA <- SingleR(test = rna.sce, ref = ref.data.HPCA,
                            labels = ref.data.HPCA$label.main, assay.type.test = "logcounts")
predictions.BED <- SingleR(test = rna.sce, ref = ref.data.BED,
                           labels = ref.data.BED$label.main, assay.type.test = "logcounts")

# Save the annotation result
rna$SingleR.HPCA <- predictions.HPCA$pruned.labels
rna$SingleR.BED <- predictions.BED$pruned.labels

table(rna$SingleR.HPCA)
table(rna$SingleR.BED)

Idents(rna)="SingleR.BED"
DimPlot(rna,reduction = "umap",label = T)

rna$SingleR.NK.merged <- ifelse(
  grepl("NK", rna$SingleR.HPCA, ignore.case = TRUE) | 
    grepl("NK", rna$SingleR.BED, ignore.case = TRUE),
  "NK cells",
  "Other"
)
#Save the P10
saveRDS(rna, "~/data_my_scATAC/Single_RNA_data/processed_rds/p10_processed.rds")


#P11 patient
counts <- Read10X(data.dir = "~/data_my_scATAC/Single_RNA_data/GSM5276943_P11")
rna <- CreateSeuratObject(counts = counts, min.cells = 3)

# Mito percentage 
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")

# Quality control
rna@meta.data$nCount_RNA_outlier_2mad <- isOutlier(log(rna$nCount_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$nFeature_RNA_outlier_2mad <- isOutlier(log(rna$nFeature_RNA), log = FALSE, type = "lower", nmads = 2)
rna@meta.data$percent_mt_outlier_2mad <- isOutlier(log1p(rna$percent.mt), log = FALSE, type = "higher", nmads = 2)

rna <- subset(rna, subset = nCount_RNA_outlier_2mad == FALSE &
                nFeature_RNA_outlier_2mad == FALSE &
                percent_mt_outlier_2mad == FALSE)

#Standard workflow
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna <- ScaleData(rna)
rna <- RunPCA(rna)
rna <- FindNeighbors(rna, dims = 1:20)
rna <- FindClusters(rna, resolution = 0.7)
rna <- RunUMAP(rna, dims = 1:20)

# Cell annotations using SingleR
rna.sce <- as.SingleCellExperiment(rna)

#HPCA
ref.data.HPCA <- celldex::HumanPrimaryCellAtlasData()
saveRDS(ref.data.HPCA, file = "~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
#Blueprint
ref.data.BED <- celldex::BlueprintEncodeData()
saveRDS(ref.data.BED, file = "~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")

ref.data.HPCA <- readRDS("~/data_my_scATAC/Single_RNA_data/HPCA_celldex.rds")
ref.data.BED <- readRDS("~/data_my_scATAC/Single_RNA_data/BluePrintEncode_celldex.rds")



predictions.HPCA <- SingleR(test = rna.sce, ref = ref.data.HPCA,
                            labels = ref.data.HPCA$label.main, assay.type.test = "logcounts")
predictions.BED <- SingleR(test = rna.sce, ref = ref.data.BED,
                           labels = ref.data.BED$label.main, assay.type.test = "logcounts")

# Save the annotation result
rna$SingleR.HPCA <- predictions.HPCA$pruned.labels
rna$SingleR.BED <- predictions.BED$pruned.labels

table(rna$SingleR.HPCA)
table(rna$SingleR.BED)

Idents(rna)="SingleR.BED"
DimPlot(rna,reduction = "umap",label = T)

rna$SingleR.NK.merged <- ifelse(
  grepl("NK", rna$SingleR.HPCA, ignore.case = TRUE) | 
    grepl("NK", rna$SingleR.BED, ignore.case = TRUE),
  "NK cells",
  "Other"
)
#Save the P11
saveRDS(rna, "~/data_my_scATAC/Single_RNA_data/processed_rds/p11_processed.rds")



##Combine all the rds files (6 patients)

dataset_paths <- list.files(path = "~/data_my_scATAC/Single_RNA_data/processed_rds", 
                            pattern = "*_processed.rds", 
                            full.names = TRUE)

#Merge all the seurat objects
seurat_list <- lapply(dataset_paths, function(f) {
  obj <- readRDS(f)
  fname <- basename(f)
  sample_id <- str_extract(fname, "^p\\d+")
  obj$Sample_ID <- toupper(sample_id)  # Transfer to the upper letter "P"
  return(obj)
})

rna <- merge(x = seurat_list[[1]], y = seurat_list[-1])

rna$Cell_barcode <- rownames(rna@meta.data)
sum(duplicated(rna$Cell_barcode))
unique(rna$Cell_barcode)

# Standard processing for merged scRNA data
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna <- ScaleData(rna)
rna <- RunPCA(rna)
rna <- FindNeighbors(rna, dims = 1:30)
rna <- FindClusters(rna, resolution = 0.7)
rna <- RunUMAP(rna, dims = 1:30)

#Visualize the UMAP with NK and other cells
DimPlot(rna, group.by = "SingleR.NK.merged", reduction = "umap", label = TRUE)
DimPlot(
  rna,
  group.by = "SingleR.NK.merged",
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  label.size = 6,
  pt.size = 1.2,
  cols = c("#fc8d62", "#a6cee3")  
) +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.position = "right",
    legend.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 18)
  ) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  ggtitle("NK cells in merged scRNA data")

#Rename the patient ID 
rna$Simple_ID_new <- recode(rna$Sample_ID,
                            "P6" = "P1",
                            "P7" = "P2",
                            "P8" = "P3",
                            "P9" = "P4",
                            "P10" = "P5",
                            "P11" = "P6")

#Save the processed merged scRNA data
saveRDS(rna, file = "~/data_my_scATAC/Single_RNA_data/rna_merged_processed.rds")



##Extract the NK cells from the merged scRNA object
rna <- readRDS("~/data_my_scATAC/Single_RNA_data/rna_merged_processed.rds")
#Subset the NK cells 
rna_NK <- subset(rna, subset = SingleR.NK.merged == "NK cells")

#Standard processing 
rna_NK <- NormalizeData(rna_NK)
rna_NK <- FindVariableFeatures(rna_NK, selection.method = "vst", nfeatures = 2000)
rna_NK <- ScaleData(rna_NK)
rna_NK <- RunPCA(rna_NK)
rna_NK <- FindNeighbors(rna_NK, dims = 1:20)
rna_NK <- FindClusters(rna_NK, resolution = 1.5)
rna_NK <- RunUMAP(rna_NK, dims = 1:20)



#Annotation for adaptive NK cells using our NK annotation dataset

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

NKdb5_ = "~/data_my_single_cell/annotation_ref/NK_annotation_6.xlsx"
gs_list = gene_sets_prepare(NKdb5_, tissue)


es.max = sctype_score(scRNAseqData = rna_NK[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


cL_resutls = do.call("rbind", lapply(unique(rna_NK@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(rna_NK@meta.data[rna_NK@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(rna_NK@meta.data$seurat_clusters==cl)), 10)
}))

NKsctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
NKsctype_scores
NKsctype_scores$type[as.numeric(as.character(NKsctype_scores$scores)) < NKsctype_scores$ncells/4] = "Unknown"
NKsctype_scores
NKsctype_scores=as.data.frame(NKsctype_scores)

#Save the annotation results 
write.csv(NKsctype_scores,"~/data_my_scATAC/Single_RNA_data/NK_annotations.csv")

# Add annotation information into the rna_NK object
rna_NK$NK_cell_type <- case_when(
  rna_NK$seurat_clusters == 4 ~ "Adaptive_NK_cells",
  rna_NK$seurat_clusters == 5 ~ "Conventional_NK_cells",
  rna_NK$seurat_clusters == 3 ~ "Inflamed_NK_cells",
  TRUE ~ paste0("NK_cluster_", rna_NK$seurat_clusters)
)


Idents(rna_NK)="NK_cell_type"

#UMAP with NK subsets 
DimPlot(
  rna_NK,
  group.by = "NK_cell_type",     
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 1.2
) +
  theme_classic(base_size = 14) +
  labs(
    title = "UMAP of NK cell clusters",
    x = "UMAP_1",
    y = "UMAP_2",
    color = "NK_Clusters"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),  
    axis.title = element_text(size = 16),               
    axis.text = element_text(size = 14),                
    legend.title = element_text(size = 14),             
    legend.text = element_text(size = 14)               
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3)) 
  )

#UMAP without annotation for NK subsets
DimPlot(
  rna_NK,
  group.by = "RNA_snn_res.1.5",     
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 1.2
) +
  theme_classic(base_size = 14) +
  labs(
    title = "UMAP of NK cell clusters",
    x = "UMAP_1",
    y = "UMAP_2",
    color = "NK_Clusters"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),  
    axis.title = element_text(size = 16),               
    axis.text = element_text(size = 14),                
    legend.title = element_text(size = 14),             
    legend.text = element_text(size = 14)               
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3)) 
  )

#Save the rna_NK scRNA object
saveRDS(rna_NK,"~/data_my_scATAC/Single_RNA_data/rna_NK.rds")

rna_NK=readRDS("~/data_my_scATAC/Single_RNA_data/rna_NK.rds")



#Use addmodulescore to score the adaptive NK cluster (Cluster 4)
#Read our annotation geneset

Idents(rna_NK) <- "seurat_clusters"

NK_annotation <- read_xlsx("~/data_my_single_cell/annotation_ref/NK_annotation_6.xlsx")
aNK_gene <- strsplit(NK_annotation$geneSymbolmore1[1], ",")[[1]]

aNK_gene<- trimws(aNK_gene)  
aNK_gene<- unique(aNK_gene)

rna_NK <- AddModuleScore(
  object = rna_NK,
  features = list(aNK_gene),
  name = "aNK_Score"
)

#Set Cluster 4 as the reference group 
ref_group <- "4"

df <- rna_NK@meta.data %>%
  dplyr::select(aNK_Score1, seurat_clusters) %>%
  dplyr::rename(score = aNK_Score1, group = seurat_clusters) %>%
  dplyr::mutate(group = as.factor(group))

mat <- pairwise.wilcox.test(df$score, df$group, p.adjust.method = "BH")$p.value

# Initialize the list
pval_list <- list()

# One situation：ref_group  in rownames
if (ref_group %in% rownames(mat)) {
  idx <- which(!is.na(mat[ref_group, ]))
  if (length(idx) > 0) {
    pval_list[[length(pval_list) + 1]] <- data.frame(
      group1 = ref_group,
      group2 = names(idx),
      p.adj = unname(mat[ref_group, idx])
    )
  }
}

# One situation:ref_group in colnames
if (ref_group %in% colnames(mat)) {
  idx <- which(!is.na(mat[, ref_group]))
  if (length(idx) > 0) {
    pval_list[[length(pval_list) + 1]] <- data.frame(
      group1 = rownames(mat)[idx],
      group2 = ref_group,
      p.adj = unname(mat[idx, ref_group])
    )
  }
}

# Combine all the results
pval_df <- do.call(rbind, pval_list)

# Add y.position
pval_df$y.position <- seq(max(df$score) * 3, by = 0.2, length.out = nrow(pval_df))

# Transder p vaule to asterisk marks
pval_df$p.signif <- cut(pval_df$p.adj,
                        breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                        labels = c("****", "***", "**", "*", "ns"))

#Visualize the violin map showing the comparsion of aNK score across all the clusters
ggviolin(
  df, x = "group", y = "score",
  fill = "group",
  palette = cors,
  trim = FALSE,           
  draw_quantiles = NULL   
) +
  geom_jitter(
    width = 0.15,
    size = 1,
    alpha = 0.8,
    shape = 16,            
    color = "black"
  ) +
  stat_pvalue_manual(
    pval_df,
    label = "p.signif",
    y.position = "y.position",
    tip.length = 0.01,
    size = 4.5
  ) +
  ylim(NA, max(pval_df$y.position) + 0.3) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),        
    axis.line = element_line(color = "black", size = 0.4),
    axis.ticks = element_line(size = 0.3),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Adaptive NK Module Score"
  )











