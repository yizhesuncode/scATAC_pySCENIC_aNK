#Prepare all the packages below for the following running

library(cowplot) 
library(ggplotify)
library(reshape2)
library(scales)
library(stats)
library(rlang)
library(viridis)
library(devtools)
library(stringr)
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggpubr)
library(Seurat)
library(dplyr)
library(DoubletFinder)
library(shinythemes)
library(ggsci)
library(immunedeconv)
library(tibble)
library(org.Hs.eg.db)
library(RcisTarget)
library(ComplexHeatmap)
library(ComplexHeatmap)
library(infercnv)
library(devtools)
library(scHCL)
library(harmony)
library(SCopeLoomR)
library(igraph)
library(readxl) 
library(dittoSeq)
library(ggpubr)
library(HGNChelper)

#Read the seurat object with already clustered from 7 ovarian patients
#Standard workflow using SCTranform; RunPCA; RunHarmony;RunUMAP,FindNeighbors, FindClusters

IPLA_NK_merge=readRDS("~/data_my_single_cell/Rds_file/IPLA_NK_merge.rds")
Idents(IPLA_NK_merge) <- "Cluster_final"
cors <- dittoColors(13) 
#Visualize the UMAP 
DimPlot(
  IPLA_NK_merge,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  cols = cors
) +
  labs(
    title = "UMAP of NK Cell Clusters from Smartseq3"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 15),       
    axis.text = element_text(size = 12),       
    legend.text = element_text(size = 14),     
    legend.title = element_text(size = 14),     
    plot.title = element_text(hjust = 0.5, size = 16),  
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  xlab("UMAP_1") + ylab("UMAP_2")


#NK cell annotation using scType 

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

#Load our annotation dataset
NKdb5_ = "~/data_my_single_cell/annotation_ref/NK_annotation_6.xlsx"
tissue = "Immune system"

gs_list = gene_sets_prepare(NKdb5_, tissue)

es.max = sctype_score(scRNAseqData = IPLA_NK_merge[["SCT"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


cL_resutls = do.call("rbind", lapply(unique(IPLA_NK_merge@meta.data$Cluster_final), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(IPLA_NK_merge@meta.data[IPLA_NK_merge@meta.data$Cluster_final==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(IPLA_NK_merge@meta.data$Cluster_final==cl)), 10)
}))

NKsctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
NKsctype_scores$type[as.numeric(as.character(NKsctype_scores$scores)) < NKsctype_scores$ncells/4] = "Unknown"
NKsctype_scores
results=as.data.frame(NKsctype_scores)

#Save the annotation results
write.csv(results,"~/IPLA_tumor_NK_Pyscenic/csv_files/Clusters_annotation_merged_aNK_annotation_IPLA.csv")

#Use addmodulescore to score the adaptive NK cluster (Cluster 1)
#Read our annotation geneset

NK_annotation=read_xlsx("~/data_my_single_cell/annotation_ref/NK_annotation_6.xlsx")
aNK_gene <- strsplit(NK_annotation$geneSymbolmore1[1], ",")[[1]]

aNK_gene<- trimws(aNK_gene)  
aNK_gene<- unique(aNK_gene)

IPLA_NK_merge <- AddModuleScore(
  object = IPLA_NK_merge,
  features = list(aNK_gene),
  name = "aNK_Score"
)
#Set Cluster 1 as the reference group 
ref_group <- "Cluster_1"
df <- IPLA_NK_merge@meta.data %>%
  dplyr::select(aNK_Score1, Cluster_final) %>%
  dplyr::rename(score = aNK_Score1, group = Cluster_final) %>%
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
pval_df$y.position <- seq(max(df$score) * 1.05, by = 0.2, length.out = nrow(pval_df))

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




###########Prepare the Pyscenic-required input ###########

#Extract the expression matrix 
Idents(IPLA_NK_merge) <- "Cluster_final"
write.csv(t(as.matrix(IPLA_NK_merge@assays$SCT@counts)),file = "~/IPLA_tumor_NK_Pyscenic/IPLA_NK_exp.csv")

#Enter into the pyscenic miniconda environment
library(reticulate)
use_condaenv("pyscenic", required = TRUE)
###########Start the pyscenic in Linux server using Shell script 

###########Prepare the files for pyscenic downstream analysis 

#Build the function to prepare multiple files required for the pyscenic downstream analysis in Python
ct
seurat_to_adata <- function(object,#seurat object
                            Dimension=c('UMAP','TSNE'),#Types of reduction
                            path){#Save directory
  seurat_obj <- object
  seurat_obj$barcode <- colnames(seurat_obj)
  if(Dimension=='UMAP'){
    cell.embeddings<- seurat_obj@reductions$umap@cell.embeddings
    seurat_obj$UMAP_1 <- cell.embeddings[,1]
    seurat_obj$UMAP_2 <- cell.embeddings[,2]
  }else{
    
    cell.embeddings<- seurat_obj@reductions$tsne@cell.embeddings
    seurat_obj$TSNE_1 <- cell.embeddings[,1]
    seurat_obj$TSNE_2 <- cell.embeddings[,2]
  }
  
  #Save metadata
  write.csv(seurat_obj@meta.data, file=paste0(path,'metadata.csv'), quote=F, row.names=F)
  #Save counts matrix
  counts_matrix <- GetAssayData(seurat_obj, assay='SCT', slot='counts')
  writeMM(counts_matrix, file=paste0(path, 'counts.mtx'))
  #Save PCA reduction
  write.csv(seurat_obj@reductions$pca@cell.embeddings, file=paste0(path,'pca.csv'), quote=F,row.names=F)
  
  #Save gene name
  write.table(data.frame('gene'=rownames(counts_matrix)),file=paste0(path,'gene_names.csv'),
              quote=F,row.names=F,col.names=F)
}

  

#Run the function seurat_to_adata
seurat_to_adata(IPLA_NK_merge,Dimension='UMAP',path = '~/IPLA_tumor_NK_Pyscenic/seurat_')

#
library(reticulate)
use_condaenv("pyscenic", required = TRUE)



##Visualize the Violin map comparing the CRCP and MTFP1 between cluster 1 and other clusters
#CRCP expression
#Create two groups：Cluster_1 vs Others
IPLA_NK_merge$CRCP_group <- ifelse(IPLA_NK_merge$Cluster_final == "Cluster_1", "Cluster_1", "Other")

# Extract the expression matrix (normalized data)
df_expr <- FetchData(IPLA_NK_merge, vars = c("CRCP", "CRCP_group"), slot = "data", assay = "SCT")

nice_colors <- c("Cluster_1" = "#fc8d62", "Other" = "#a6cee3")

# Visualize the comparison by violin map
ggplot(df_expr, aes(x = CRCP_group, y = CRCP, fill = CRCP_group)) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    color = "gray30",        
    linewidth = 0.6,
    alpha = 0.85            
  ) +
  geom_jitter(
    width = 0.15,
    height = 0.05,
    shape = 21,              
    color = "black",
    fill = "black",
    size = 1.6,
    stroke = 0.3,
    alpha = 0.5
  ) +
  stat_compare_means(
    method = "wilcox.test", #Wilcoxon T test was used to compare
    comparisons = list(c("Cluster_1", "Other")),
    label = "p.signif",
    label.y = max(df_expr$CRCP, na.rm = TRUE) * 1.08,
    tip.length = 0.02,
    size = 6
  ) +
  scale_fill_manual(values = nice_colors) +
  theme_classic(base_size = 16) +
  labs(
    x = "",
    y = expression(italic("CRCP") * " Expression"),
    title = expression(italic("CRCP") * " Expression")
  ) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  )


#MTFP1 expression
# Cluster_1 vs Other
df_expr <- FetchData(IPLA_NK_merge, vars = c("MTFP1", "Cluster_final"))
df_expr$Group <- ifelse(df_expr$Cluster_final == "Cluster_1", "Cluster_1", "Other")


nice_colors <- c("Cluster_1" = "#fc8d62", "Other" = "#a6cee3")

ggplot(df_expr, aes(x = Group, y = MTFP1, fill = Group)) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    color = "gray30",
    linewidth = 0.6,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.6,
    color = "black",
    alpha = 0.6
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Cluster_1", "Other")),
    label = "p.signif",
    label.y = max(df_expr$MTFP1, na.rm = TRUE) * 1.08,
    tip.length = 0.02,
    size = 6
  ) +
  scale_fill_manual(values = nice_colors) +
  theme_classic(base_size = 16) +
  labs(
    x = "",
    y = expression(italic("MTFP1") * " Expression"),
    title = expression(italic("MTFP1") * " Expression")
  ) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  )








