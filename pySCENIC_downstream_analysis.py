
#import dependencies
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import operator as op
import matplotlib.pyplot as plt
import seaborn as sns
import loompy as lp
import json
import base64
import zlib
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from MulticoreTSNE import MulticoreTSNE as TSNE
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
import seaborn as sns
from pyscenic.binarization import binarize
from adjustText import adjust_text



sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.set_figure_params(dpi=300, fontsize=10, dpi_save=600)

#Construct adata in python
#Build a function which can transfer the seurat file to adata
def seurat_to_adata(counts,#counts.mtx
                   meta,#metadata.csv
                   gene_name,#gene_names.csv
                   pca,#pca.csv
                   reduction1,#UMAP_1
                   reduction2):#UMAP_2
    # Load the count matrix
    X = io.mmread(counts)
    # Create adata
    adata = anndata.AnnData(X=X.transpose().tocsr())
    cell_meta = pd.read_csv(meta)#metadata
    with open(gene_name, 'r') as f:gene_names = f.read().splitlines()
    adata.obs = cell_meta
    adata.obs.index = adata.obs['barcode']
    adata.var.index = gene_names
    pca = pd.read_csv(pca)
    pca.index = adata.obs.index
    adata.obsm['X_pca'] = pca.to_numpy()
    adata.obsm['X_tsne'] = np.vstack((adata.obs[reduction1].to_numpy(), adata.obs[reduction2].to_numpy())).T
    return adata

#Run the function using our files generated from R
IPLA_NK_merge = seurat_to_adata(
    counts=os.path.expanduser('~/IPLA_tumor_NK_Pyscenic/seurat_counts.mtx'),
    meta=os.path.expanduser('~/IPLA_tumor_NK_Pyscenic/seurat_metadata.csv'),
    gene_name=os.path.expanduser('~/IPLA_tumor_NK_Pyscenic/seurat_gene_names.csv'),
    pca=os.path.expanduser('~/IPLA_tumor_NK_Pyscenic/seurat_pca.csv'),
    reduction1='UMAP_1',
    reduction2='UMAP_2'
)


# Set the overall style
sc.set_figure_params(dpi=150, frameon=False, figsize=(5, 5))
sc.settings.set_figure_params(dpi_save=300)  

# Customize the colors
custom_palette = sns.color_palette("tab10", n_colors=10)

#Save the adata file for IPLA object
IPLA_NK_merge.obs = IPLA_NK_merge.obs.astype(str)
IPLA_NK_merge.write('IPLA_NK_merge.h5ad')


# Read the loom file generated from PySCENIC workflow

IPLA_NK_SCENIC = '/disk2/user/yizhsu/IPLA_tumor_NK_Pyscenic/IPLA_NK_merge_SCENIC.loom'

# Extract the AUC for TFs
lf = lp.connect(IPLA_NK_SCENIC, mode='r+', validate=False )
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()
auc_mtx

#Calculate RSS
rss_cellType = regulon_specificity_scores(auc_mtx, IPLA_NK_merge.obs.Cluster_final)
rss_cellType.to_csv("IPLA_NK_rss_clusters.csv", index=True)


celltype=["Cluster_0","Cluster_1","Cluster_2","Cluster_3",
          "Cluster_4","Cluster_5","Cluster_6","Cluster_7","Cluster_8"]


#Visualize the top5 TF in each cluster
fig = plt.figure(figsize=(15, 13)) 

for c, num in zip(celltype, range(1, len(celltype)+1)):
    x = rss_cellType.T[c]
    ax = fig.add_subplot(3, 3, num)  # 图的排列
    
    plot_rss(rss_cellType, c, top_n=5, max_n=None, ax=ax)

    ax.set_ylim(x.min() - (x.max() - x.min()) * 0.05,
                x.max() + (x.max() - x.min()) * 0.08)  
    for t in ax.texts:
        t.set_fontsize(11)  

    ax.set_ylabel('')
    ax.set_xlabel('')

    adjust_text(
        ax.texts,
        autoalign='xy',
        ha='right',
        va='bottom',
        only_move={'points': 'y', 'text': 'xy'},
        expand_text=(1.05, 1.1),  #
        arrowprops=dict(arrowstyle='-', color='black', lw=0.4),
        precision=0.001
    )
    
fig.supxlabel('Regulon', fontsize='x-large')
fig.supylabel('Regulon specificity score (RSS)', fontsize='x-large')

plt.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.08)

plt.rcParams.update({
    'figure.autolayout': False,  
    'figure.titlesize': 'large',
    'axes.labelsize': 'medium',
    'axes.titlesize': 'large',
    'xtick.labelsize': 'medium',
    'ytick.labelsize': 'medium'
})

plt.savefig("rank_RSS_top5_2.pdf", dpi=600, bbox_inches="tight")
plt.show()


#Visualize the RSS matrix by heatmap
#Calculate the Z-score
file_path = os.path.expanduser("~/IPLA_tumor_NK_Pyscenic/IPLA_NK_rss_clusters.csv")
rss_cellType = pd.read_csv(file_path, index_col=0)


rss_cellType_Z = pd.DataFrame(index=rss_cellType.index )
for col in list(rss_cellType.columns):
    rss_cellType_Z[col] = ( rss_cellType[col] - rss_cellType[col].mean()) / rss_cellType[col].std(ddof=0)
    rss_cellType_Z.sort_index(inplace=True)

rss_cellType_Z
rss_cellType_Z.to_csv("IPLA_NK_rss_clusters_SD.csv", index=True)

file_path = os.path.expanduser("~/IPLA_tumor_NK_Pyscenic/IPLA_NK_rss_clusters_SD.csv")
rss_cellType_Z =pd.read_csv(file_path, index_col=0)




#RSS Heatmap for all the differentially enriched TFs

colors = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#660099"
]
#Read RSS file with z-score
file_path = os.path.expanduser("~/IPLA_tumor_NK_Pyscenic/IPLA_NK_rss_clusters_SD.csv")
rss_cellType_Z =pd.read_csv(file_path, index_col=0)

#Read RSS file
file_path_2 = os.path.expanduser("~/IPLA_tumor_NK_Pyscenic/IPLA_NK_rss_clusters.csv")
rss_cellType = pd.read_csv(file_path_2, index_col=0)

#
sns.set(font_scale=1.2)
# 
cluster_list = rss_cellType_Z.index.tolist()
# Color projection to clusters
colors = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#660099"
]
cluster_palette = dict(zip(cluster_list, colors))  
cluster_info = pd.Series(rss_cellType_Z.index, index=rss_cellType_Z.index)
row_colors = cluster_info.map(cluster_palette)
sns.set(font_scale=1.2)
g = sns.clustermap(
    rss_cellType_Z,
    annot=False,
    square=False,
    linecolor='black',
    yticklabels=True,       
    xticklabels=False,       # TF not shown
    vmin=-2, vmax=2,
    row_colors=row_colors,
    cmap="RdBu_r",
    figsize=(10, 4),
    row_cluster=False,
    cbar_pos=(0.02, 0.8, 0.02, 0.18)  # 调整 colorbar 位置（避免压右边标签）
)
# Adjustment
g.fig.subplots_adjust(right=0.93)
plt.tight_layout()
plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=11)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=10)
plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=11)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.show()


#Visualize the 21 selectd TFs RSS by heatmap
#The 21 TFs
#Read RSS file with z-score
file_path = os.path.expanduser("~/IPLA_tumor_NK_Pyscenic/IPLA_NK_rss_clusters_SD.csv")
rss_cellType_Z =pd.read_csv(file_path, index_col=0)

TF_plot = ["PRDM1(+)","SOX6(+)","TBX3(+)","MEIS1(+)","RARA(+)",
          "SP1(+)","STAT2(+)","TCF12(+)","SP2(+)","PURA(+)",
          "SP3(+)","SPI1(+)","ZFX(+)","SP4(+)","BCL11A(+)",
          "EGR1(+)","ZNF281(+)","RREB1(+)","CTCF(+)","FOS(+)",
          "JUNB(+)"]
          
 # Exclude the possible space
TF_plot = [r.strip() for r in TF_plot]
 
# 
cluster_list = rss_cellType_Z.index.tolist()

colors = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#660099"
]
cluster_palette = dict(zip(cluster_list, colors))  
cluster_info = pd.Series(rss_cellType_Z.index, index=rss_cellType_Z.index)
row_colors = cluster_info.map(cluster_palette)


sns.set(font_scale=1.2)
g = sns.clustermap(
    rss_cellType_Z[TF_plot],
    annot=False,
    square=False,
    linecolor='black',
    yticklabels=True,       
    xticklabels=True,       
    vmin=-2, vmax=2,
    row_colors=row_colors,
    row_cluster=False,
    cmap="RdBu_r",
    figsize=(10, 4),
    cbar_pos=(0.02, 0.8, 0.02, 0.18)  
)

g.fig.subplots_adjust(right=0.93)
plt.tight_layout()
plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=11)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=10)
plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=11)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.show()


#Heatmap for 21 TFs in all the cells 

TF = ["PRDM1(+)","SOX6(+)","TBX3(+)","MEIS1(+)","RARA(+)",
          "SP1(+)","STAT2(+)","TCF12(+)","SP2(+)","PURA(+)",
          "SP3(+)","SPI1(+)","ZFX(+)","SP4(+)","BCL11A(+)",
          "EGR1(+)","ZNF281(+)","RREB1(+)","CTCF(+)","FOS(+)",
          "JUNB(+)"]#选择需要呈现的TF
celltype=["Cluster_0","Cluster_1","Cluster_2","Cluster_3",
          "Cluster_4","Cluster_5","Cluster_6","Cluster_7","Cluster_8"]


IPLA_NK_SCENIC = '/disk2/user/yizhsu/IPLA_tumor_NK_Pyscenic/IPLA_NK_merge_SCENIC.loom'
lf = lp.connect(IPLA_NK_SCENIC, mode='r+', validate=False )
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()
auc_mtx

#z-score for each cells
auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
auc_mtx_Z = (auc_mtx - auc_mtx.mean()) / auc_mtx.std(ddof=0)


colors = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#660099"
]
colorsd = dict( zip( celltype, colors ))
colormap = [ colorsd[x] for x in IPLA_NK_merge.obs.Cluster_final]

sns.set(font_scale=1.2)
g = sns.clustermap(auc_mtx_Z[TF], annot=False,  square=False,  linecolor='black',
    yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=colormap,
    cmap="YlGnBu", figsize=(21,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')

g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xticklabels(),
    fontsize=30,
    rotation=90,
    ha='center'
)

plt.savefig("/disk2/user/yizhsu/IPLA_tumor_NK_Pyscenic/figures/TF-heatmap-overlapped_ATAC_all_cells.pdf", dpi=600, bbox_inches = "tight")




#Visualize the two TF regulons STAT2 and PRDM1 in UMAP

IPLA_NK_SCENIC = '/disk2/user/yizhsu/IPLA_tumor_NK_Pyscenic/IPLA_NK_merge_SCENIC.loom'
lf = lp.connect(IPLA_NK_SCENIC, mode='r+', validate=False )
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()
auc_mtx

IPLA_NK_merge = seurat_to_adata(
    counts=os.path.expanduser('~/IPLA_tumor_NK_Pyscenic/seurat_counts.mtx'),
    meta=os.path.expanduser('~/IPLA_tumor_NK_Pyscenic/seurat_metadata.csv'),
    gene_name=os.path.expanduser('~/IPLA_tumor_NK_Pyscenic/seurat_gene_names.csv'),
    pca=os.path.expanduser('~/IPLA_tumor_NK_Pyscenic/seurat_pca.csv'),
    reduction1='UMAP_1',
    reduction2='UMAP_2'
)

#Merge AUC in scRNA object 
add_scenic_metadata(IPLA_NK_merge, auc_mtx)# Add AUC in scRNA object 

#Choose the TFs we would like to show
sc.set_figure_params(frameon=False, dpi=150, fontsize=8)
sc.pl.tsne(IPLA_NK_merge, color=['Regulon(STAT2(+))','Regulon(PRDM1(+))','Cluster_final'], ncols=3, cmap = "RdPu")

clusters = IPLA_NK_merge.obs['Cluster_final'].unique().tolist()
sorted_clusters = sorted(clusters, key=lambda x: int(x.split('_')[1]))

colors = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#660099"
]

cluster_palette = dict(zip(sorted_clusters, colors))
sc.pl.tsne(
    IPLA_NK_merge,
    color=['Regulon(PRDM1(+))','Regulon(STAT2(+))','Cluster_final'],
    ncols=3,
    cmap="RdPu",
    vmax=[0.18, 0.3, None],
    palette=cluster_palette,
    size=35,
    show=False  # 👈
)

fig = plt.gcf()

for ax in fig.axes:
    ax.title.set_fontsize(12)              
    ax.tick_params(labelsize=10)           
    legend = ax.get_legend()
    if legend is not None:
        for text in legend.get_texts():
            text.set_fontsize(12)          
        legend.set_title(legend.get_title().get_text(), prop={'size': 13})  

plt.tight_layout()
plt.show()





#Visualize the selected TF regulons by violin map
df_obs = IPLA_NK_merge.obs
signature_column_names = list(df_obs.select_dtypes('number').columns)
signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
df_scores = df_obs[signature_column_names + ['Cluster_final']]
#Calulate the Z-score
df_results = ((df_scores.groupby(by='Cluster_final').mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})

regulon_to_show = ['Regulon(PRDM1(+))', 'Regulon(SOX6(+))', 'Regulon(TBX3(+))','Regulon(MEIS1(+))',
                   'Regulon(RARA(+))','Regulon(SP1(+))','Regulon(STAT2(+))','Regulon(TCF12(+))','Regulon(SP2(+))','Regulon(PURA(+))',
                   'Regulon(SP3(+))','Regulon(SPI1(+))','Regulon(ZFX(+))','Regulon(SP4(+))','Regulon(BCL11A(+))',
                   'Regulon(EGR1(+))','Regulon(ZNF281(+))','Regulon(RREB1(+))','Regulon(CTCF(+))','Regulon(FOS(+))',
                   'Regulon(JUNB(+))']


# Retain the results only from the selected regulons
df_filtered = df_results[df_results['regulon'].isin(regulon_to_show)]

# Create the heatmao 
df_heatmap = pd.pivot_table(
    data=df_filtered.sort_values('Z', ascending=False),
    index='Cluster_final', columns='regulon', values='Z'
)

lf = lp.connect(IPLA_NK_SCENIC, mode='r+', validate=False )
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()
auc_mtx
auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
auc_mtx_Z = (auc_mtx - auc_mtx.mean()) / auc_mtx.std(ddof=0)

aucell_adata = sc.AnnData(X=auc_mtx_Z.sort_index())
aucell_adata.obs = df_obs

plt.rcParams.update({
    'xtick.labelsize': 10,   
    'ytick.labelsize': 10,   
    'axes.labelsize': 10,    
    'font.size': 10          
})

sc.pl.stacked_violin(
    aucell_adata,
    df_heatmap.columns,
    groupby='Cluster_final',
    cmap='PRGn',
    stripplot=False,
    dendrogram=False,
    rotation=45,  # 
    figsize=(len(df_heatmap.columns) * 0.5, 6),
    show=False
)

plt.xticks(rotation=45, ha='right')     
plt.subplots_adjust(bottom=0.35)        
plt.tight_layout()
plt.show()








