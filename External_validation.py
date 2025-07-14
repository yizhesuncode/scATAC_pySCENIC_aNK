#Import all the dependencies 

import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator
import matplotlib as mpl

#Read the processed NK data from Netskar H, et al., Nat Immunol, 2024
sce_all_NK = sc.read_h5ad(r"C:\\Users\\sunyi\\Desktop\\adata_all_nk_after_mapping.h5ad")


##Show the MTFP1 expression across all the NK clusters


#Set the font type as Arial
mpl.rcParams['font.family'] = 'Arial'

# Extract the expression of MTFP1, perform log1p transferring
mtfp1_expr = sce_all_NK[:, "MTFP1"].X
if not isinstance(mtfp1_expr, np.ndarray):
    mtfp1_expr = mtfp1_expr.toarray().flatten()
mtfp1_expr_log = np.log1p(mtfp1_expr)

# Construct the dataFrame
df = pd.DataFrame({
    "MTFP1": mtfp1_expr_log,
    "subset": sce_all_NK.obs["subset"].values
})

# Order the clusters by the expression of MTFP1 (high to low)
subset_order = (
    df.groupby("subset")["MTFP1"]
    .mean()
    .sort_values(ascending=False)
    .index.tolist()
)
df["subset"] = pd.Categorical(df["subset"], categories=subset_order, ordered=True)

# Set the comparsion groups (Adaptive vs. others)
target_group = "Adaptive"
comparisons = [(target_group, group) for group in subset_order if group != target_group]

#Visualize and beautify the plot
palette = sns.color_palette("tab20", n_colors=len(subset_order))
sns.set_context("notebook", font_scale=1.5)

plt.figure(figsize=(14, 9))
ax = sns.barplot(data=df, x="subset", y="MTFP1", ci="sd", palette=palette, edgecolor="black")
sns.stripplot(data=df, x="subset", y="MTFP1", color="black", size=2, jitter=True, alpha=0.4)

# Mark the significance
annotator = Annotator(ax, pairs=comparisons, data=df, x="subset", y="MTFP1")
annotator.configure(test="Mann-Whitney", text_format="star", comparisons_correction=None)
annotator.apply_and_annotate()

ax.set_title(r"$\it{MTFP1}$ Expression in NK Cell Subsets (vs Adaptive)", fontsize=22, fontweight="bold")
ax.set_ylabel("log1p(MTFP1)", fontsize=22)
ax.set_xlabel("NK Cell Subset", fontsize=22)
ax.tick_params(axis='x', labelrotation=45, labelsize=22)
ax.tick_params(axis='y', labelsize=22)
sns.despine()

plt.tight_layout()
plt.show()





##Show the CRCP expression across all the NK clusters

#Set the font type as Arial
mpl.rcParams['font.family'] = 'Arial'

#Extract the expression of MTFP1, perform log1p transferring
crcp_expr = sce_all_NK[:, "CRCP"].X
if not isinstance(crcp_expr, np.ndarray):
    crcp_expr = crcp_expr.toarray().flatten()
crcp_expr_log = np.log1p(crcp_expr)

# Dataframe construction
df = pd.DataFrame({
    "CRCP": crcp_expr_log,
    "subset": sce_all_NK.obs["subset"].values
})

# Order the clusters by the expression of MTFP1 (high to low)
subset_order = (
    df.groupby("subset")["CRCP"]
    .mean()
    .sort_values(ascending=False)
    .index.tolist()
)
df["subset"] = pd.Categorical(df["subset"], categories=subset_order, ordered=True)

# Set the comparison groups (TiCD56dim vs others)
target_group = "TrCD56dim"
comparisons = [(target_group, group) for group in subset_order if group != target_group]

# Visualize and beautify the plot
palette = custom_palette
sns.set_context("notebook", font_scale=1.5)

plt.figure(figsize=(14, 9))
ax = sns.barplot(data=df, x="subset", y="CRCP", ci="sd", palette=palette, edgecolor="black")
sns.stripplot(data=df, x="subset", y="CRCP", color="black", size=2, jitter=True, alpha=0.4)

# Mark the significance 
annotator = Annotator(ax, pairs=comparisons, data=df, x="subset", y="CRCP")
annotator.configure(test="Mann-Whitney", text_format="star", comparisons_correction=None)
annotator.apply_and_annotate()

ax.set_title(r"$\it{CRCP}$ Expression in NK Cell Subsets (vs TrCD56dim)", fontsize=22, fontweight="bold")
ax.set_ylabel("log1p(CRCP)", fontsize=22)
ax.set_xlabel("NK Cell Subset", fontsize=22)
ax.tick_params(axis='x', labelrotation=45, labelsize=22)
ax.tick_params(axis='y', labelsize=22)

sns.despine()

plt.tight_layout()
plt.show()



###End



