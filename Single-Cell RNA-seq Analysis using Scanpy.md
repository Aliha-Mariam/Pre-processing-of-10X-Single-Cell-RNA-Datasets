# 🧬 Single-Cell RNA-seq Analysis using Scanpy


#  Project Overview

This project demonstrates a **complete single-cell RNA sequencing (scRNA-seq) analysis pipeline** using Scanpy.

## 🔬 Dataset Description

* **Biological source**: Bone marrow mononuclear cells (healthy donors)
* **Technology**: 10X Multiome (Gene Expression + Chromatin Accessibility)
* **Total Cells**: 8,785
* **Total Genes**: 36,601

## Objectives

* Perform quality control and filtering
* Detect and remove doublets
* Normalize and preprocess data
* Identify highly variable genes
* Perform dimensionality reduction (PCA, UMAP)
* Cluster cells (Leiden algorithm)
* Annotate cell types
* Identify marker genes

---


# ⚙️ Installation

```bash
pip install scanpy anndata scrublet leidenalg celltypist decoupler omnipath
```

---

# 📥 Data Loading

The dataset is downloaded and stored using **pooch**, then loaded into an **AnnData object**, which is the main data structure used in Scanpy.

## 🔹 Code

```python
import anndata as ad
import scanpy as sc
import pooch

sc.set_figure_params(dpi=50, facecolor="white")

def download_sample(sample_id, known_hash):
    path = pooch.retrieve(
        path=pooch.os_cache("scverse_tutorials"),
        url=f"https://exampledata.scverse.org/tutorials/neurips-2021/{sample_id}_filtered_feature_bc_matrix.h5",
        known_hash=known_hash,
    )
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    return adata

samples = {
    "s1d1": "md5:a99285913ea3f3d22600d3d2f8a88e34",
    "s1d3": "md5:825f7f7578e3dc0b8955f5a97a402338",
}

adatas = {k: download_sample(k, v) for k, v in samples.items()}
adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
```
This function downloads each sample securely using hash verification and loads it into AnnData format. It also ensures gene names are unique to avoid downstream conflicts.Multiple samples are downloaded and merged into a single AnnData object. Each cell is labeled by sample origin, and duplicate cell IDs are corrected for consistency.


---

#  Quality Control (QC)

Quality control ensures removal of **low-quality or damaged cells**.

### Metrics Used:

* `n_genes_by_counts` → gene diversity per cell
* `total_counts` → sequencing depth
* `pct_counts_mt` → mitochondrial content (indicator of cell stress)

## 🔹 Code

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True)
```
Violin plots visualize distributions of QC metrics across cells, helping detect outliers and poor-quality cells.

---
<img width="1502" height="476" alt="pct" src="https://github.com/user-attachments/assets/3c2a622d-eb5d-4fc6-ab42-0efba4a22fcf" />


```python
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
```

<img width="428" height="401" alt="pct2" src="https://github.com/user-attachments/assets/f84bf9cb-1e99-4117-9d30-d9b0df5adbc6" />

This scatter plot shows relationship between sequencing depth and gene diversity, with mitochondrial percentage indicating cell health.

---

# 🔍 Filtering

Removes:

* Cells with very low gene counts
* Genes expressed in very few cells

## 🔹 Code

```python
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
```

---

#  Doublet Detection

Doublets = two cells captured as one → leads to incorrect clustering.

## 🔹 Tool

* Scrublet

## 🔹 Code

```python
sc.external.pp.scrublet(adata, batch_key="sample")

adata = adata[~adata.obs["predicted_doublet"]].copy()
```

---

#  Normalization

* Normalizes sequencing depth across cells
* Applies log transformation to stabilize variance

## 🔹 Code

```python
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

---

#  Feature Selection

## 🔹 Code

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    batch_key="sample"
)
sc.pl.highly_variable_genes(adata)
```
<img width="725" height="377" alt="pct 3" src="https://github.com/user-attachments/assets/c586d4b5-bd7f-40e3-955b-7988cb929a41" />

---

This step selects genes with the highest biological variation across cells, which are most informative for downstream analysis.

# 📉 Dimensionality Reduction (PCA)
## 🔹 Code

```python
sc.tl.pca(adata)

sc.pl.pca_variance_ratio(adata,
                         log=True,
                         save="_pca_variance.png")
```
PCA reduces high-dimensional gene expression into principal components capturing major variance.
This plot shows how much variance is explained by each principal component.

<img width="374" height="404" alt="pct 4" src="https://github.com/user-attachments/assets/cf543132-1c01-41bf-ada7-58afff295ae9" />

---

# 🌐 UMAP Visualization

## 🔹 Code

```python
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.pl.umap(adata,
           color="sample",
           save="_umap_sample.png")
```
<img width="444" height="376" alt="umap" src="https://github.com/user-attachments/assets/7eb891b6-454f-4122-a2f7-2e6501abc792" />
A neighborhood graph is constructed and projected into 2D space using UMAP to visualize cell relationships.This visualizes how cells cluster based on their sample origin.

---

#  Clustering (Leiden)

Groups cells into clusters based on similarity.

## 🔹 Code

```python
sc.tl.leiden(adata)

sc.pl.umap(adata,
           color="leiden",
           save="_umap_clusters.png")
```
<img width="498" height="379" alt="clustering" src="https://github.com/user-attachments/assets/c9ada67a-e0d1-44b3-a88a-1af28ef84da0" />

Clusters are visualized in UMAP space to identify distinct cell populations.

---

#  Cell Type Annotation

---

## 🔹 1. Marker-Based Annotation

Uses known gene markers to identify cell types.

```python
marker_genes = {
    "CD14+ Mono": ["FCN1", "CD14"],
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    # Note: DMXL2 should be negative
    "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
    "Erythroblast": ["MKI67", "HBA1", "HBB"],
    # Note HBM and GYPA are negative markers
    "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
    "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
    "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
    # Note IGHD and IGHM are negative markers
    "B cells": ["MS4A1", "ITGB1", "COL4A4", "PRDM1", "IRF4", "PAX5", "BCL11A", "BLK", "IGHD", "IGHM"],
    "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
    "Plasmablast": ["XBP1", "PRDM1", "PAX5"],  # Note PAX5 is a negative marker
    "CD4+ T": ["CD4", "IL7R", "TRBC2"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],

def group_max(adata: sc.AnnData, groupby: str) -> str:
    import pandas as pd

    agg = sc.get.aggregate(adata, by=groupby, func="mean")
    return pd.Series(agg.layers["mean"].sum(1), agg.obs[groupby]).idxmax()
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res0_02")
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res0_5")

```
<img width="2227" height="748" alt="dabang" src="https://github.com/user-attachments/assets/615ecb63-2815-4116-8a6d-113ce255891a" />

---

## 🔹 2. Automatic Annotation (CellTypist)

Uses machine learning to predict cell types.

```python
import celltypist as ct

predictions = ct.annotate(
    adata,
    model="Immune_All_Low.pkl",
    majority_voting=True
)

adata = predictions.to_adata()
sc.pl.umap(adata, color="majority_voting", ncols=1)
```
<img width="625" height="376" alt="dabang2" src="https://github.com/user-attachments/assets/875308d6-6999-43af-851d-0abc6c872747" />

---

## 🔹 3. Enrichment-Based Annotation (Decoupler)

Uses marker databases (PanglaoDB) to infer cell identity.

```python
import decoupler as dc

# 1. Fetch
markers = dc.get_resource("PanglaoDB")

# 2. Debug: Look for the column that specifies Human vs Mouse
# Common columns are 'taxon', 'ncbitaxonid', 'species', or 'organism'
possible_org_cols = ['taxon', 'ncbitaxonid', 'species', 'organism']
org_col = next((c for c in possible_org_cols if c in markers.columns), None)

if org_col:
    # Filter for Human (Taxon 9606 or name 'human')
    is_human = markers[org_col].astype(str).str.lower().isin(['human', '9606', 'hsapiens'])
    
    # Check for canonical markers - if column missing, skip this filter
    if "canonical_marker" in markers.columns:
        markers = markers[is_human & (markers["canonical_marker"])]
    else:
        markers = markers[is_human]
else:
    print("Could not find an organism column. Columns present are:", markers.columns.tolist())

# 3. Deduplicate (PanglaoDB uses 'genesymbol' or 'target' for the gene name)
gene_col = 'genesymbol' if 'genesymbol' in markers.columns else 'target'
markers = markers[~markers.duplicated(["cell_type", gene_col])]

markers.head()
)
```

---

# Differential Expression Analysis

Identifies genes that distinguish clusters.

```python
sc.tl.rank_genes_groups(adata, groupby="leiden")
sc.tl.filter_rank_genes_groups(adata, min_fold_change=1.5)
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_res0_5", standard_scale="var", n_genes=5)
```
<img width="2840" height="683" alt="degs" src="https://github.com/user-attachments/assets/a4ca65c8-3f6b-430d-a416-f8ec6cd094c2" />

Top marker genes per cluster are visualized using dot plots.

---

# 📊 Final Results Visualization

```python
cluster3_genes = ["LYZ", "ACTB", "S100A6", "S100A4", "CST3"]
sc.pl.umap(adata, color=[*cluster3_genes, "leiden_res0_5"], legend_loc="on data", frameon=False, ncols=3)
)
```
<img width="1179" height="710" alt="final" src="https://github.com/user-attachments/assets/39f0ae7b-e4d3-4cb4-b289-148046b287f1" />

Final UMAP shows clearly separated biological clusters representing different immune cell populations.

---

# Common Issues

* Duplicate gene names → fixed
* AnnData warnings → safe
* Dependency conflicts may occur
* CellTypist requires proper normalization

---

#  Key Results

* High-quality dataset obtained
* Doublets removed
* Clear clustering achieved
* Cell types identified:

  * B cells
  * T cells
  * NK cells
  * Monocytes
* Marker genes successfully identified

---

#  Conclusion

This pipeline successfully:

* Processed scRNA-seq data
* Identified biologically meaningful clusters
* Validated results using multiple annotation approaches

---

#  Future Work

* Batch correction
* Sub-clustering
* Pathway analysis
* Multi-omics integration

