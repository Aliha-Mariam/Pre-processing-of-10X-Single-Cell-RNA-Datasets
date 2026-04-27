# 🧬 Single-Cell RNA-seq Analysis using Scanpy

---

# 📌 Project Overview

This project demonstrates a **complete single-cell RNA sequencing (scRNA-seq) analysis pipeline** using Scanpy.

## 🔬 Dataset Description

* **Biological source**: Bone marrow mononuclear cells (healthy donors)
* **Technology**: 10X Multiome (Gene Expression + Chromatin Accessibility)
* **Total Cells**: 8,785
* **Total Genes**: 36,601

## 🎯 Objectives

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

## 🔹 Explanation

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

## ⚠️ Notes

* Duplicate gene names and cell IDs are corrected to avoid downstream errors

---

# 🧪 Quality Control (QC)

## 🔹 Explanation

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

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True)
```

## 📊 Add Figures (IMPORTANT)

```python
sc.pl.violin(adata,
             ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
             save="_qc_violin.png")

sc.pl.scatter(adata,
              "total_counts",
              "n_genes_by_counts",
              color="pct_counts_mt",
              save="_qc_scatter.png")
```

### 👉 Add in README:

```
## 📊 Quality Control Plots
![QC Violin](figures/qc_violin.png)
![QC Scatter](figures/qc_scatter.png)
```

---

# 🔍 Filtering

## 🔹 Explanation

Removes:

* Cells with very low gene counts
* Genes expressed in very few cells

## 🔹 Code

```python
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)
```

---

# ⚠️ Doublet Detection

## 🔹 Explanation

Doublets = two cells captured as one → leads to incorrect clustering.

## 🔹 Tool

* Scrublet

## 🔹 Code

```python
sc.external.pp.scrublet(adata, batch_key="sample")

adata = adata[~adata.obs["predicted_doublet"]].copy()
```

---

# 📊 Normalization

## 🔹 Explanation

* Normalizes sequencing depth across cells
* Applies log transformation to stabilize variance

## 🔹 Code

```python
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

---

# 🎯 Feature Selection

## 🔹 Explanation

Selects **highly variable genes (HVGs)** → most informative genes

## 🔹 Code

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    batch_key="sample"
)
```

---

# 📉 Dimensionality Reduction (PCA)

## 🔹 Explanation

Reduces high-dimensional gene space into principal components.

## 🔹 Code

```python
sc.tl.pca(adata)

sc.pl.pca_variance_ratio(adata,
                         log=True,
                         save="_pca_variance.png")
```

## 📊 Add Figure

```
## 📊 PCA Variance
![PCA](figures/pca_variance.png)
```

---

# 🌐 UMAP Visualization

## 🔹 Explanation

UMAP projects cells into 2D space preserving structure.

## 🔹 Code

```python
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.pl.umap(adata,
           color="sample",
           save="_umap_sample.png")
```

## 📊 Add Figure

```
## 🌐 UMAP Visualization
![UMAP](figures/umap_sample.png)
```

---

# 🔗 Clustering (Leiden)

## 🔹 Explanation

Groups cells into clusters based on similarity.

## 🔹 Code

```python
sc.tl.leiden(adata)

sc.pl.umap(adata,
           color="leiden",
           save="_umap_clusters.png")
```

## 📊 Add Figure

```
## 🔗 Clustering
![Clusters](figures/umap_clusters.png)
```

---

# 🧬 Cell Type Annotation

---

## 🔹 1. Marker-Based Annotation

### Explanation

Uses known gene markers to identify cell types.

```python
marker_genes = {
    "B cells": ["MS4A1"],
    "T cells": ["CD4", "CD8A"],
    "NK cells": ["GNLY", "NKG7"],
}

sc.pl.dotplot(adata,
              marker_genes,
              groupby="leiden",
              save="_marker_dotplot.png")
```

---

## 🔹 2. Automatic Annotation (CellTypist)

### Explanation

Uses machine learning to predict cell types.

```python
import celltypist as ct

predictions = ct.annotate(
    adata,
    model="Immune_All_Low.pkl",
    majority_voting=True
)

adata = predictions.to_adata()
```

---

## 🔹 3. Enrichment-Based Annotation (Decoupler)

### Explanation

Uses marker databases (PanglaoDB) to infer cell identity.

```python
import decoupler as dc

markers = dc.get_resource("PanglaoDB")

dc.run_mlm(
    adata,
    net=markers.rename(columns=dict(cell_type="source", genesymbol="target")),
    weight=None
)
```

---

# 🧪 Differential Expression Analysis

## 🔹 Explanation

Identifies genes that distinguish clusters.

```python
sc.tl.rank_genes_groups(adata, groupby="leiden")
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, save="_deg.png")
```

---

# 📊 Final Results Visualization

```python
sc.pl.umap(
    adata,
    color=["leiden", "majority_voting"],
    save="_final_annotation.png"
)
```

---

# ⚠️ Common Issues

* Duplicate gene names → fixed
* AnnData warnings → safe
* Dependency conflicts may occur
* CellTypist requires proper normalization

---

# ✅ Key Results

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

# 🚀 Conclusion

This pipeline successfully:

* Processed scRNA-seq data
* Identified biologically meaningful clusters
* Validated results using multiple annotation approaches

---

# 📌 Future Work

* Batch correction
* Sub-clustering
* Pathway analysis
* Multi-omics integration

