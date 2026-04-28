# 🧬 Single-Cell RNA-seq Analysis using AnnData & Scanpy

##  Overview

This project demonstrates how to process and explore single-cell RNA sequencing (scRNA-seq) data using **AnnData** and **Scanpy**. It covers data loading, preprocessing, visualization, and structural understanding of scRNA-seq datasets.

---

##  Step 1: Import Libraries

```python
import anndata
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pooch
```
All required libraries are imported to set up the environment:

- `anndata`: handles `.h5ad` structured datasets  
- `scanpy`: single-cell RNA-seq analysis tools  
- `numpy`: numerical computation  
- `matplotlib`: plotting and visualization  
- `pooch`: reproducible dataset downloading  

---

##  Step 2: Download Dataset

```python
datapath = pooch.retrieve(
    path=pooch.os_cache("scverse_tutorials"),
    url="https://exampledata.scverse.org/tutorials/scverse-getting-started-anndata-pbmc3k_processed.h5ad",
    known_hash="md5:b80deb0997f96b45d06f19c694e46243",
)
```
The dataset is downloaded from an online repository and stored locally. The MD5 hash verification ensures that the file is not corrupted or modified. This guarantees reproducibility of results.
 

---

## Step 3: Load AnnData Object

```python
adata = anndata.read_h5ad(datapath)
adata

```
The dataset is loaded into an AnnData object which organizes gene expression data, metadata, and embeddings into a structured format. The output shows number of cells and genes confirming successful loading.
<img width="616" height="155" alt="sc1" src="https://github.com/user-attachments/assets/c68f6515-d5c3-4d50-89a9-f77a0d84a58f" />

  

---

## Step 4: Expression Matrix

```python
adata.X
```

- Rows = cells  
- Columns = genes  
- Stores log-normalized gene expression  
- Sparse format improves memory efficiency  

---

##  Sparse Matrix Inspection

```python
print(adata.X.data)
print(adata.X.indices)
print(adata.X.nnz / np.prod(adata.X.shape))
```

- Extracts non-zero values  
- Shows indices of expression values  
- Calculates sparsity of dataset  

---

## Step 5: Layers

```python
adata.layers
```

- Stores multiple versions of data  
- Example: raw counts, normalized counts  

---

##  CPM Normalization Layer

```python
adata.layers["counts_per_million"] = adata.layers["raw"].copy()

sc.pp.normalize_total(
    adata,
    target_sum=10**6,
    layer="counts_per_million"
)
```

- Creates CPM-normalized dataset  
- Adjusts sequencing depth differences  
- Preserves raw data separately  

---

##  Step 6: Gene Expression Visualization

```python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3))
genes_of_interest = ["CD8A", "CD4", "KLRB1"]
```

A figure with two subplots is created to compare gene expression patterns across different data representations using selected marker genes.

---

### Normalized Expression

```python
sc.pl.matrixplot(
    adata,
    groupby="louvain_cell_types",
    var_names=genes_of_interest,
    layer="counts_per_million",
    ax=ax1,
    show=False
)
```

Shows normalized gene expression across different cell clusters. It helps identify which genes are active in specific cell types.

### Raw Expression

```python
sc.pl.matrixplot(
    adata,
    groupby="louvain_cell_types",
    var_names=genes_of_interest,
    layer="raw",
    ax=ax2,
    show=False
)
ax2.set_title("raw counts")
plt.tight_layout()
```

Uses raw counts instead of normalized values, allowing comparison between raw and normalized gene expression.
### Result Output

<img width="990" height="240" alt="cp normal" src="https://github.com/user-attachments/assets/e0330bef-e583-4149-91e9-7a7765391543" />

---

## 🧾 Step 7: Cell Metadata

```python
adata.obs
```

The obs table contains metadata for each cell, including gene counts, total counts, mitochondrial percentage, and cluster assignments. This information is essential for quality control and downstream biological interpretation.
<img width="441" height="299" alt="sc2" src="https://github.com/user-attachments/assets/3a05999a-f88f-4ab6-8f84-22ef09424034" />


---

### Low Quality Cells

```python
adata.obs["is_low_quality"] = adata.obs["percent_mito"] > 0.03
```
A new column is created to identify low-quality cells based on mitochondrial gene expression. Cells with high mitochondrial content are often damaged or stressed
<img width="520" height="305" alt="sc6" src="https://github.com/user-attachments/assets/732b686d-3fe1-4093-a292-28b41c4edc67" />


---

##  Step 8: Gene Metadata

```python
adata.var
```

- Contains gene-level information  
- Includes gene names and statistics
<img width="313" height="303" alt="sc7" src="https://github.com/user-attachments/assets/7ad28dea-6825-4f83-9642-73497254fb54" />


---

## Step 9: Subsetting Data

```python
adata_small = adata[:5, ["LYZ", "FOS", "MALAT1"]]
```

```python
adata_high_quality = adata[~adata.obs["is_low_quality"], :]
```

A subset of the dataset is created by selecting specific cells and genes. AnnData ensures that all related structures are updated automatically, maintaining consistency across the dataset.
Low-quality cells are removed using boolean indexing. This improves data quality for downstream analysis.

---
<img width="543" height="294" alt="sc8" src="https://github.com/user-attachments/assets/911fae45-d53a-48de-abdb-606e9ee60bbd" />

##  Step 10: Embeddings

```python
adata.obsm
```
The obsm attribute contains low-dimensional embeddings such as PCA and UMAP, which help visualize complex high-dimensional gene expression data.

###  Visualization

```python
plt.figure(figsize=(12, 4))

# PCA
plt.subplot(1, 3, 1)
plt.scatter(
    x=adata.obsm["X_pca"][:, 0],  # PCA dim 1
    y=adata.obsm["X_pca"][:, 1],  # PCA dim 2
    c=adata.obs["louvain_cell_types"] == "B cells",  # B cell flag
    s=3,
    linewidth=0,
    cmap="coolwarm",
)
plt.title("PCA")
plt.axis("off")
plt.gca().set_aspect("equal")

# t-SNE
plt.subplot(1, 3, 2)
plt.scatter(
    x=adata.obsm["X_tsne"][:, 0],  # t-SNE dim 1
    y=adata.obsm["X_tsne"][:, 1],  # t-SNE dim 2
    c=adata.obs["louvain_cell_types"] == "B cells",  # B cell flag
    s=3,
    linewidth=0,
    cmap="coolwarm",
)
plt.title("t-SNE")
plt.axis("off")
plt.gca().set_aspect("equal")

# UMAP
plt.subplot(1, 3, 3)
plt.scatter(
    x=adata.obsm["X_umap"][:, 0],  # UMAP dim 1
    y=adata.obsm["X_umap"][:, 1],  # UMAP dim 2
    c=adata.obs["louvain_cell_types"] == "B cells",  # B cell flag
    s=3,
    linewidth=0,
    cmap="coolwarm",
)
plt.title("UMAP")
plt.axis("off")
plt.gca().set_aspect("equal")
)
```

This code creates a side-by-side comparison of three dimensionality reduction methods—PCA, t-SNE, and UMAP—using single-cell RNA-seq data. Each subplot shows the same cells projected into 2D space using different techniques stored in adata.obsm.

PCA captures the main global variation in the data, t-SNE focuses on local similarity between cells, and UMAP preserves both local and global structure, making clusters more interpretable. In all plots, each point represents a single cell, and cells are colored based on whether they belong to “B cells” using the louvain_cell_types annotation. Axis labels are removed for a cleaner look, and equal aspect ratio ensures proper shape representation of clusters.
<img width="945" height="350" alt="maps" src="https://github.com/user-attachments/assets/7213da22-2f48-46b3-84ab-7bc2b83412d5" />



---

##  Step 11: Distance Matrix

```python
plt.imshow(adata.obsp["distances_all"])
plt.colorbar(label="Euclidean distance in PCA space")
plt.show()
```
<img width="524" height="418" alt="mapy" src="https://github.com/user-attachments/assets/69c3cf29-a250-4be7-a088-01c0e695f1c3" />


- Displays pairwise cell distances  
- Helps visualize clustering structure
- Here, we cannot see much structure. This is because the order of the cells is random! Therefore, let’s re-order the cells such that they are sorted by cell type before plotting:
```python
reorder_by_celltype = np.argsort(adata.obs["louvain_cell_types"])
plt.imshow(adata[reorder_by_celltype, :].obsp["distances_all"])
plt.colorbar(label="Euclidean distance in PCA space")
plt.show()
```
<img width="524" height="418" alt="mappy2" src="https://github.com/user-attachments/assets/68c43075-054f-4dc4-ac98-bb00017b53e6" />

---

##  Step 12: Unstructured Data

```python
adata.uns
```

### Explanation
- Stores global analysis metadata  
- Includes PCA variance and clustering settings  

---

## Step 13: Views vs Copies

```python
adata_view = adata[:5, 5:10]
adata_copy = adata_view.copy()
```
<img width="677" height="131" alt="sc9" src="https://github.com/user-attachments/assets/153f8af9-d6b7-4283-ad51-662b24ad0680" />


- Views save memory (linked data)  
- Copies create independent dataset  

---

##  Final Summary

| Component | Description |
|----------|-------------|
| adata.X | Gene expression matrix |
| adata.obs | Cell metadata |
| adata.var | Gene metadata |
| adata.layers | Data versions (raw, normalized) |
| adata.obsm | Embeddings (UMAP, PCA) |
| adata.obsp | Distance matrices |
| adata.uns | Global metadata |

---

##  Conclusion

This workflow demonstrates how to:

- Load and explore scRNA-seq data  
- Normalize and visualize gene expression  
- Understand AnnData structure  
- Perform basic QC and filtering  

This forms the foundation for advanced bioinformatics analyses like clustering, differential expression, and pathway enrichment.

---

