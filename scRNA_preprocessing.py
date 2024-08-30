import pandas as pd
import scanpy as sc
import torch
import scvi
import matplotlib.pyplot as plt

### Step 1: Load and preprocess the data

# Load the raw count matrix
adata = sc.read_mtx("./data/GSE136831_RawCounts_Sparse.mtx").T # 45947: Number of genes / 312928: Number of cells.

# Assign Ensembl Gene IDs from GSE136831_AllCells.GeneIDs.txt to var_names
gene_names = pd.read_csv("./data/GSE136831_AllCells.GeneIDs.txt", sep="\t")["HGNC_EnsemblAlt_GeneID"].values
adata.var_names = gene_names

# Assign cell barcodes from GSE136831_AllCells.cellBarcodes.txt to obs_names
cell_barcodes = pd.read_csv("./data/GSE136831_AllCells.cellBarcodes.txt", header=None).iloc[:, 0].values
adata.obs_names = cell_barcodes

# Incorporate metadata from GSE136831_AllCells.Samples.CellType.MetadataTable_2.txt into the obs dataframe
metadata = pd.read_csv("./data/GSE136831_AllCells.Samples.CellType.MetadataTable_2.txt", sep="\t", index_col="CellBarcode_Identity")
adata.obs = metadata.loc[adata.obs_names]

# Ensure that the order of the metadata matches the order of the obs_names
assert all(adata.obs.index == adata.obs_names), "Mismatch between cell barcodes in adata and metadata."

# Save the preprocessed AnnData object as a checkpoint
adata.write('./result/processed_adata.h5ad')

### Step 2: Doublet Removal (optional but recommended)

# Set the desired GPU
torch.cuda.set_device(0)

# Setup and train the SCVI model
scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)
vae.train()

# Train SOLO for doublet detection
solo = scvi.external.SOLO.from_scvi_model(vae)
solo.train()

# Predict doublets using SOLO
df = solo.predict()
df['prediction'] = solo.predict(soft=False)
df['dif'] = df.doublet - df.singlet

# Identify and remove doublets
doublets = df[(df.prediction == 'doublet') & (df.dif > 1)]
adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
adata = adata[~adata.obs.doublet]

### Step 3: Quality Control

# Identify mitochondrial and ribosomal genes
adata.var['mt'] = adata.var.index.str.startswith('MT-')
ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)
adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)

# Calculate quality control metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

# Filter genes and cells based on QC metrics
sc.pp.filter_genes(adata, min_cells=3)
upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
adata = adata[adata.obs.n_genes_by_counts < upper_lim]
adata = adata[adata.obs.pct_counts_mt < 20]
adata = adata[adata.obs.pct_counts_ribo < 2]

### Step 4: Normalization

# Normalize to 10,000 UMI per cell and log transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # Save the original data

# Save the processed AnnData object
adata.write('./result/GSE136831_double_3.h5ad')

### Step 5: Dimensionality Reduction

# Identify highly variable genes and perform PCA
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt', 'pct_counts_ribo'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)

### Step 6: Clustering

# Compute the neighborhood graph and perform UMAP
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)
sc.pl.umap(adata)

# Perform Leiden clustering
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=['leiden'])

### Step 7: Marker Gene Identification

# Identify marker genes for each cluster using the Wilcoxon test
sc.tl.leiden(adata, resolution=1)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]

# Extract top 5 markers for each group
group_markers = {}
for group in markers['group'].unique():
    group_markers[group] = markers[markers['group'] == group].nlargest(5, 'logfoldchanges')['names'].tolist()

# Print top 5 markers for each group
for group, top_markers in group_markers.items():
    print(f"Group {group} markers: {', '.join(top_markers)}")

### Step 8: Cell Type Annotation

# Define cell types for each cluster
cell_type = {
    '0': 'Macrophage', '1': 'Macrophage_Alveolar', '2': 'ncMonocyte', '3': 'cMonocyte', '4': 'Club',
    '6': 'pDC', '7': 'Basal', '8': 'ATII', '9': 'NK', '11': 'cDC2',
    '12': 'B_Plasma', '13': 'T_Cytotoxic', '14': 'Ciliated', '15': 'VE_Venous', '16': 'ATI',
    '18': 'B', '19': 'T', '5': 'Goblet', '10': 'Mesothelial', '17': 'Pericyte',
    '20': 'VE_Peribronchial', '21': 'T_Regulatory', '22': 'Mast', '23': 'ILC_B', '24': 'DC_Langerhans'
}

# Map cell types to clusters and visualize with UMAP
adata.obs['cell type'] = adata.obs.leiden.map(cell_type)

# Final UMAP visualization with cell type annotations
fig, ax = plt.subplots(figsize=(12, 8))
sc.pl.umap(adata, color=['cell type'], frameon=False, ax=ax)
plt.tight_layout()
plt.savefig('./result/filename.png', dpi=300)
