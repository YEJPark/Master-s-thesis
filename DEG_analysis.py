# Single-cell RNA-seq Differential Expression Analysis and Heatmap Generation

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your data
adata = sc.read('./result/GSE136831_double_3.h5ad')

# Remove cells with 'Multiplet' category in 'CellType_Category'
adata = adata[~(adata.obs['CellType_Category'] == 'Multiplet')]

def generate_heatmap_for_conditions(condition1, condition2):
    """
    Generates and saves a heatmap of differentially expressed genes (DEGs) between two conditions.

    Parameters:
    - condition1: First condition to compare (e.g., 'Control').
    - condition2: Second condition to compare (e.g., 'COPD').
    """

    # Filter the data for the specific conditions based on 'Disease_Identity'
    adata_filtered = adata[adata.obs['Disease_Identity'].isin([condition1, condition2])].copy()

    # Run PCA to reduce dimensions (required by Scanpy)
    sc.pp.pca(adata_filtered, svd_solver='arpack')

    # Differential Expression Analysis for the conditions
    categories = adata_filtered.obs['CellType_Category'].unique()
    sc.tl.rank_genes_groups(adata_filtered, 'CellType_Category', method='wilcoxon', corr_method='benjamini-hochberg')
    
    # Compute the dendrogram for better clustering visualization
    sc.tl.dendrogram(adata_filtered, groupby='CellType_Category')

    # Extracting top genes and their statistics
    result = adata_filtered.uns['rank_genes_groups']
    genes = result['names']
    logfoldchanges = result['logfoldchanges']
    pvals = result['pvals']
    pvals_adj = result['pvals_adj']

    data = []
    for category in categories:
        selected_genes = []
        for i in range(len(genes[category])):
            if pvals_adj[category][i] < 0.01 and abs(logfoldchanges[category][i]) > 2:
                selected_genes.append((genes[category][i], logfoldchanges[category][i], pvals[category][i], pvals_adj[category][i]))

        # Sort genes by absolute LogFC
        selected_genes = sorted(selected_genes, key=lambda x: abs(x[1]), reverse=True)[:20]

        for gene in selected_genes:
            data.append([category, gene[0], gene[1], gene[2], gene[3]])

    # Create a DataFrame for the selected genes
    df = pd.DataFrame(data, columns=['CellType_Category', 'Gene', 'LogFC', 'P-value', 'Adjusted P-value'])
    filename_csv = f"./result/results_DEGs_filtered_{condition1}_vs_{condition2}.csv"
    df.to_csv(filename_csv, index=False)

    # Prepare the colormap
    vmin = df['LogFC'].min()
    vmax = df['LogFC'].max()
    
    vmin = vmin if vmin < 0 else -np.abs(vmin)  # Ensure vmin is negative
    vmax = vmax if vmax > 0 else np.abs(vmax)   # Ensure vmax is positive

    midpoint = 1 - vmax / (vmax + abs(vmin))
    colors = [(0, "blue"), (midpoint, "white"), (1, "red")]
    cmap = plt.cm.colors.LinearSegmentedColormap.from_list("", colors)
    
    # Sort the DataFrame and pivot for heatmap
    sorted_df = df.sort_values(by=['CellType_Category'])
    heatmap_data = sorted_df.pivot(index='Gene', columns='CellType_Category', values='LogFC')

    # Set the style for seaborn
    sns.set_style("whitegrid", {'axes.grid': True, 'grid.color': 'lightgray'})
    
    # Plot the heatmap using seaborn
    sns.set(font_scale=1.7)
    plt.figure(figsize=(10, 50), facecolor='white')
    ax = sns.heatmap(heatmap_data, cmap=cmap, center=0, linecolor='lightgray', linewidths=0.5)
    ax.set_facecolor('white')
    
    colorbar = ax.collections[0].colorbar
    colorbar.set_label('logFC')
    
    plt.title(f"Heatmap of DEGs: {condition1} vs {condition2}")
    plt.savefig(f"./result/gene_heatmap_filtered_{condition1}_vs_{condition2}.png", dpi=300, bbox_inches='tight')

    # Generate and save Scanpy's heatmap with gene labels
    filename_png = f"./result/DEGs_heatmap_filtered_{condition1}_vs_{condition2}.png"
    sc.pl.rank_genes_groups_heatmap(adata_filtered, n_genes=20, groups=categories, use_raw=False, swap_axes=True, vmin=vmin, vmax=vmax, cmap='bwr', figsize=(18, 10), show_gene_labels=True, show=False, save=filename_png)
    
# Generate heatmap for Control vs COPD
generate_heatmap_for_conditions('Control', 'COPD')

# Generate heatmap for Control vs IPF
generate_heatmap_for_conditions('Control', 'IPF')
