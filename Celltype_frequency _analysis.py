import scanpy as sc
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the processed AnnData object
adata = sc.read('GSE136831_double_3.h5ad')

# Function to classify cells into broader categories
def find_class(cell_type):
    cell_class_dict = {
        'Endothelial': ['Lymphatic', 'VE_Venous', 'VE_Capillary_A', 'VE_Arterial', 'VE_Capillary_B', 'VE_Peribronchial'],
        'Epithelial': ['Basal', 'ATII', 'ATI', 'Ciliated', 'Club', 'Goblet', 'Mesothelial', 'Aberrant_Basaloid', 'Ionocyte', 'PNEC'],
        'Lymphoid': ['NK', 'B', 'T', 'B_Plasma', 'ILC_A', 'T_Cytotoxic', 'T_Regulatory', 'ILC_B'],
        'Myeloid': ['ncMonocyte', 'Macrophage_Alveolar', 'cMonocyte', 'Macrophage', 'cDC2', 'cDC1', 'DC_Langerhans', 'Mast', 'pDC', 'DC_Mature'],
        'Stromal': ['Fibroblast', 'Myofibroblast', 'SMC', 'Pericyte']
    }
    
    for cls, types in cell_class_dict.items():
        if cell_type in types:
            return cls
    return None

# Apply the classification function to assign broader cell classes
adata.obs['class'] = adata.obs['cell type'].apply(find_class)

# Counting cells
num_tot_cells = adata.obs.groupby(['Library_Identity']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.doublet))

cell_type_counts = adata.obs.groupby(['Library_Identity', 'Disease_Identity', 'Manuscript_Identity', 'CellType_Category']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis=1) > 0].reset_index()

cell_type_counts['total_cells'] = cell_type_counts.Library_Identity.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.doublet / cell_type_counts.total_cells

# Get unique CellType_Categories
unique_cell_type_categories = cell_type_counts['CellType_Category'].unique()

# Set style and increase the context scale for larger fonts
sns.set_style("white")
sns.set_context("talk")

# For each unique CellType_Category, create a separate plot
for category in unique_cell_type_categories:
    subset = cell_type_counts[cell_type_counts['CellType_Category'] == category]
    
    # Only keep the Manuscript_Identities that are relevant for this category
    relevant_manuscript_identities = subset['Manuscript_Identity'].unique()
    
    plt.figure(figsize=(15, 12))
    ax = sns.boxplot(data=subset, x='Manuscript_Identity', y='frequency', hue='Disease_Identity', dodge=True, order=relevant_manuscript_identities)
    
    # Customize x-ticks and legend
    plt.xticks(rotation=35, rotation_mode='anchor', ha='right', fontsize=25, fontweight='bold')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, title='Condition')
    ax.set_xlabel('Cell Type', fontsize=40, fontweight='bold')
    ax.set_title(f'{category}', fontsize=55, fontweight='bold')

    sns.despine()  # remove the top and right borders

    # Save and show the plot
    plt.tight_layout()
    plt.savefig(f"result_plot_{category}_4.png")
    plt.show()

# For whole data analysis (optional)
# Counting cells across the whole dataset
num_tot_cells = adata.obs.groupby(['Library_Identity']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.doublet))

cell_type_counts = adata.obs.groupby(['Library_Identity', 'Disease_Identity', 'Manuscript_Identity']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis=1) > 0].reset_index()

cell_type_counts['total_cells'] = cell_type_counts.Library_Identity.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.doublet / cell_type_counts.total_cells

# Plotting the whole dataset
plt.figure(figsize=(15, 8))
ax = sns.boxplot(data=cell_type_counts, x='Manuscript_Identity', y='frequency', hue='Disease_Identity', dodge=True)
plt.xticks(rotation=35, rotation_mode='anchor', ha='right')

# Customize x-ticks and legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels, title='Condition')
ax.set_xlabel('Cell Type')

plt.tight_layout()
plt.savefig("result_plot_whole.png")
plt.show()
