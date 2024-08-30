import scanpy as sc
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read('./result/GSE136831_double_3.h5ad')

 
# Counting cells
adata.obs.Library_Identity.unique().tolist()

# Disease_Identity
num_tot_cells = adata.obs.groupby(['Library_Identity']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.doublet))
num_tot_cells

cell_type_counts = adata.obs.groupby(['Library_Identity', 'Disease_Identity', 'Manuscript_Identity', 'CellType_Category']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis=1) > 0].reset_index()

cell_type_counts['total_cells'] = cell_type_counts.Library_Identity.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.doublet / cell_type_counts.total_cells

# Get unique CellType_Categories
unique_cell_type_categories = cell_type_counts['CellType_Category'].unique()

# Set style and increase the context scale for larger fonts
sns.set_style("white")
sns.set_context("talk")  # this context mode increases font size

# For each unique CellType_Category, create a separate plot
for category in unique_cell_type_categories:
    subset = cell_type_counts[cell_type_counts['CellType_Category'] == category]
    
    # Only keep the Manuscript_Identities that are relevant for this category
    relevant_manuscript_identities = subset['Manuscript_Identity'].unique()
    
    plt.figure(figsize=(15, 12))
    ax = sns.boxplot(data=subset, x='Manuscript_Identity', y='frequency', hue='Disease_Identity', dodge=True, order=relevant_manuscript_identities)
    
    #plt.xticks(rotation=35, rotation_mode='anchor', ha='right')
    plt.xticks(rotation=35, rotation_mode='anchor', ha='right', fontsize=25, fontweight='bold')


    # Change the legend title and x-axis label
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, title='Condition')
    ax.set_xlabel('Cell Type', fontsize=40, fontweight='bold')  
    #ax.set_title(f'Cell Class: {category}', fontsize=0, fontweight='bold')
    ax.set_title(f'{category}', fontsize=55, fontweight='bold')

    sns.despine()  # remove the top and right borders

    plt.tight_layout()
    plt.savefig(f"./result/result_plot_{category}_4.png")
    plt.show()
