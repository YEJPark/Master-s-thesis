import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind, mannwhitneyu, shapiro
import matplotlib.pyplot as plt
import numpy as np
import os

# Load the data
file_path = 'D:/RNAseq.SRP095361.FPKM.txt'  # Update the file path to your data file's location
data = pd.read_csv(file_path, sep='\t')  # Replace '\t' with the actual delimiter if different

# Define your hub genes
hub_genes = ['IL1B', 'PDGFRB', 'CDH5', 'VWF', 'COL1A1', 'THY1']  # Update this list with your actual hub genes

# Filter data for hub genes
filtered_data = data[data['genename'].isin(hub_genes)]

# Reshape the data for plotting
melted_data = pd.melt(filtered_data, id_vars=['genename', 'group'], value_vars=['FPKM'], var_name='Sample', value_name='Expression')

# Add a new column for Condition based on the 'group' column
melted_data['Condition'] = melted_data['group']

# Create a directory for the plots
output_dir = 'hub_genes_plots/IPF_2'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Perform plotting and statistical tests for each hub gene
for gene in hub_genes:
    gene_data = melted_data[melted_data['genename'] == gene]
    plt.figure(figsize=(6, 6))
    sns.boxplot(data=gene_data, x='Condition', y='Expression', width=0.5)
    plt.title(f'{gene} Expression in IPF vs Control', fontsize=16, fontweight='bold')
    plt.xlabel('')
    plt.ylabel('Expression Level', fontsize=14, fontweight='bold')
    plt.xticks(fontsize=14, fontweight='bold')  # Updated line for x-axis tick labels

    # Perform Shapiro-Wilk test for normality
    control_expr = gene_data[gene_data['Condition'] == 'Control']['Expression']
    ipf_expr = gene_data[gene_data['Condition'] != 'Control']['Expression']
    control_normal = shapiro(control_expr).pvalue > 0.05
    ipf_normal = shapiro(ipf_expr).pvalue > 0.05

    # Choose the statistical test based on normality
    if control_normal and ipf_normal:
        # Both distributions are normal, use T-test
        t_stat, p_val = ttest_ind(ipf_expr, control_expr)
    else:
        # At least one distribution is not normal, use Mann-Whitney U test
        p_val = mannwhitneyu(ipf_expr, control_expr).pvalue

    # Annotate the p-value
    height = max(ipf_expr.max(), control_expr.max())
    stars = 'ns' if p_val >= 0.05 else ('*' if p_val < 0.05 else ('**' if p_val < 0.01 else '***'))
    plt.text(0.5, height, stars, ha='center', va='bottom')

    plt.tight_layout()

    # Save the plot
    plot_filename = f'{gene}_expression_plot.png'
    plt.savefig(os.path.join(output_dir, plot_filename))
    plt.close()

print("Boxplots have been generated for the hub genes.")
