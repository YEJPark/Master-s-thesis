import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind, wilcoxon, shapiro
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import os

# Load the data
file_path = 'D:/GSE239897/GSE239897_TissueCoreProject_fpkm.txt'
data = pd.read_csv(file_path, sep='\t')

# Define your hub genes
hub_genes = ['CD2', 'CD247', 'ITK', 'COL1A1', 'KLRK1', 'CD86', 'CDH5', 'CD28']

# Filter data for hub genes
filtered_data = data[data['gene_short_name'].isin(hub_genes)]

# Extract columns for COPD and Control conditions
copd_cols = [col for col in data.columns if 'COPD' in col]
ctrl_cols = [col for col in data.columns if 'CTRL' in col]

# Reshape the data for plotting
melted_data = pd.melt(filtered_data, id_vars=['gene_short_name'],
                      value_vars=copd_cols + ctrl_cols,
                      var_name='Sample', value_name='Expression')

# Add a new column for Condition based on the 'Sample' column
melted_data['Condition'] = melted_data['Sample'].apply(lambda x: 'COPD' if 'COPD' in x else 'Control')

# Create a directory for the plots
output_dir = 'hub_genes_plots'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

from scipy.stats import mannwhitneyu


# Perform plotting and statistical tests for each hub gene
for gene in hub_genes:
    gene_data = melted_data[melted_data['gene_short_name'] == gene]
    plt.figure(figsize=(6, 6))
    sns.boxplot(data=gene_data, x='Condition', y='Expression', width=0.5)
    plt.title(gene + ' Expression in COPD vs Control', fontsize=16, fontweight='bold')
    plt.xlabel('')
    plt.ylabel('Expression Level', fontsize=14, fontweight='bold')

    plt.xticks(fontsize=14, fontweight='bold')  # Update this line

    # Perform Shapiro-Wilk test for normality
    copd_expr = gene_data[gene_data['Condition'] == 'COPD']['Expression']
    ctrl_expr = gene_data[gene_data['Condition'] == 'Control']['Expression']
    copd_normal = shapiro(copd_expr).pvalue > 0.05
    ctrl_normal = shapiro(ctrl_expr).pvalue > 0.05

    # Choose the statistical test based on normality
    if copd_normal and ctrl_normal:
        # Both distributions are normal, use T-test
        t_stat, p_val = ttest_ind(copd_expr, ctrl_expr)
    else:
        # At least one distribution is not normal, use Mann-Whitney U test
        p_val = mannwhitneyu(copd_expr, ctrl_expr).pvalue

    # Annotate the p-value
    height = max(copd_expr.max(), ctrl_expr.max())
    stars = ''
    if p_val < 0.001:
        stars = '***'
    elif p_val < 0.01:
        stars = '**'
    elif p_val < 0.05:
        stars = '*'
    elif p_val >= 0.05:
        stars = 'ns'

    plt.text(0.5, height, stars, ha='center', va='bottom')
    plt.tight_layout()

    # Save the plot
    plt.savefig(os.path.join(output_dir, f'{gene}_expression_plot.png'))
    plt.close()
