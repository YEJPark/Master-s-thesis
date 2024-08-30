import gseapy as gp
import pandas as pd

# Load DEGs from CSV
filename_csv = "./result/results_DEGs_filtered_Control_vs_COPD.csv"
df = pd.read_csv(filename_csv)
degs = df["Gene"].tolist()

# GO Enrichment Analysis
enr = gp.enrichr(
    gene_list=degs,
    gene_sets='GO_Biological_Process_2021',
    outdir='Enrichment_GO',
    cutoff=0.02  # adjusted p-value
)

# Check if any significant results were found before saving or plotting
if not enr.results.empty:
    print(enr.results)
    # Save the results to CSV files
    enr.results.to_csv("./result/GO_enrichment_results_IPF.csv", index=False)
    # Visualization for GO
    gp.barplot(enr.res2d.head(5), title='GO Enrichment Analysis', ofname="GO_enrichment_barplot_COPD.png")
else:
    print("No significant GO enrichment terms found with p-value < 0.05")

# KEGG Enrichment Analysis
enr_kegg = gp.enrichr(
    gene_list=degs,
    gene_sets='KEGG_2021_Human',
    outdir='Enrichment_KEGG',
    cutoff=0.02  # adjusted p-value
)

if not enr_kegg.results.empty:
    print(enr_kegg.results)
    # Save the results to CSV files
    enr_kegg.results.to_csv("./result/KEGG_enrichment_results_IPF.csv", index=False)
    # Visualization for KEGG
    gp.barplot(enr_kegg.res2d.head(5), title='KEGG Enrichment Analysis', ofname="KEGG_enrichment_barplot_COPD.png")
else:
    print("No significant KEGG enrichment terms found with p-value < 0.05")
