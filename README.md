# Prion GWAS × CRISPR Reproducible Analysis

Reproduces Bayesian integration and enrichment analyses combining:
- GWAS gene-level p-values (MAGMA) from Mead et al., Nat Commun 2022 (CJD)
- CRISPR internalization screen gene list with log2-ratio and p-values

## Setup
```
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt

mkdir -p data outputs
# Put your inputs into ./data:
#  - MAGMA_genes_with_Entrez.xlsx  (columns: Gene, MAGMA Pvalue)
#  - CRISPR_hits.tsv               (columns: Gene Symbol, log2 Ratio, pValue[, directionality])
```

## Step 1: Bayesian combination
```
python bayes_combine.py   --gwas data/MAGMA_genes_with_Entrez.xlsx --gwas-gene-col Gene --gwas-p-col "MAGMA Pvalue"   --crystal data/CRISPR_hits.tsv --cr-gene-col "Gene Symbol" --cr-p-col pValue --cr-effect-col "log2 Ratio"   --W-gwas 0.2 --W-cr 0.2 --pi 0.01   --out outputs/combined_results.csv
```

## Step 2: Enrichment tests & plots
```
python enrich_and_plots.py   --gwas data/MAGMA_genes_with_Entrez.xlsx   --cr data/CRISPR_hits.tsv   --combined outputs/combined_results.csv   --outdir outputs
```

Outputs:
- outputs/combined_results.csv
- outputs/Supplemental_Table_S6_enrichment.csv
- outputs/Supplemental_Table_S7_families.csv
- outputs/Supplemental_Figure_S10.png
