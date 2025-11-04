#!/usr/bin/env bash
set -euo pipefail
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
mkdir -p data outputs
echo "Place MAGMA_genes_with_Entrez.xlsx and CRISPR_hits.tsv into ./data, then run the two commands in README.md."
