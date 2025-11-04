#!/usr/bin/env python3
import argparse, math
import numpy as np, pandas as pd
from scipy.stats import norm

def infer_z_from_p(p_two_sided):
    p = np.clip(np.asarray(p_two_sided, dtype=float), 1e-300, 1.0)
    return norm.isf(p / 2.0)

def wakefield_log_abf(z, V, W):
    z = np.asarray(z, dtype=float)
    V = np.asarray(V, dtype=float)
    W = float(W)
    V = np.where((~np.isfinite(V)) | (V <= 0), np.nan, V)
    logBF = np.full_like(z, np.nan, dtype=float)
    valid = np.isfinite(z) & np.isfinite(V) & (V > 0) & np.isfinite(W) & (W > 0)
    if np.any(valid):
        logBF[valid] = -0.5 * np.log1p(W / V[valid]) + 0.5 * (z[valid] ** 2) * ( W / (V[valid] + W) )
    return logBF

def posterior_from_logbf(logBF, pi):
    logit_pi = math.log(pi) - math.log(1 - pi)
    log_odds = logBF + logit_pi
    return 1.0 / (1.0 + np.exp(-np.clip(log_odds, -50, 50)))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--gwas-gene-col", default="Gene")
    ap.add_argument("--gwas-p-col", default="MAGMA Pvalue")
    ap.add_argument("--crystal", required=True)
    ap.add_argument("--cr-gene-col", default="Gene Symbol")
    ap.add_argument("--cr-p-col", default="pValue")
    ap.add_argument("--cr-effect-col", default="log2 Ratio")
    ap.add_argument("--W-gwas", type=float, default=0.2)
    ap.add_argument("--W-cr", type=float, default=0.2)
    ap.add_argument("--pi", type=float, default=0.01)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # Load GWAS
    if args.gwas.lower().endswith((".xlsx", ".xls")):
        gwas = pd.read_excel(args.gwas)
    else:
        try:
            gwas = pd.read_csv(args.gwas)
        except Exception:
            gwas = pd.read_csv(args.gwas, sep="\t")
    gwas = gwas.rename(columns={args.gwas_gene_col: "gene", args.gwas_p_col: "p_gwas"})
    gwas["gene"] = gwas["gene"].astype(str).str.upper()
    gwas["p_gwas"] = pd.to_numeric(gwas["p_gwas"], errors="coerce")
    gwas = gwas[np.isfinite(gwas["p_gwas"]) & (gwas["p_gwas"] > 0)].copy()
    gwas["z_gwas"] = infer_z_from_p(gwas["p_gwas"])
    gwas["V_gwas"] = 1.0

    # Load CRISPR
    try:
        cr = pd.read_csv(args.crystal)
    except Exception:
        cr = pd.read_csv(args.crystal, sep="\t")
    cr = cr.rename(columns={args.cr_gene_col: "gene", args.cr_p_col: "p_cr", args.cr_effect_col: "effect"})
    cr["gene"] = cr["gene"].astype(str).str.upper()
    cr["p_cr"] = pd.to_numeric(cr["p_cr"], errors="coerce")
    cr["effect"] = pd.to_numeric(cr["effect"], errors="coerce")
    cr = cr[np.isfinite(cr["p_cr"]) & (cr["p_cr"] > 0)].copy()
    sign = np.sign(cr["effect"].values)
    sign[sign == 0] = 1
    cr["z_cr"] = infer_z_from_p(cr["p_cr"].values) * sign
    cr["V_cr"] = 1.0

    df = pd.merge(gwas[["gene","z_gwas","V_gwas"]], cr[["gene","z_cr","V_cr"]], on="gene", how="inner")
    if df.empty:
        raise SystemExit("No overlapping genes between GWAS and CRISPR tables.")

    logBF_gwas = wakefield_log_abf(df["z_gwas"].values, df["V_gwas"].values, args.W_gwas)
    logBF_cr = wakefield_log_abf(df["z_cr"].values, df["V_cr"].values, args.W_cr)
    logBF_combined = logBF_gwas + logBF_cr
    PPA = posterior_from_logbf(logBF_combined, args.pi)

    out = df.assign(
        logBF_gwas=logBF_gwas,
        logBF_crystal=logBF_cr,
        logBF_combined=logBF_combined,
        BF_combined=np.exp(np.clip(logBF_combined, -50, 50)),
        PPA=PPA
    ).sort_values("PPA", ascending=False).reset_index(drop=True)
    out.to_csv(args.out, index=False)

if __name__ == "__main__":
    main()
