#!/usr/bin/env python3
import argparse, numpy as np, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--cr", required=True)
    ap.add_argument("--combined", required=True)
    ap.add_argument("--outdir", default="outputs")
    ap.add_argument("--families", default="LRIG,SMAD,BMP,BMPR,SOX,GATA,FGF,FGFR,MAPK,TGF,TYRO,LRP")
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)

    # GWAS
    if args.gwas.lower().endswith((".xlsx", ".xls")):
        gwas = pd.read_excel(args.gwas)
    else:
        try:
            gwas = pd.read_csv(args.gwas)
        except Exception:
            gwas = pd.read_csv(args.gwas, sep="\t")
    gwas = gwas.rename(columns={"Gene":"gene", "MAGMA Pvalue": "p"})
    gwas["gene"] = gwas["gene"].astype(str).str.upper()
    gwas["p"] = pd.to_numeric(gwas["p"], errors="coerce")
    gwas = gwas[np.isfinite(gwas["p"]) & (gwas["p"] > 0)].copy()
    gwas["mlog10p"] = -np.log10(gwas["p"])

    # CRISPR
    try:
        cr = pd.read_csv(args.cr)
    except Exception:
        cr = pd.read_csv(args.cr, sep="\t")
    cr = cr.rename(columns={"Gene Symbol":"gene", "pValue":"p_cr", "log2 Ratio":"log2_ratio"})
    cr["gene"] = cr["gene"].astype(str).str.upper()
    cr["p_cr"] = pd.to_numeric(cr["p_cr"], errors="coerce")
    cr = cr[np.isfinite(cr["p_cr"]) & (cr["p_cr"] > 0)].copy()

    gwas_genes = set(gwas["gene"]); cr_hits = set(cr["gene"])
    overlap = sorted(list(gwas_genes & cr_hits))
    bg = sorted(list(gwas_genes - set(overlap)))
    gmap = gwas.set_index("gene")["p"].to_dict()
    mlog = gwas.set_index("gene")["mlog10p"].to_dict()

    hit_p = np.array([gmap[g] for g in overlap])
    hit_m = np.array([mlog[g] for g in overlap])
    bg_m = np.array([mlog[g] for g in bg])

    mw_u, mw_p = stats.mannwhitneyu(hit_m, bg_m, alternative="greater")
    ks_D, ks_p = stats.ks_2samp(hit_m, bg_m, alternative="greater")

    rng = np.random.default_rng(args.seed)
    pool = np.array(list(mlog.values())); n = len(hit_m)
    perm = [rng.choice(pool, size=n, replace=False).mean() for _ in range(args.permutations)]
    perm = np.array(perm)
    obs_mean = hit_m.mean()
    emp_p = (np.sum(perm >= obs_mean) + 1) / (args.permutations + 1)

    s6 = pd.DataFrame([{
        "n_gwas_genes": len(gwas_genes),
        "n_crispr_list": len(cr_hits),
        "n_overlap": len(overlap),
        "mean_mlog10p_hits": float(obs_mean),
        "mean_mlog10p_bg": float(bg_m.mean()),
        "MW_p_greater": mw_p,
        "KS_p_greater": ks_p,
        "perm_p_greater": emp_p
    }])
    s6_path = f"{args.outdir}/Supplemental_Table_S6_enrichment.csv"
    s6.to_csv(s6_path, index=False)

    fams = [x.strip() for x in args.families.split(",") if x.strip()]
    rows = []
    for fam in fams:
        fam_genes = [g for g in gwas_genes if fam in g]
        if not fam_genes: continue
        fam_m = np.array([mlog[g] for g in fam_genes])
        nonfam = [g for g in gwas_genes if fam not in g]
        nonfam_m = np.array([mlog[g] for g in nonfam])
        mw_u_f, mw_p_f = stats.mannwhitneyu(fam_m, nonfam_m, alternative="greater")
        perm_means = [rng.choice(nonfam_m, size=len(fam_m), replace=False).mean() for _ in range(args.permutations)]
        perm_means = np.array(perm_means)
        emp_f = (np.sum(perm_means >= fam_m.mean()) + 1) / (args.permutations + 1)
        rows.append({
            "family": fam,
            "n_genes_in_GWAS": len(fam_genes),
            "mean_-log10p": float(fam_m.mean()),
            "MannWhitney_p": mw_p_f,
            "perm_p": emp_f
        })
    s7 = pd.DataFrame(rows).sort_values("perm_p")
    s7_path = f"{args.outdir}/Supplemental_Table_S7_families.csv"
    s7.to_csv(s7_path, index=False)

    comb = pd.read_csv(args.combined)
    plt.figure(figsize=(6,5))
    sc = plt.scatter(comb["z_gwas"], comb["z_crystal"], c=comb["PPA"], s=30)
    plt.colorbar(sc, label="Posterior Probability of Association (PPA)")
    plt.xlabel("GWAS z-score"); plt.ylabel("CRISPR screen z-score")
    plt.title("Supplemental Figure S10. GWAS vs CRISPR z with PPA")
    plt.tight_layout()
    fig_path = f"{args.outdir}/Supplemental_Figure_S10.png"
    plt.savefig(fig_path, dpi=300)

    print("Wrote:", s6_path, s7_path, fig_path)

if __name__ == "__main__":
    main()
