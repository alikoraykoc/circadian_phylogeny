#!/usr/bin/env python3
"""
13_consensus.py

Merge per-site evidence across methods into one table, translate positions to
human reference residue numbering, and flag high-confidence sites (supported by
>= 2 independent methods). Figures are written as vector PDF with Plotly.

Inputs (per gene, as available):
  - PCOC posteriors         results/pcoc/<gene>/...
  - TDG09 LRTs              results/tdg09/<gene>.tdg09.out
  - Contrast-FEL JSON       results/selection/<gene>.contrastfel.json
  - trimmed->human map      results/trim/<gene>.colnumbering.txt + REF row
"""
import os
import json
import pandas as pd
import plotly.graph_objects as go

PROJ = os.environ.get("PROJ", os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
RES = os.path.join(PROJ, "results")
REF_SPECIES = os.environ.get("REF_SPECIES", "Homo_sapiens")
GENES = os.environ.get("GENES", "CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

PCOC_THRESHOLD = 0.90          # EDIT: use the value calibrated by pcoc_sim (step 09)
MIN_METHODS = 2                # a site is "high confidence" if backed by >= this many

def load_pcoc(gene):
    # EDIT: parse the PCOC per-site posterior table for this gene.
    # return dict {trimmed_col: posterior_PCOC}
    return {}

def load_tdg09(gene):
    # EDIT: parse LRT / p-value per site from the TDG09 output.
    return {}

def load_contrastfel(gene):
    p = os.path.join(RES, "selection", f"{gene}.contrastfel.json")
    if not os.path.exists(p):
        return {}
    with open(p) as fh:
        _ = json.load(fh)
    # EDIT: pull per-site q-values for differential dN/dS from the JSON.
    return {}

def trimmed_to_human(gene):
    # EDIT: use results/trim/<gene>.colnumbering.txt (retained original columns)
    # plus the REF_SPECIES row of the alignment to map trimmed column -> human residue.
    return {}

def main():
    rows = []
    for g in GENES:
        pcoc = load_pcoc(g)
        tdg  = load_tdg09(g)
        cfel = load_contrastfel(g)
        h    = trimmed_to_human(g)
        cols = set(pcoc) | set(tdg) | set(cfel)
        for c in cols:
            hits = {
                "pcoc": pcoc.get(c, 0) >= PCOC_THRESHOLD,
                "tdg09": bool(tdg.get(c, False)),
                "contrastfel": bool(cfel.get(c, False)),
            }
            n = sum(hits.values())
            rows.append({
                "gene": g,
                "trimmed_col": c,
                "human_residue": h.get(c),
                **hits,
                "n_methods": n,
                "high_confidence": n >= MIN_METHODS,
            })

    df = pd.DataFrame(rows).sort_values(["high_confidence", "n_methods"], ascending=False)
    out_csv = os.path.join(RES, "consensus", "consensus_sites.csv")
    df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv} ({df['high_confidence'].sum()} high-confidence sites).")

    # ---- summary figure: high-confidence sites per gene (vector PDF) ----
    if not df.empty:
        agg = df[df.high_confidence].groupby("gene").size().reindex(GENES, fill_value=0)
        fig = go.Figure(go.Bar(x=list(agg.index), y=list(agg.values)))
        fig.update_layout(title="High-confidence convergent sites per gene",
                          xaxis_title="Gene", yaxis_title="Sites (>= 2 methods)",
                          template="simple_white")
        fig.write_image(os.path.join(RES, "figures", "consensus_by_gene.pdf"))
        print("Wrote results/figures/consensus_by_gene.pdf")

if __name__ == "__main__":
    main()
