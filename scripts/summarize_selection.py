#!/usr/bin/env python3
"""Summarize the selection track: Contrast-FEL sites and RELAX gene-level tests.

This is the only track that asks whether the SELECTIVE REGIME differs in diurnal
lineages, rather than whether the same residue appeared repeatedly. The two
dissociate routinely, so the convergence null constrains this result not at all.

Contrast-FEL, per site: is dN/dS different between diurnal branches and the rest?
It emits its own q-values within a gene, but the study tests 18 genes, so a
second correction across the pooled site list is applied here. Reporting only
HyPhy's within-gene q-values would understate the multiple-testing burden by
roughly the number of genes.

RELAX, per gene: is selection on diurnal branches RELAXED (K < 1, dN/dS pulled
toward 1) or INTENSIFIED (K > 1, pushed away from 1)? K is the interesting
quantity and its direction matters, so it is reported alongside the p-value
rather than collapsed into significance.

Inputs : results/selection/<gene>.contrastfel.json
         results/selection/<gene>.relax.json
Outputs: results/selection/contrastfel_sites.csv     (all sites, both q-values)
         results/selection/relax_summary.csv
         results/figures/selection_<...>.pdf
"""
import csv
import json
import os

import plotly.graph_objects as go
from plotly.subplots import make_subplots

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEL = os.path.join(PROJ, "results", "selection")
FIG = os.path.join(PROJ, "results", "figures")

GENES = ("CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
         "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

FDR = 0.05


def bh(pvals):
    """Benjamini-Hochberg q-values, returned in the input order."""
    n = len(pvals)
    if not n:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = n - rank + 1
        val = min(prev, pvals[i] * n / k)
        q[i] = val
        prev = val
    return q


def load_contrastfel(gene):
    p = os.path.join(SEL, f"{gene}.contrastfel.json")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return []
    d = json.load(open(p))
    mle = d.get("MLE", {})
    hdr = [h[0] for h in mle.get("headers", [])]
    rows = mle.get("content", {}).get("0", [])
    idx = {name: i for i, name in enumerate(hdr)}

    def col(r, *names):
        for n in names:
            if n in idx and idx[n] < len(r):
                return r[idx[n]]
        return None

    out = []
    for site, r in enumerate(rows, start=1):
        out.append(dict(
            gene=gene, site=site,
            alpha=col(r, "alpha"),
            beta_foreground=col(r, "beta (Foreground)"),
            beta_background=col(r, "beta (background)"),
            subs_foreground=col(r, "subs (Foreground)"),
            p_value=col(r, "P-value (overall)"),
            q_within_gene=col(r, "Q-value (overall)"),
            permutation_p=col(r, "Permutation p-value"),
        ))
    return out


def load_relax(gene):
    p = os.path.join(SEL, f"{gene}.relax.json")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return None
    d = json.load(open(p))
    tr = d.get("test results", {})
    k = tr.get("relaxation or intensification parameter")
    pv = tr.get("p-value")
    if k is None or pv is None:
        return None
    return dict(gene=gene, K=k, p_value=pv,
                direction=("intensified" if k > 1 else "relaxed"),
                LRT=tr.get("LRT"))


def main():
    sites, relax = [], []
    for g in GENES:
        sites.extend(load_contrastfel(g))
        r = load_relax(g)
        if r:
            relax.append(r)

    if not sites and not relax:
        raise SystemExit(f"No parsable selection output in {SEL}. Run scripts/12_selection_hyphy.sh")

    # ---- Contrast-FEL -------------------------------------------------------
    if sites:
        usable = [s for s in sites if isinstance(s["p_value"], (int, float))]
        qs = bh([s["p_value"] for s in usable])
        for s, q in zip(usable, qs):
            s["q_across_genes"] = q
        for s in sites:
            s.setdefault("q_across_genes", None)

        with open(os.path.join(SEL, "contrastfel_sites.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(sites[0].keys()))
            w.writeheader()
            w.writerows(sites)

        hits = [s for s in usable if s["q_across_genes"] is not None
                and s["q_across_genes"] <= FDR]
        within = [s for s in usable
                  if isinstance(s["q_within_gene"], (int, float))
                  and s["q_within_gene"] <= FDR]
        print(f"Contrast-FEL: {len(usable)} sites tested across "
              f"{len({s['gene'] for s in usable})} genes")
        print(f"  significant within gene   (q <= {FDR}): {len(within)}")
        print(f"  significant across genes  (q <= {FDR}): {len(hits)}")
        if hits:
            print(f"\n  {'gene':9s} {'site':>5s} {'dN/dS fg':>9s} {'dN/dS bg':>9s} {'q':>9s}")
            for s in sorted(hits, key=lambda x: x["q_across_genes"])[:25]:
                a = s["alpha"] or 0
                fg = (s["beta_foreground"] / a) if a else float("nan")
                bg = (s["beta_background"] / a) if a else float("nan")
                print(f"  {s['gene']:9s} {s['site']:5d} {fg:9.3f} {bg:9.3f} "
                      f"{s['q_across_genes']:9.2e}")

    # ---- RELAX --------------------------------------------------------------
    if relax:
        qs = bh([r["p_value"] for r in relax])
        for r, q in zip(relax, qs):
            r["q_value"] = q
        relax.sort(key=lambda r: r["p_value"])
        with open(os.path.join(SEL, "relax_summary.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(relax[0].keys()))
            w.writeheader()
            w.writerows(relax)
        print(f"\nRELAX: {len(relax)} genes")
        print(f"  {'gene':9s} {'K':>7s} {'direction':>13s} {'p':>9s} {'q':>9s}")
        for r in relax:
            print(f"  {r['gene']:9s} {r['K']:7.3f} {r['direction']:>13s} "
                  f"{r['p_value']:9.4f} {r['q_value']:9.4f}")
        sig = [r for r in relax if r["q_value"] <= FDR]
        print(f"  significant after correction (q <= {FDR}): "
              f"{', '.join(r['gene'] for r in sig) if sig else 'none'}")

    _figure(sites, relax)


def _figure(sites, relax):
    os.makedirs(FIG, exist_ok=True)
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.13,
                        subplot_titles=("RELAX: selection intensity K per gene",
                                        "Contrast-FEL: dN/dS, diurnal vs rest"))
    if relax:
        rs = sorted(relax, key=lambda r: r["K"])
        fig.add_trace(go.Bar(x=[r["K"] for r in rs], y=[r["gene"] for r in rs],
                             orientation="h", name="K",
                             marker=dict(color=["#c0392b" if r["K"] > 1 else "#2c6fbb"
                                                for r in rs])),
                      row=1, col=1)
        # K = 1 is the null: identical selection intensity on both branch sets.
        fig.add_vline(x=1.0, line=dict(dash="dash", width=1), row=1, col=1)

    if sites:
        pts = [s for s in sites if isinstance(s["alpha"], (int, float)) and s["alpha"]]
        fg = [s["beta_foreground"] / s["alpha"] for s in pts]
        bg = [s["beta_background"] / s["alpha"] for s in pts]
        sig = [isinstance(s.get("q_across_genes"), float) and s["q_across_genes"] <= FDR
               for s in pts]
        fig.add_trace(go.Scatter(
            x=bg, y=fg, mode="markers", name="site",
            marker=dict(size=[8 if s else 4 for s in sig],
                        color=["#c0392b" if s else "rgba(120,120,120,0.35)" for s in sig]),
            text=[f"{s['gene']} site {s['site']}" for s in pts],
            hovertemplate="%{text}<br>bg %{x:.3f}<br>fg %{y:.3f}<extra></extra>"),
            row=1, col=2)
        m = max([v for v in fg + bg if v == v] + [1.0])
        fig.add_trace(go.Scatter(x=[0, m], y=[0, m], mode="lines",
                                 line=dict(dash="dash", width=1, color="#888"),
                                 showlegend=False), row=1, col=2)

    fig.update_xaxes(title_text="K (below 1 relaxed, above 1 intensified)", row=1, col=1)
    fig.update_xaxes(title_text="dN/dS, background branches", row=1, col=2)
    fig.update_yaxes(title_text="dN/dS, diurnal branches", row=1, col=2)
    fig.update_layout(template="simple_white", width=1250, height=620, showlegend=False,
                      title="Selection on circadian genes in diurnal mammal lineages")
    out = os.path.join(FIG, "selection_summary.pdf")
    fig.write_image(out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
