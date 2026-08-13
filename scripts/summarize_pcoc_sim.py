#!/usr/bin/env python3
"""Aggregate pcoc_sim benchmark output into a power/FPR calibration.

Usage: python scripts/summarize_pcoc_sim.py [SET]      (default: ER_gain)

Answers two questions the raw detection sweep cannot:

1. Could this design have detected convergence at all? pcoc_sim plants a known
   convergent shift on the real tree at the real transition branches, so the
   recovered sensitivity is the power of the actual analysis, per gene.
2. What posterior threshold should the detection sweep use? CLAUDE.md forbids a
   fixed 0.8 default, so the threshold is picked here as the lowest one whose
   false-positive rate stays within budget across every gene.

pcoc_sim reports three methods. Only PCOC is used for inference; PC and OC are
its two components and are carried through for diagnostics, since PC alone is
known to be uninformative on this data (its posterior is a constant 0.5).

Inputs : results/pcoc_sim/<SET>/<gene>/RUN_*/Tree_1/BenchmarkResults.tsv
Outputs: results/pcoc_sim/<SET>/summary_by_gene.csv
         results/pcoc_sim/<SET>/summary_by_distance.csv
         results/pcoc_sim/<SET>/threshold_choice.txt
         results/figures/pcoc_sim_power_<SET>.pdf
"""
import glob
import os
import sys

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SET = sys.argv[1] if len(sys.argv) > 1 else "ER_gain"
SIM_ROOT = os.path.join(PROJ, "results", "pcoc_sim", SET)
FIG_DIR = os.path.join(PROJ, "results", "figures")

# Maximum tolerated false-positive rate when choosing the threshold. Sites are
# tested in the tens of thousands across the 18 genes, so even 1 percent would
# swamp any real signal; 0.1 percent keeps the expected false count near zero.
MAX_FPR = 0.001

# Power may be traded away by this much to gain a stricter threshold.
POWER_TOL = 0.01


def load() -> pd.DataFrame:
    frames = []
    for f in sorted(glob.glob(os.path.join(SIM_ROOT, "*", "RUN_*", "Tree_1",
                                           "BenchmarkResults.tsv"))):
        gene = f.split(os.sep)[-4]
        d = pd.read_csv(f, sep="\t")
        d["gene"] = gene
        frames.append(d)
    if not frames:
        sys.exit(f"No BenchmarkResults.tsv under {SIM_ROOT}. Run run_pcoc_sim.sh first.")
    d = pd.concat(frames, ignore_index=True)
    # Specificity is over the noOneChange (null) sites, so FPR is its complement.
    d["FPR"] = 1.0 - d["Specificity"]
    return d


def main() -> None:
    d = load()
    pcoc = d[d["Method"] == "PCOC"].copy()

    genes = sorted(pcoc["gene"].unique())
    print(f"set={SET}  genes={len(genes)}  "
          f"couples/gene={pcoc.groupby('gene')['SimuCoupleID'].nunique().median():.0f}")

    # ---- per gene and threshold -------------------------------------------
    by_gene = (pcoc.groupby(["gene", "Threshold"])
               .agg(power=("Sensitivity", "mean"),
                    power_min=("Sensitivity", "min"),
                    fpr=("FPR", "mean"),
                    fpr_max=("FPR", "max"),
                    n_couples=("SimuCoupleID", "nunique"),
                    n_events=("NumberOfConvergentEvents", "first"))
               .reset_index())
    by_gene.to_csv(os.path.join(SIM_ROOT, "summary_by_gene.csv"), index=False)

    # ---- power as a function of profile divergence -------------------------
    # Power is expected to fall off as the ancestral and convergent profiles get
    # closer, so a single mean would hide the range where the method is blind.
    pcoc["dist_bin"] = pd.cut(pcoc["DistanceSimuCouple"],
                              bins=[0, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0])
    by_dist = (pcoc.groupby(["Threshold", "dist_bin"], observed=True)
               .agg(power=("Sensitivity", "mean"),
                    fpr=("FPR", "mean"),
                    n=("Sensitivity", "size"))
               .reset_index())
    by_dist.to_csv(os.path.join(SIM_ROOT, "summary_by_distance.csv"), index=False)

    # ---- threshold choice --------------------------------------------------
    # The threshold is shared across all 18 genes, so the weakest gene sets it:
    # aggregate on the WORST per-gene FPR and the WORST per-gene power.
    worst = by_gene.groupby("Threshold").agg(fpr_max=("fpr_max", "max"),
                                             power_min=("power", "min"),
                                             power_mean=("power", "mean")).reset_index()
    ok = worst[worst["fpr_max"] <= MAX_FPR]

    # Among thresholds meeting the FPR budget, take the HIGHEST one that keeps
    # power within POWER_TOL of the best available. Picking the lowest instead
    # would only make sense if there were a power/FPR tradeoff to buy back; when
    # FPR is already 0 across the range, a permissive threshold buys nothing and
    # costs robustness to the model misspecification the simulation cannot
    # capture (real alignments have indels and partial convergence, simulated
    # ones do not).
    if ok.empty:
        chosen = float(worst["Threshold"].max())
        note = "no threshold met the FPR budget; falling back to the strictest"
    else:
        best = ok["power_min"].max()
        keep = ok[ok["power_min"] >= best - POWER_TOL]
        chosen = float(keep["Threshold"].max())
        spread = ok["power_min"].max() - ok["power_min"].min()
        if spread <= POWER_TOL and ok["fpr_max"].max() <= MAX_FPR:
            note = ("power and FPR are saturated across the whole tested range, "
                    "so the calibration does not identify a threshold; the "
                    "strictest value is taken because it costs no power")
        else:
            note = "highest threshold retaining power within tolerance of the best"

    lines = [f"PCOC posterior threshold calibration, set {SET}",
             f"FPR budget: {MAX_FPR}",
             "",
             f"{'thr':>6s} {'worst_gene_FPR':>15s} {'min_gene_power':>15s} {'mean_power':>11s}"]
    for _, r in worst.iterrows():
        mark = "  <== chosen" if r["Threshold"] == chosen else ""
        lines.append(f"{r['Threshold']:6.2f} {r['fpr_max']:15.4f} "
                     f"{r['power_min']:15.3f} {r['power_mean']:11.3f}{mark}")
    lines += ["", f"CHOSEN_THRESHOLD={chosen}", f"basis: {note}", "",
              "Per-gene power at the chosen threshold:"]
    at = by_gene[by_gene["Threshold"] == chosen].sort_values("power")
    for _, r in at.iterrows():
        lines.append(f"  {r['gene']:10s} power={r['power']:.3f} "
                     f"(worst couple {r['power_min']:.3f})  FPR={r['fpr']:.4f}  "
                     f"events={int(r['n_events'])}")
    txt = "\n".join(lines)
    with open(os.path.join(SIM_ROOT, "threshold_choice.txt"), "w") as fh:
        fh.write(txt + "\n")
    print()
    print(txt)

    # ---- figure ------------------------------------------------------------
    os.makedirs(FIG_DIR, exist_ok=True)
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.12,
                        subplot_titles=("Power vs posterior threshold, per gene",
                                        "Power vs ancestral/convergent profile distance"))
    for g in genes:
        s = by_gene[by_gene["gene"] == g].sort_values("Threshold")
        fig.add_trace(go.Scatter(x=s["Threshold"], y=s["power"], name=g,
                                 mode="lines+markers", legendgroup=g),
                      row=1, col=1)
    for thr in sorted(by_dist["Threshold"].unique()):
        s = by_dist[by_dist["Threshold"] == thr]
        fig.add_trace(go.Scatter(x=[iv.mid for iv in s["dist_bin"]], y=s["power"],
                                 name=f"thr {thr:g}", mode="lines+markers",
                                 showlegend=True),
                      row=1, col=2)
    fig.update_xaxes(title_text="PCOC posterior threshold", row=1, col=1)
    fig.update_xaxes(title_text="euclidean distance between profiles", row=1, col=2)
    fig.update_yaxes(title_text="sensitivity (power)", range=[-0.02, 1.02], row=1, col=1)
    fig.update_yaxes(title_text="sensitivity (power)", range=[-0.02, 1.02], row=1, col=2)
    fig.update_layout(template="simple_white", width=1200, height=520,
                      title=f"PCOC power calibration on the real tree and scenario ({SET})")
    out = os.path.join(FIG_DIR, f"pcoc_sim_power_{SET}.pdf")
    fig.write_image(out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
