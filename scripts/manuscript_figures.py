#!/usr/bin/env python3
"""Manuscript figures 1 to 3. Plotly, vector PDF, per project convention.

Palette is the validated categorical default (slots 1 to 4). Validated with the
dataviz validator rather than by eye:
  4 slots, light: lightness band PASS, chroma PASS, CVD adjacent dE 9.1 PASS,
  normal-vision 22.9 PASS, contrast WARN on aqua and yellow.
The contrast warning obligates relief, so every bar in Figure 2 carries a visible
direct label. Figures 1 and 3 use slots 1 and 2 only, which pass all-pairs
including contrast.

Usage: python scripts/manuscript_figures.py
"""
import csv
import os
from collections import defaultdict

import plotly.graph_objects as go
import plotly.io as pio
from ete3 import Tree

# kaleido renders a "Loading [MathJax]" placeholder into the output unless MathJax
# is disabled. None of these figures use TeX, so switch it off.
try:
    pio.kaleido.scope.mathjax = None
except Exception:
    pass

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(PROJ, "results", "figures")
# results/ is gitignored, so a figure written only there is invisible to anyone
# reading the repository. Every figure is mirrored into shared_results/figures,
# which is tracked, so the two cannot drift apart the way they did before.
SHARED_FIG = os.path.join(PROJ, "shared_results", "figures")


def write_fig(fig, name):
    for d in (FIG, SHARED_FIG):
        os.makedirs(d, exist_ok=True)
        out = os.path.join(d, name)
        fig.write_image(out)
        print("wrote", out)
os.makedirs(FIG, exist_ok=True)

S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8880"
SURFACE = "#fcfcfb"

BASE = dict(
    plot_bgcolor=SURFACE, paper_bgcolor=SURFACE,
    font=dict(family="Helvetica, Arial, sans-serif", size=12, color=INK),
    margin=dict(l=70, r=30, t=54, b=54),
)


def axis(**kw):
    d = dict(showgrid=False, zeroline=False, linecolor=MUTED, linewidth=1,
             ticks="outside", tickcolor=MUTED, tickfont=dict(color=INK2, size=11))
    d.update(kw)
    return d


# ---------------------------------------------------------------- Figure 2
def figure2():
    """Positive control recovery by method and difficulty level."""
    p = os.path.join(PROJ, "results", "spikein", "spikein_scorecard.csv")
    if not os.path.exists(p):
        print("fig2: no scorecard, skipped"); return
    rec = defaultdict(dict)
    for r in csv.DictReader(open(p)):
        rec[r["method"]][r["level"]] = (int(r["recovered"]), int(r["planted"]))

    levels = ["all", "most", "few"]
    labels = ["All 10 lineages", "7 of 10", "3 of 10"]
    order = ["PCOC", "TDG09", "parsimony (flagged)", "parsimony (phenotype-specific)"]
    short = ["PCOC", "TDG09", "Parsimony", "Parsimony, phenotype-specific"]
    colors = [S1, S2, S3, S4]

    # False positive counts belong in the legend. Without them the chart reads as
    # "TDG09 beats PCOC at 3 of 10", when TDG09 reaches that recall by flagging a
    # quarter of all testable sites.
    fp = {}
    for r in csv.DictReader(open(p)):
        fp[r["method"]] = int(r["false_positives"])

    fig = go.Figure()
    for m, name, c in zip(order, short, colors):
        if m not in rec:
            continue
        name = f"{name} ({fp.get(m, 0)} false pos.)"
        ys, txt = [], []
        for lv in levels:
            got, want = rec[m].get(lv, (0, 0))
            ys.append(100 * got / want if want else 0)
            txt.append(f"{got}/{want}")
        fig.add_bar(x=labels, y=ys, name=name, marker_color=c, text=txt,
                    textposition="outside", textfont=dict(size=10, color=INK2),
                    cliponaxis=False,
                    marker_line=dict(width=2, color=SURFACE))

    fig.update_layout(
        **BASE, barmode="group", bargap=0.30, bargroupgap=0.04,
        title=dict(text="Recovery of planted convergent sites",
                   x=0, xanchor="left", font=dict(size=15)),
        xaxis=axis(title="Convergent lineages planted"),
        yaxis=axis(title="Planted sites recovered (%)", range=[0, 112],
                   showgrid=True, gridcolor="#eceae4"),
        legend=dict(orientation="h", y=-0.20, x=0, font=dict(size=10)),
        width=760, height=440,
    )
    write_fig(fig, "fig2_positive_control.pdf")


# ---------------------------------------------------------------- Figure 3
def figure3():
    """Observed against null same-residue rate, per gene."""
    p = os.path.join(PROJ, "shared_results", "pcoc_sim_calibration",
                     "observed_convergent_substitutions.csv")
    if not os.path.exists(p):
        print("fig3: no parsimony table, skipped"); return
    rows = []
    for r in csv.DictReader(open(p)):
        try:
            obs = float(r["obs_same_residue_rate"]); nul = float(r["null_same_residue_rate"])
        except (ValueError, KeyError):
            continue
        n = int(r["sites_2plus_events_changed"])
        if n == 0:
            continue
        rows.append((r["gene"], obs, nul, n))
    rows.sort(key=lambda t: t[1] - t[2])

    genes = [t[0] for t in rows]
    fig = go.Figure()
    for g, obs, nul, n in rows:
        fig.add_trace(go.Scatter(x=[nul, obs], y=[g, g], mode="lines",
                                 line=dict(color="#d9d6ce", width=2),
                                 showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter(
        x=[t[2] for t in rows], y=genes, mode="markers", name="Null expectation",
        marker=dict(color=S2, size=9, line=dict(width=2, color=SURFACE))))
    fig.add_trace(go.Scatter(
        x=[t[1] for t in rows], y=genes, mode="markers", name="Observed",
        marker=dict(color=S1, size=9, line=dict(width=2, color=SURFACE))))

    fig.update_layout(
        **BASE,
        title=dict(text="Same-residue rate among sites changing in 2 or more diurnal lineages",
                   x=0, xanchor="left", font=dict(size=15)),
        xaxis=axis(title="Proportion of changed sites reaching the same residue",
                   range=[0, 0.8], showgrid=True, gridcolor="#eceae4"),
        yaxis=axis(title="", autorange="reversed"),
        legend=dict(orientation="h", y=-0.13, x=0, font=dict(size=10)),
        width=720, height=560,
    )
    write_fig(fig, "fig3_parsimony_rates.pdf")


# ---------------------------------------------------------------- Figure 1
def figure1():
    """Species tree with diel states and reconstructed transitions."""
    sp = os.path.join(PROJ, "data", "species_tree.nwk")
    diel_p = os.path.join(PROJ, "data", "diel_activity.csv")
    tb = os.path.join(PROJ, "shared_results", "scenario", "transition_branches_ER.csv")
    if not all(os.path.exists(x) for x in (sp, diel_p, tb)):
        print("fig1: missing inputs, skipped"); return

    t = Tree(sp, format=1)
    diel = {r["species"]: r["activity"] for r in csv.DictReader(open(diel_p))}

    # ape numbers tips 1..Ntip in tree.tip.label order, internals Ntip+1 upward in
    # preorder. Rebuild that mapping so the R transition table can be joined.
    leaves = t.get_leaves()
    ape = {}
    for i, lf in enumerate(leaves, start=1):
        ape[i] = lf
    n = len(leaves) + 1
    for node in t.traverse("preorder"):
        if not node.is_leaf():
            ape[n] = node; n += 1

    trans = {}
    for r in csv.DictReader(open(tb)):
        child = int(r["child"]); trans[child] = r["direction"].strip('"')

    # layout: x = cumulative distance from root, y = leaf order
    ypos, ynext = {}, [0.0]
    for node in t.traverse("postorder"):
        if node.is_leaf():
            ypos[node] = ynext[0]; ynext[0] += 1.0
        else:
            kids = [ypos[c] for c in node.children]
            ypos[node] = sum(kids) / len(kids)
    xpos = {t: 0.0}
    for node in t.traverse("preorder"):
        if node.up is not None:
            xpos[node] = xpos[node.up] + (node.dist or 0.0)
    # x is distance from the root, so tips sit at the RIGHT (the present) and the
    # root at the left. Axis ticks are relabelled as time before present rather
    # than reversing the axis, which would put the present on the left.
    depth = max(xpos[lf] for lf in leaves)

    def seg(nodes, color, width):
        xs, ys = [], []
        for nd in nodes:
            if nd.up is None:
                continue
            xs += [xpos[nd.up], xpos[nd], None]
            ys += [ypos[nd], ypos[nd], None]
            xs += [xpos[nd.up], xpos[nd.up], None]
            ys += [ypos[nd.up], ypos[nd], None]
        return go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=color, width=width),
                          hoverinfo="skip", showlegend=False)

    inv = {id(v): k for k, v in ape.items()}
    tr_nodes = [nd for nd in t.traverse() if inv.get(id(nd)) in trans]
    plain = [nd for nd in t.traverse() if nd not in tr_nodes]

    fig = go.Figure()
    fig.add_trace(seg(plain, "#c9c6bf", 1.2))
    gains = [nd for nd in tr_nodes if trans[inv[id(nd)]] == "gain"]
    revs = [nd for nd in tr_nodes if trans[inv[id(nd)]] == "reversal"]
    fig.add_trace(seg(gains, S1, 3.0))
    fig.add_trace(seg(revs, S2, 3.0))

    for state, color, nm in (("diurnal", S1, "Diurnal"), ("nocturnal", INK2, "Nocturnal")):
        sel = [lf for lf in leaves if diel.get(lf.name) == state]
        fig.add_trace(go.Scatter(
            x=[xpos[lf] for lf in sel], y=[ypos[lf] for lf in sel],
            mode="markers", name=nm,
            marker=dict(color=color, size=8, line=dict(width=1.5, color=SURFACE)),
            text=[lf.name.replace("_", " ") for lf in sel], hoverinfo="text"))

    # legend entries for the branch classes
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines", name="Gain of diurnality",
                             line=dict(color=S1, width=3)))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines", name="Reversal to nocturnality",
                             line=dict(color=S2, width=3)))

    step = 25
    ticks = [v for v in range(0, int(depth) + step, step) if v <= depth]
    fig.update_layout(
        **BASE,
        title=dict(text="Diel activity and reconstructed transitions across 60 mammals",
                   x=0, xanchor="left", font=dict(size=15)),
        xaxis=axis(title="Time before present (My)",
                   tickvals=[depth - v for v in ticks],
                   ticktext=[str(v) for v in ticks],
                   range=[-2, depth * 1.42]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False,
                   range=[len(leaves) + 0.5, -1.5]),
        legend=dict(orientation="h", y=-0.055, x=0, font=dict(size=10)),
        width=700, height=1040,
    )
    # Species labels. 60 tips over 1040px leaves room for 7pt italic binomials.
    for lf in leaves:
        fig.add_annotation(
            x=xpos[lf], y=ypos[lf], text=f"<i>{lf.name.replace('_', ' ')}</i>",
            xanchor="left", yanchor="middle", xshift=7, showarrow=False,
            font=dict(size=7, color=INK if diel.get(lf.name) == "diurnal" else INK2))
    write_fig(fig, "fig1_tree.pdf")


if __name__ == "__main__":
    figure2()
    figure3()
    figure1()
