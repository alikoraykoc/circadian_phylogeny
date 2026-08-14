#!/usr/bin/env python3
"""
13_consensus.py

Merge per-site evidence across methods onto shared coordinates and flag sites
supported by more than one method. Figures are vector PDF via Plotly.

The consensus rule is the point of the whole design: a site counts only if at
least MIN_METHODS independent methods flag it. That matters here more than
usual, because the methods disagree wildly. PCOC and Contrast-FEL each report
ZERO sites across the dataset, while TDG09 reports 885 of 3820 tested sites at
FDR <= 0.05, which is 23 percent. A quarter of all variable sites showing
phenotype-associated shifts is not a plausible biological result; TDG09 asks
whether amino acid FREQUENCIES differ between the two branch sets, and diurnal
species are phylogenetically clustered, so any site whose residues track taxonomy
differs without any convergence being involved. The consensus requirement is what
stops that from reaching the results.

Coordinates
-----------
Each method reports positions in a different frame:
  PCOC and TDG09  trimmed protein alignment column
  Contrast-FEL    codon index in the CLEANED codon alignment
Both are mapped to the untrimmed protein column, then to the REF_SPECIES residue
number, so a claim can be checked against a real sequence:
  trimmed col -> untrimmed col   via results/trim/<gene>.colnumbering.txt
  hyphy site  -> untrimmed col   via results/codon/<gene>.codonmap.csv
  untrimmed col -> ref residue   by counting non-gap characters in the REF row

A parsing trap worth keeping
----------------------------
TDG09 writes per-site `# Result{...}` comment lines whose `lrt` field is ALWAYS
0.0, including at sites whose two model likelihoods differ by 9.7 log units.
Those are a progress log written before the test is evaluated. The real numbers
are in the `FullResults:` block. Parsing the comments yields a clean all-zero
null that looks entirely plausible.
"""
import csv
import glob
import json
import os

import plotly.graph_objects as go

PROJ = os.environ.get("PROJ", os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
RES = os.path.join(PROJ, "results")
REF_SPECIES = os.environ.get("REF_SPECIES", "Homo_sapiens")
GENES = os.environ.get(
    "GENES",
    "CLOCK NPAS2 ARNTL PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 "
    "RORA RORB RORC CSNK1D CSNK1E FBXL3 BHLHE40 BHLHE41").split()

# From results/pcoc_sim/ER_gain/threshold_choice.txt. The calibration does not
# identify a threshold (power and FPR are saturated from 0.70 to 0.99), so the
# strictest value is used because it costs no power. Moot in practice: the
# highest posterior anywhere in the real data is 0.098.
PCOC_THRESHOLD = 0.99
FDR_THRESHOLD = 0.05
MIN_METHODS = 2


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def load_pcoc(gene):
    """{trimmed_col (1-based): PCOC posterior}. Blank cells underflowed to 0."""
    f = glob.glob(os.path.join(RES, "pcoc", "ER_gain", gene, "RUN_*",
                               f"{gene}.trim.results.tsv"))
    if not f:
        return {}
    out = {}
    with open(f[0]) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        col = next((c for c in rd.fieldnames if c.startswith("PCOC")), None)
        site_col = next((c for c in rd.fieldnames if c.lower().startswith("site")), None)
        for i, r in enumerate(rd, start=1):
            v = (r.get(col) or "").strip()
            site = int(r[site_col]) if site_col and r.get(site_col) else i
            out[site] = float(v) if v else 0.0
    return out


def load_tdg09(gene):
    """{trimmed_col: FDR} from the FullResults block, NOT the # Result lines."""
    p = os.path.join(RES, "tdg09", f"{gene}.tdg09.out")
    if not os.path.exists(p):
        return {}
    out, inblock = {}, False
    for line in open(p):
        if line.startswith("FullResults:"):
            inblock = True
            continue
        if not inblock:
            continue
        if line.startswith("- ["):
            parts = [x.strip() for x in line.strip()[3:-1].split(",")]
            if len(parts) < 2 or parts[-1] == "NA":
                continue
            try:
                out[int(parts[0])] = float(parts[-1])
            except ValueError:
                continue
        elif line.strip() and not line.startswith("#"):
            break
    return out


def load_contrastfel(gene):
    """{hyphy_site (1-based): q-value} for differential dN/dS."""
    p = os.path.join(RES, "selection", f"{gene}.contrastfel.json")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return {}
    d = json.load(open(p))
    mle = d.get("MLE", {})
    hdr = [h[0] for h in mle.get("headers", [])]
    if "Q-value (overall)" not in hdr:
        return {}
    qi = hdr.index("Q-value (overall)")
    return {i: r[qi] for i, r in enumerate(mle.get("content", {}).get("0", []), start=1)
            if qi < len(r)}


def trimmed_to_untrimmed(gene):
    """trimAl colnumbering: trimmed col (1-based) -> untrimmed col (1-based)."""
    p = os.path.join(RES, "trim", f"{gene}.colnumbering.txt")
    if not os.path.exists(p):
        return {}
    txt = open(p).read()
    if "#ColumnsMap" in txt:
        txt = txt.split("#ColumnsMap", 1)[1]
    cols = [int(x) for x in txt.replace("\n", " ").split(",") if x.strip().isdigit()]
    # trimAl writes 0-based original indices
    return {i: c + 1 for i, c in enumerate(cols, start=1)}


def codon_to_untrimmed(gene):
    """hyphy site -> untrimmed protein column (they are the same coordinate).

    The codon alignments were built on the UNTRIMMED protein alignment, so codon
    i corresponds to untrimmed protein column i; the map only accounts for the
    all-gap codon columns dropped before HyPhy could choke on them.
    """
    p = os.path.join(RES, "codon", f"{gene}.codonmap.csv")
    if not os.path.exists(p):
        return {}
    with open(p) as fh:
        return {int(r["hyphy_site"]): int(r["original_codon"]) for r in csv.DictReader(fh)}


def untrimmed_to_ref(gene):
    """untrimmed col -> REF_SPECIES residue number, by counting non-gap chars."""
    p = os.path.join(PROJ, "data", "alignments", f"{gene}_aligned.fa")
    if not os.path.exists(p):
        return {}
    seqs = read_fasta(p)
    row = seqs.get(REF_SPECIES)
    if row is None:
        return {}
    out, n = {}, 0
    for i, ch in enumerate(row, start=1):
        if ch not in "-.":
            n += 1
            out[i] = n
    return out          # gapped columns are absent: no residue in the reference


def main():
    os.makedirs(os.path.join(RES, "consensus"), exist_ok=True)
    rows = []
    tallies = []
    for g in GENES:
        pcoc = load_pcoc(g)
        tdg = load_tdg09(g)
        cfel = load_contrastfel(g)
        t2u = trimmed_to_untrimmed(g)
        c2u = codon_to_untrimmed(g)
        u2r = untrimmed_to_ref(g)

        # everything keyed on untrimmed column
        by_untrimmed = {}
        for col, post in pcoc.items():
            u = t2u.get(col)
            if u:
                by_untrimmed.setdefault(u, {})["pcoc_posterior"] = post
        for col, q in tdg.items():
            u = t2u.get(col)
            if u:
                by_untrimmed.setdefault(u, {})["tdg09_fdr"] = q
        for site, q in cfel.items():
            u = c2u.get(site)
            if u:
                by_untrimmed.setdefault(u, {})["contrastfel_q"] = q

        for u, d in sorted(by_untrimmed.items()):
            hits = {
                "pcoc": d.get("pcoc_posterior", 0.0) >= PCOC_THRESHOLD,
                "tdg09": d.get("tdg09_fdr", 1.0) <= FDR_THRESHOLD,
                "contrastfel": d.get("contrastfel_q", 1.0) <= FDR_THRESHOLD,
            }
            n = sum(hits.values())
            if n == 0:
                continue                      # only record sites some method flagged
            rows.append(dict(gene=g, untrimmed_col=u, ref_residue=u2r.get(u),
                             pcoc_posterior=round(d.get("pcoc_posterior", 0.0), 4),
                             tdg09_fdr=(round(d["tdg09_fdr"], 6) if "tdg09_fdr" in d else None),
                             contrastfel_q=(round(d["contrastfel_q"], 6) if "contrastfel_q" in d else None),
                             **hits, n_methods=n, high_confidence=n >= MIN_METHODS))
        tallies.append(dict(gene=g,
                            pcoc_sites=sum(1 for v in pcoc.values() if v >= PCOC_THRESHOLD),
                            tdg09_sites=sum(1 for v in tdg.values() if v <= FDR_THRESHOLD),
                            contrastfel_sites=sum(1 for v in cfel.values() if v <= FDR_THRESHOLD),
                            mapped_positions=len(by_untrimmed)))

    out_csv = os.path.join(RES, "consensus", "consensus_sites.csv")
    if rows:
        with open(out_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(sorted(rows, key=lambda r: (-r["n_methods"], r["gene"])))
    with open(os.path.join(RES, "consensus", "method_tallies.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tallies[0].keys()))
        w.writeheader()
        w.writerows(tallies)

    hc = [r for r in rows if r["high_confidence"]]
    print(f"{'gene':9s} {'PCOC':>6s} {'TDG09':>7s} {'C-FEL':>7s}")
    for t in tallies:
        print(f"{t['gene']:9s} {t['pcoc_sites']:6d} {t['tdg09_sites']:7d} {t['contrastfel_sites']:7d}")
    print(f"\nsites flagged by at least one method: {len(rows)}")
    print(f"HIGH CONFIDENCE (>= {MIN_METHODS} methods):  {len(hc)}")
    if hc:
        print(f"\n  {'gene':9s} {'ref_res':>8s} {'methods':>8s}")
        for r in hc[:20]:
            m = ",".join(k for k in ("pcoc", "tdg09", "contrastfel") if r[k])
            print(f"  {r['gene']:9s} {str(r['ref_residue']):>8s} {m:>8s}")
    else:
        print("\n  No site is supported by two or more methods. TDG09's hits stand alone,")
        print("  and with PCOC and Contrast-FEL both at zero none can be corroborated.")
    print(f"\nwrote {out_csv}")

    _figure(tallies, len(hc))


def _figure(tallies, n_hc):
    genes = [t["gene"] for t in tallies]
    fig = go.Figure()
    for key, name in (("pcoc_sites", "PCOC"), ("tdg09_sites", "TDG09"),
                      ("contrastfel_sites", "Contrast-FEL")):
        fig.add_trace(go.Bar(x=genes, y=[t[key] for t in tallies], name=name))
    fig.update_layout(
        template="simple_white", barmode="group", width=1150, height=520,
        title=(f"Sites flagged per method (high-confidence, {MIN_METHODS}+ methods: {n_hc})"),
        xaxis_title="gene", yaxis_title="sites flagged")
    os.makedirs(os.path.join(RES, "figures"), exist_ok=True)
    out = os.path.join(RES, "figures", "consensus_by_gene.pdf")
    fig.write_image(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
