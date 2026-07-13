#!/usr/bin/env bash
# 06_concordance.sh
# Annotate every branch of the Upham tree with gene concordance factor (gCF) and
# likelihood-based site concordance factor (sCFL). Low-gCF branches are where
# gene-tree/species-tree discordance concentrates, i.e. the branches most exposed
# to hemiplasy-driven false convergence. These values gate the scenario in step 08.
#
# Two separate IQ-TREE runs:
#   1) gCF from the 18 gene trees (--gcf)
#   2) sCFL from the partitioned locus alignments (--scfl)
# Results are merged by branch ID into a single table.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

# 1) Gene concordance factors
echo "== gCF from gene trees =="
iqtree -t "$SPTREE" \
        --gcf "$RES/genetrees_qc/all_gene_trees.treefile" \
        --prefix "$RES/concordance/gcf" -T 4 -redo

# 2) Likelihood-based site concordance factors
echo "== sCFL from locus alignments =="
iqtree -t "$SPTREE" \
        -p "$RES/supermatrix/parts" --scfl 100 \
        --prefix "$RES/concordance/scfl" -T 4 -redo

# 3) Merge gCF and sCFL stat tables by branch ID
echo "== Merging gCF + sCFL by branch ID =="
python3 - "$RES/concordance/gcf.cf.stat" "$RES/concordance/scfl.cf.stat" \
          "$RES/concordance/concord.cf.stat" << 'PYEOF'
import sys, csv

gcf_file, scfl_file, out_file = sys.argv[1], sys.argv[2], sys.argv[3]

def parse_cf_stat(path):
    rows = {}
    with open(path) as f:
        header = None
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if header is None:
                header = parts
                continue
            row = dict(zip(header, parts))
            try:
                bid = int(row["ID"])
                rows[bid] = row
            except (ValueError, KeyError):
                continue
    return rows

gcf_data = parse_cf_stat(gcf_file)
scfl_data = parse_cf_stat(scfl_file)

all_ids = sorted(set(gcf_data.keys()) | set(scfl_data.keys()))

with open(out_file, "w") as f:
    f.write("ID\tgCF\tgCF_N\tgDF1\tgDF2\tgDFP\tgN\tsCFL\tsCFL_N\tsCF_DF1\tsCF_DF2\tsCF_N\n")
    for bid in all_ids:
        g = gcf_data.get(bid, {})
        s = scfl_data.get(bid, {})
        gcf = g.get("gCF", "NA")
        gcf_n = g.get("gCF_N", g.get("gN", "NA"))
        gdf1 = g.get("gDF1", "NA")
        gdf2 = g.get("gDF2", "NA")
        gdfp = g.get("gDFP", "NA")
        gn = g.get("gN", "NA")
        # sCFL columns vary by IQ-TREE version; try common names
        scfl = s.get("sCFL", s.get("sCF", "NA"))
        scfl_n = s.get("sCFL_N", s.get("sCF_N", "NA"))
        sdf1 = s.get("sDF1", s.get("sCF_DF1", "NA"))
        sdf2 = s.get("sDF2", s.get("sCF_DF2", "NA"))
        sn = s.get("sN", s.get("sCF_N", "NA"))
        f.write(f"{bid}\t{gcf}\t{gcf_n}\t{gdf1}\t{gdf2}\t{gdfp}\t{gn}\t{scfl}\t{scfl_n}\t{sdf1}\t{sdf2}\t{sn}\n")

print(f"Merged {len(all_ids)} branches -> {out_file}")
PYEOF

# Expose two trees for step 08:
#   concord.cf.tree    internal labels = support/gCF (human-readable)
#   concord.cf.branch  internal labels = integer branch IDs matching concord.cf.stat
# Step 08 matches transition branches by bipartition, so it needs the branch-ID
# tree; the .stat IDs come from the gCF run, so copy that run's .cf.branch.
cp "$RES/concordance/gcf.cf.tree"   "$RES/concordance/concord.cf.tree"
cp "$RES/concordance/gcf.cf.branch" "$RES/concordance/concord.cf.branch"

echo
echo "Outputs:"
echo "  concord.cf.tree      Upham tree with gCF annotations (from gcf run)"
echo "  concord.cf.stat      merged gCF + sCFL per-branch table (feeds step 08)"
echo "  gcf.cf.stat          raw gCF output"
echo "  scfl.cf.stat         raw sCFL output"
