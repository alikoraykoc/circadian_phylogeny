# Gene-level rate association (RERconverge)

- `rer_gene_association.csv` step 11 output, parametric p-values from
  `correlateWithBinaryPhenotype`.
- `rer_permulation_pvalues_ssm.csv` **the file to quote.** Same association
  re-tested against 1000 permulations in species-subset-match mode, which
  rebuilds the null against each gene's own taxon set. Gene occupancy ranges from
  40 to 60 taxa here, so this mode is the appropriate one.
- `rer_permulation_pvalues_cc.csv` complete-case permulations, run as a
  robustness check. Agrees with ssm throughout.

Parametric p-values are not calibrated under phylogenetic non-independence. On
this tree they were anticonservative in 10 of 18 genes. BHLHE40, the only gene
whose parametric p-value approached significance at 0.051, has a permulation
p-value of 0.114. No gene reaches 0.05 under either permulation mode; the
smallest adjusted values are 0.728 (ssm) and 0.764 (cc).

Method: Saputra, Kowalczyk et al., Molecular Biology and Evolution 2021.
Script: `scripts/11b_rerconverge_permulations.R`.
