**Table 1. Gene set, taxon occupancy and alignment dimensions.**

Trimmed columns are those retained by trimAl `-automated1`. A column is variable when it holds more than one amino acid, and parsimony-informative when at least two amino acids each occur in at least two sequences. The parsimony-informative columns are exactly the sites TDG09 could test; it returns NA at every other column. Codons analysed are those remaining after removal of all-gap columns.

| Gene | Module | Taxa | Untrimmed columns | Trimmed columns | Variable columns | Parsimony-informative sites | Codons analysed |
|---|---|---|---|---|---|---|---|
| CLOCK | Positive arm | 59 | 991 | 841 | 308 | 177 | 991 |
| NPAS2 | Positive arm | 58 | 1145 | 814 | 642 | 351 | 1120 |
| ARNTL | Positive arm | 60 | 778 | 620 | 136 | 50 | 750 |
| PER1 | Negative arm | 60 | 1643 | 1233 | 581 | 409 | 1584 |
| PER2 | Negative arm | 60 | 1726 | 1208 | 971 | 801 | 1704 |
| PER3 | Negative arm | 56 | 2122 | 897 | 829 | 637 | 1640 |
| CRY1 | Negative arm | 59 | 670 | 582 | 126 | 80 | 666 |
| CRY2 | Negative arm | 60 | 665 | 592 | 251 | 211 | 661 |
| NR1D1 | Auxiliary loop | 60 | 760 | 610 | 321 | 144 | 731 |
| NR1D2 | Auxiliary loop | 60 | 660 | 562 | 281 | 156 | 660 |
| RORA | Auxiliary loop | 58 | 602 | 466 | 69 | 25 | 602 |
| RORB | Auxiliary loop | 60 | 665 | 459 | 117 | 34 | 609 |
| RORC | Auxiliary loop | 60 | 790 | 518 | 266 | 210 | 753 |
| CSNK1D | Post-translational | 40 | 792 | 401 | 63 | 12 | 553 |
| CSNK1E | Post-translational | 52 | 564 | 509 | 197 | 174 | 564 |
| FBXL3 | Post-translational | 60 | 442 | 426 | 56 | 34 | 442 |
| BHLHE40 | Output repressors | 60 | 666 | 407 | 200 | 132 | 599 |
| BHLHE41 | Output repressors | 58 | 768 | 582 | 343 | 183 | 720 |

**Table 2. Convergent scenarios, calibrated detection power and outcome.**

Branch and leaf counts are on the 60-taxon species tree (118 branches). Power and false positive rate are the worst value across the 18 genes at the calibrated posterior threshold of 0.99. All sites above threshold were subsequently rejected as alignment artefacts.

| Reconstruction | Direction | Events | Convergent branches | Percent of branches | Convergent leaves | Minimum gene power | Worst gene FPR | Sites above threshold |
|---|---|---|---|---|---|---|---|---|
| ER | Gains of diurnality | 10 | 55 | 46.6 | 41 | 1.000 | 0.0000 | 0 |
| ER | Reversals to nocturnality | 5 | 19 | 16.1 | 13 | 0.996 | 0.0000 | 1 |
| ARD | Gains of diurnality | 6 | 64 | 54.2 | 55 | 0.990 | 0.0000 | 2 |
| ARD | Reversals to nocturnality | 10 | 44 | 37.3 | 29 |  |  | not run |

**Table 4. Gene and site concordance at every internal transition branch.**

The remaining 10 of the 15 transitions are on tip branches, where a concordance factor is undefined because a tip appears in every gene tree. A branch is flagged as a hemiplasy risk only when gCF and sCFL are both below 50 percent, since either alone is usually gene-tree estimation error rather than genuine conflict. No transition is flagged.

| Event | Direction | Descendant taxa | gCF | sCFL | Decisive gene trees | Weak axis | Flagged |
|---|---|---|---|---|---|---|---|
| 1 | gain | 24 | 38.9 | 75.7 | 18 | gCF | False |
| 3 | reversal | 8 | 50.0 | 42.0 | 18 | sCFL | False |
| 8 | gain | 8 | 55.6 | 60.8 | 18 | none | False |
| 2 | reversal | 2 | 56.2 | 48.1 | 16 | sCFL | False |
| 11 | gain | 4 | 58.8 | 55.4 | 17 | none | False |

**Table 5. Convergent sites as a function of the phenotype-specificity threshold.**

Sites counted are those that are both statistically unusual (q <= 0.05) and phenotype-specific at the given threshold. The 0.25 cutoff was fixed by judgment before the positive control existed and never adjusted; the row between 0.25 and 0.20 is the largest gap observed among the statistically unusual sites, so every threshold above it returns zero.

| phenotype-specificity threshold | sites passing both criteria |
|---|---|
| 0.3 | 0 |
| 0.25 (used here) | 0 |
| 0.232 | 1 |
| 0.2 | 2 |
| 0.15 | 2 |
| 0.1 | 4 |
| 0 | 9 |
