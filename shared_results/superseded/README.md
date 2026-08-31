# Superseded results

Everything in this directory is kept for provenance only.
**Do not quote any number from these files.** Each entry below names the file
that replaced it.

The project was re-run twice after two defects were fixed, so several quantities
have three historical values. All three generations support the same conclusion,
but only the current generation should be cited.

| generation | what was wrong | fixed in |
|---|---|---|
| pre-rooting-fix | per-gene trees were rooted inconsistently with the species tree, so transition branches sat at the wrong nodes | commit `72b42b2` |
| v1, post-rooting-fix | correct rooting, but the parsimony test took one arbitrary resolution per site and used an unmatched null | test hardening |
| current | in use | |

## parsimony_runs/

Superseded runs of `scripts/observed_convergent_substitutions.py`.
Current file: `../pcoc_sim_calibration/observed_convergent_substitutions.csv`.

| file | sites changed in 2+ diurnal lineages | same residue | observed rate | null rate |
|---|---|---|---|---|
| `observed_convergent_substitutions.prereroot.csv` | 639 | 206 | 0.322 | 0.336 |
| `observed_convergent_substitutions.v1.csv` | 645 | **213** | 0.330 | 0.364 |
| current, not here | 648 | **201** | 0.310 | 0.337 |

If you have met **213** or **206** anywhere, this is where it came from. The
current value is **201**.

The 213 to 201 drop is the test hardening, not a change in the data. Resolving
parsimony ties 10 ways per site, instead of accepting one arbitrary resolution,
drops sites whose apparent convergence depended on which equally parsimonious
history happened to be chosen. The site count rose slightly, 645 to 648, because
the matched null changed at the same time.

`observed_convergent_sites.v1.csv` is the matching per-site table: 213 rows and
5 columns. The current one has 201 rows and adds `diurnal_gap`,
`phenotype_specific`, `p_site` and `q_site`, the columns that carry the
phenotype-specificity filter.

`transition_opportunity.prereroot.csv` is the pre-rooting-fix opportunity
estimate. Current file: `../pcoc_sim_calibration/transition_opportunity.csv`.

**Telling them apart from the inside:** superseded substitution tables carry a
`p_count_confounded` column. The current table carries `n_perm`,
`n_resolutions` and `null_draws_needing_relaxed_bl` instead.

## partial_convergence_prereroot/

The partial convergence k-curve, computed 2026-08-13 before the rooting fix.
Current directory: `../partial_convergence/`.

This one was actively misleading and is the reason this archive exists. Until
now the superseded curve held the plain directory name `partial_convergence/`
while the corrected re-run sat under `partial_convergence_rerun/`, which is the
opposite of what a reader would assume. The names have been swapped.

**Telling them apart from the inside:** the corrected curve has 55 convergent
branches at k=10, matching the 55 in Table 2 for ER gains of diurnality. The
superseded curve has 52.

| converging lineages of 10 | power at 0.9, current | power at 0.9, superseded |
|---|---|---|
| 4 | 0.000 | 0.007 |
| 5 | 0.012 | 0.020 |
| 6 | 0.062 | 0.088 |
| 7 | 0.258 | 0.250 |
| 8 | 0.577 | 0.548 |

The re-run confirmed the published curve rather than overturning it; the shape
and the conclusion are unchanged. That was not knowable in advance, which is why
the re-run was done.
