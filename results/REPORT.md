# Results

## Fit to the SAPT0 scan (973 points)

RMSE in kJ/mol over all points. Term-resolved fits each SAPT0 component on its own; total-only fits the sum.

| fit | elst | exch | ind | disp | total |
|---|---|---|---|---|---|
| term_resolved | 25.28 | 2.90 | 1.80 | 0.08 | 25.13 |
| total_only | 41.05 | 39.96 | 46.78 | 26.94 | 13.62 |

## Leave one solvent out

Total-energy RMSE on the held-out solvent (kJ/mol).

| held out | term_resolved | total_only |
|---|---|---|
| acetone | 374.07 | 345.99 |
| dmc | 126.23 | 92.43 |
| dme | 210.62 | 195.99 |
| ec | 226.80 | 198.18 |
| methanol | skipped | skipped |
| water | skipped | skipped |

Skipped means the atom types of that solvent do not occur in the remaining training data.

## Na+(H2O)n many-body test

Mean of force field minus SAPT0 over relaxed and perturbed clusters (kJ/mol). The force field is summed pairwise over the waters.

| n | term_resolved total | term_resolved ind | total_only total | total_only ind |
|---|---|---|---|---|
| 1 | +0.8 | +2.1 | -2.7 | +21.1 |
| 2 | -7.4 | +2.4 | -13.2 | +39.5 |
| 3 | -19.4 | -0.7 | -26.4 | +53.9 |
| 4 | -41.2 | -1.7 | -46.3 | +68.3 |
| 5 | -61.1 | -7.6 | -63.6 | +76.2 |
| 6 | -101.7 | -16.6 | -95.8 | +57.3 |

## Na+ in SPC/E water, pooled over seeds

Mean over seeds with the standard error over seeds (each seed has 4 Na+). D is the MSD slope between 10 and 80 ps, without and with the Yeh-Hummer correction. Experiment: Na-O 0.24 nm, coordination 4 to 8, D = 1.33e-9 m2/s at 298 K.

| parameters | seeds | Na-O peak (nm) | first minimum (nm) | coordination | D raw (1e-9 m2/s) | D corrected (1e-9 m2/s) |
|---|---|---|---|---|---|---|
| term_resolved | 4 | 0.239 +/- 0.001 | 0.315 +/- 0.001 | 6.44 +/- 0.03 | 0.139 +/- 0.007 | 0.426 +/- 0.006 |
| total_only | 4 | 0.231 +/- 0.001 | 0.304 +/- 0.002 | 5.96 +/- 0.02 | 0.043 +/- 0.005 | 0.332 +/- 0.006 |

## Na+ in SPC/E water, per seed

Experiment: Na-O 0.24 nm, coordination 4 to 8, D = 1.33e-9 m2/s at 298 K. D is shown without and with the Yeh-Hummer finite-size correction.

| parameters | seed | Na-O peak (nm) | coordination | D raw (1e-9 m2/s) | D corrected (1e-9 m2/s) | block error (1e-9) |
|---|---|---|---|---|---|---|
| params_term_resolved.json | 1 | 0.242 | 6.5 | 0.12 | 0.41 | 0.01 |
| params_term_resolved.json | 2 | 0.237 | 6.5 | 0.15 | 0.43 | 0.01 |
| params_term_resolved.json | 3 | 0.237 | 6.4 | 0.14 | 0.43 | 0.03 |
| params_term_resolved.json | 4 | 0.237 | 6.4 | 0.15 | 0.44 | 0.02 |
| params_total_only.json | 1 | 0.233 | 6.0 | 0.05 | 0.34 | 0.00 |
| params_total_only.json | 2 | 0.228 | 6.0 | 0.04 | 0.33 | 0.00 |
| params_total_only.json | 3 | 0.233 | 5.9 | 0.03 | 0.32 | 0.01 |
| params_total_only.json | 4 | 0.233 | 6.0 | 0.06 | 0.34 | 0.01 |
