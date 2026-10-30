# ionff: results in one page

Indranil Saha, 30 October 2026. Energies in kJ/mol. Numbers come from `results/fit_summary.json`, `results/cluster_summary.json`, `results/loso_neutralised.json` and `results/md_pooled.json`.

## Question

SAPT0 splits an interaction energy into electrostatics, exchange, induction and dispersion. Does fitting the same four analytic terms to each component, rather than to their sum, give a better Na+ force field?

## Setup

- Data: SAPT0/aug-cc-pVDZ, Na+ around six solvents (water, methanol, dimethyl ether, ethylene carbonate, dimethyl carbonate, acetone), 973 scan points used.
- Model: Slater Coulomb, Born-Mayer exchange, Tang-Toennies damped induction and dispersion. Fit by differential evolution plus bounded least squares, term-resolved or total-only.
- Tests: leave one solvent out, Na+(H2O)n for n = 1 to 6, and OpenMM MD of Na+ in SPC/E water (4 seeds).

## Results

| Na+(H2O)6, force field minus SAPT0 (SAPT0 total -403.2) | elst | exch | ind | disp | total |
|---|---|---|---|---|---|
| term-resolved | -119.7 | +34.0 | -16.6 | +0.7 | -101.7 |
| total-only | -44.9 | +36.6 | +57.3 | -144.7 | -95.8 |

1. Fitting by component exposes cancellation. The total-only fit reaches its total through induction (+57.3) and dispersion (-144.7) errors that cancel.
2. Correct components did not give a better total. Both fits miss the n = 6 total by about 100. The term-resolved error is mostly electrostatic and grows with n, not the induction non-additivity that was expected. The cause is untested.
3. Transfer fails as fitted. Held-out total RMSE is 92 to 374 against a reference spread of 28 to 52. After shifting type charges to neutrality it falls to 14 to 50, but the shift needs the new molecule's SAPT0 data, so it is not transfer.
4. MD in SPC/E water, term-resolved parameters: Na-O 0.239 nm (experiment 0.24), coordination 6.44 (experiment 4 to 8), D corrected 0.426e-9 m2/s against 1.33e-9. D is 3 to 4 times too low, and the finite-size correction is 67 to 87 percent of it.

## Limits

- Gas phase SAPT0, rigid monomers, pair-additive model.
- Cluster relaxation capped at 200 steps. n = 2, 5 and 6 are not all converged minima, so the n = 6 row is not at a minimum.
- 15 (term-resolved) and 17 (total-only) fitted parameters sit on a search bound.
- The MD box has net charge +4 e, with 4 Na+ and 4 seeds.
- The C++17 header in `cxx/` is tested against Python reference values but is not used by the fit or the MD.
- `results/REPORT.md` predates the neutralised transfer test. This page is current.

## Next

1. SAPT0 of Na+ against single waters from the cluster geometries, to separate non-additive electrostatics from a misfit of the single-water terms.
2. Assign charges of a held-out molecule without its SAPT0 data and repeat the transfer test.
3. Converge the clusters and rerun n = 2, 5 and 6.
4. Larger box and more Na+ before comparing D with experiment.
