# Adsorption-penalty tolerance sensitivity

Deterministic sensitivity analysis for SI Section S9, with final SI S3 parameters and no refitting. Tolerances: 0.05, 0.10 and 0.20. Background NaCl is fixed for each concentration--salinity map calculation, following the manuscript regime-map construction. The archived Figure 4 curves are reused for the separate low-coverage point-selection comparison.

## Files and units

- `threshold_regime_boundaries.csv`: 5,400 rows (3 compounds x 3 tolerances x 600 NaCl values). `Surface_activity_exit_mg_L` is the exact root of max(Lambda_el,Lambda_st)=tolerance; `Limiting_penalty` identifies the term reaching the tolerance first. `Equal_penalty_concentration_mg_L` and `Penalty_at_equality` locate the nonzero electrostatic--steric equality locus; these columns are identical across tolerances. Blank equality coordinates mean that no positive electrostatic-dominance interval exists at that salinity, not that data were lost. Coordinates may lie outside the plotted concentration interval.
- `threshold_selected_concentrations.csv`: exact surface-activity exit coordinates at 1, 10, 100, 200, 300, 500 and 1,000 mM NaCl.
- `threshold_regime_summary.csv`: class fractions and reclassification percentages on the same 700 concentration x 600 salt log-uniform grid (10^-4 to 10^4 mg/L; 1 to 1,000 mM). These are domain-specific numerical summaries, not environmental occurrence probabilities. `Terminal_salt_for_distinct_electrostatic_regime_mM` is a continuous-root estimate of the salinity at which the equality point reaches the tolerance, so the labelled intermediate band disappears; it is not an independently measured critical concentration.
- `threshold_transfer_pointwise.csv`: all 2,240 archived Figure 4 curve points, parameters, derived penalties, fixed theoretical references and eligibility at each cutoff. No points are redigitized or new curve fits performed.
- `threshold_transfer_summary.csv`: eligible-point count, R-squared and RMSE against the fixed reference exp(-Delta Lambda), and maximum absolute and relative deviations. R-squared = 1 - sum((ratio-reference)^2)/sum((ratio-mean(ratio))^2); RMSE = sqrt(mean((ratio-reference)^2)). Neither uses logarithmic coordinates. All eligible points contribute; displayed-point subsampling is not used for these statistics.
- `verification.json`: fixed parameters, domains, grid dimensions, input-file SHA-256 checksums, onset-equation residual and successful reproduction of the original 0.10 Figure 4g count and statistics.

K1 is in m3/mol; K2 in m6/mol2; Gamma_m in mol/m2. Concentrations used by the equations are mol/m3; multiply by molecular weight (g/mol) to obtain mg/L. Map molecular weights follow the existing plotting workflow (PFOA 414, PFHxA 313, PFOS 499 g/mol). The Figure 4 comparison preserves its archived curves without converting their concentrations again.

## Reproduction

Python with NumPy, pandas, SciPy and Matplotlib is required. From a manuscript project, run:

```text
python SI/code/analyze_regime_threshold_sensitivity.py --project . --make-figure
```

From the repository preparation layout, run:

```text
python code/sensitivity/analyze_regime_threshold_sensitivity.py --project results/threshold_sensitivity --fig4-source data/fig4 --make-figure
```

The latter command writes its output under `results/threshold_sensitivity/SI/`. Figure exports are PDF, SVG and PNG. The main-text maps and the production transport calculations are not overwritten. This tolerance controls labels and low-coverage point selection, not the adsorption or transport equations. No uncertainty intervals or experimental replicates are implied by the numerical grid.
