# Adsorption slopes and dilute limits

The six plotted conditions match main-text Figure 4a-f. No parameters are refitted.

- `adsorption_slope_curves.csv`: concentrations, MA surface excess, and analytical local derivatives dGamma/dc for both models; derivatives use c in mol/m3 and have units m.
- `adsorption_slope_dilute_limits.csv`: analytical dilute-limit coefficients for all 13 archived PFAS/NaCl parameter series, in m. Their SL/MA ratios agree with the existing empirical dilute-limit source records.
- `adsorption_slopes_verification.json`: numerical verification against central finite differences, analytical limits, and archived Figure 4 surface excesses.
- `adsorption_slope_calibration_ranges.csv`: minimum and maximum retained surface-tension concentrations, separate shading limits starting at the second-lowest distinct concentration, counts and source references for each of the six SL fits.
- `adsorption_slope_calibration_observations.csv`: the individual surface-tension observations and model-calculated surface excesses and local slopes at those concentrations. The calculated slopes are not experimental measurements.
- `../code/plot_adsorption_slopes.py`: reproduction script. Run with Python from any working directory.

Background NaCl is fixed. MA bulk counterion concentration is NaCl + dissolved PFAS concentration, following Figure 4; the accompanying counterions vanish from this sum in the dilute limit. The figure uses the same anion molecular weights for the mass-concentration axes as Figure 4 (413, 313, 499 g/mol). The slope in m is independent of this mass-axis conversion.

Finite-concentration derivatives are local slopes, not Gamma/c. Both converge to k_aw as c tends to zero. Calculations are deterministic model outputs, not new independent measurements. The figure is exported as PDF, SVG and PNG. The three rows use paired equal y ranges, stable model colors, line-style redundancy and direct plateau annotations.

Grey bands run from the second-lowest to the highest distinct retained surface-tension concentration for the corresponding Szyszkowski--Langmuir fits. This display convention assigns the slope between adjacent observations to the upper concentration of the interval. It is not a claim that a finite-difference slope is the exact local derivative at that endpoint. Top-edge ticks retain all observation concentrations, including the first point. All original observations remain in the fit; no fitting or model calculations change. These ranges refer to SL calibration, not necessarily the joint mass-action calibration. For example, PFOA at 10 mM is an independent surface-tension evaluation of the mass-action model. Molar observation concentrations are converted with the existing figure-axis molecular weights above, avoiding the different acid/salt mass conventions in the calibration figure. Axes extend to show the complete retained ranges: 10,000 mg/L for PFOA, 100,000 mg/L for PFHxA and 1,000 mg/L for PFOS. The shaded ranges indicate the adjacent-slope display convention and do not imply equal statistical influence of all observations or agreement of the two adsorption predictions throughout the bands.
