# Figure 4 source data

Updated on 19 October 2026 using the final mass-action parameters in SI S3.

- a-f: model curves and the experimental observations shown in the figure; measurements retain the original PFOA/PFOS source keys.
- g: all 1,158 calculated points satisfying both steric penalties <= 0.1. Every fifteenth eligible point in each series and the last point are displayed (86 symbols). The statistics use all eligible points and the fixed reference exp(-Delta Lambda), without fitting an exponential.
- h: 13 dilute-limit prediction ratios. An equal-weight linear regression of log10(ratio) against log10(I / 1 mM) gives intercept -1.79448185874, slope 0.655642357189, R-squared 0.941533601253, and RMSE in log10 coordinates 0.124806561002.
- `fig4_parameters.csv` records all 13 series. PFOS uses K1=1553 m3/mol, K2=3.8825 m6/mol2 and Gamma_m=3.0695e-6 mol/m2, consistent with the finalized SI. Langmuir parameters were retained; no adsorption-parameter fit was rerun.

The source-data rows are deterministic model calculations, not independent experimental replicates.

Reproduction: run `plot_selected_pfas_master_curve.py --output-dir <project figures directory>` from the analysis workspace.
