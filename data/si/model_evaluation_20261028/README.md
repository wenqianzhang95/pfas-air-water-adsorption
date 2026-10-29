# Evaluation with the currently declared manuscript parameter sets

The script `SI/code/evaluate_model_predictions.py` evaluates the current SI
parameter sets without refitting, updating original observations or redrawing
figures. Run it from the Nature_Water paper workspace.

- `parameter_sets_used.csv`: the declared parameters copied from the provenance table.
- `surface_tension_pointwise.csv`: all 114 archived Fig. 2 observations and
  direct predictions at their concentrations; residuals are prediction minus observation.
- `surface_tension_summary.csv`: per-series and selected pooled unweighted RMSDs.
- `binding_pointwise.csv`: 36 distinct numerical MD records and predictions.
- `binding_summary.csv`: per-compound and pooled unweighted RMSEs in percentage points.
- `exported_surface_curve_diagnostics.csv`: errors obtained by interpolating
  historical plotted polylines; diagnostics only, not the reported current-parameter errors.
- `current_vs_exported_surface_curves.csv`: differences between the current
  declared model and the historical plotted surface-tension curves.
- `current_vs_exported_binding_curves.csv`: historical plotted binding constants
  and their errors alongside the current declared constants.
- `verification.json`: source hashes, definitions and numerical solver checks.

Calibration/prediction designations and literature attribution follow the
current manuscript. The PFOS 50 mM input is specifically marked as a comparison
read from an excluded-data source by the original exporter; no pooled PFOS
calibration statistic is reported.

The plotted historical curves are not all generated with the currently
declared parameters. This discrepancy is preserved and documented, not hidden
by overwriting their source data. A separate independent Gouy--Chapman test
was not performed; numerical solver verification is not physical validation.
