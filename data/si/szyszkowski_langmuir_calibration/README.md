Figure S3 calibration source data

parameters.csv contains all 13 recorded calibrations (PFOA: 10/50/100/200/300/500 mM; PFHxA: 100/200/300/500 mM; PFOS: 1/10/50 mM). all_series_parameters.csv supplies the same parameter values for the SI table. No parameters were refitted.

observations.csv contains the 98 retained surface-tension observations; curves.csv contains the plotted predictions. The original 70 observations are unchanged; 28 observations were added for the four previously omitted series. PFOA/PFHxA fits and points come from the archived all-ionic-strength calibration files. PFOS fits match the recorded parameter sets used for Figure 4; RMSE values were reproduced from the source observations.

PFOS selection: SCH2019_NaCl1mM, 5 points (highest-concentration sixth point excluded in the recorded fit); HUA2022_NaCl10mM_Set1, 10 points; STE2022_NaCl50mM, 5 points, KPFOS potassium-salt series. Surface-tension calibration observations are distinct from the direct surface-excess data reserved for evaluation.

Axis conversions retain the source conventions: PFOA 414, PFHxA 314, PFOS 500.075 g/mol. PFOS curves span the full displayed axis; PFOA/PFHxA curves span each retained calibration range.

Colours use the aurora-8 palette in palette.json, with an ordered, shared seven-concentration mapping across all panels. Marker shapes supply a redundant cue.

Reproduce: python SI/code/plot_szyszkowski_langmuir_calibrations.py
