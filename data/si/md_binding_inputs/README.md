# MD counterion binding source data

The snapshots in this directory preserve the input arrays used for the MD
points in main-text Figure 2d-f (the figure file is `fig3_surface_tension.pdf`).
The canonical manuscript figure and the working figure from which the point
CSV was exported have the same SHA-256:
`5702a76e4f8d1a527d7abd46ca317dc4d9c779ec981f0eb383bdb4a41e4312a2`.

Original inputs, relative to the research workspace:

- `PFAS_plotResults_package/data/PFOA/plotResults_panels_data.json`
- `PFAS_plotResults_package/data/PFHxA/plotResults_panels_data.json`
- `PFAS_plotResults_package/data/PFOS/plotResults_panels_data.json`
- `NatureWater_Overleaf/figures/fig3_surface_tension_source_data.csv`

`main_binding_plot_points.csv` contains only the `Binding MD` rows of the
last file. All 37 plotted records are checked against the source JSON arrays,
including NaCl concentration, local counterion concentration, binding ratio,
and both error magnitudes. Surface excesses come from the corresponding JSON
entries and are converted from mol/m2 to micromol/m2. At the author's request,
the SI table displays the original 499 and 997 mM labels as 500 and 1000 mM,
respectively. The CSV retains the original concentrations in `NaCl_mM` and
provides the displayed labels in `NaCl_display_mM`. The archived inputs and
all binding-ratio and local counterion-concentration values are unchanged.

The PFHxA input contains an exact duplicate at source indices 7 and 11
(1-based). It is plotted twice at the same coordinates and listed once in the
36-row numerical-source portion of the SI table. Both indices are retained in
the full-precision CSV.
The PFOS input contains a point excluded by the existing figure-generation
script. It is retained in `md_counterion_binding_source_records.csv` with its
exclusion status, but is not added to the displayed-data table. No exclusion,
fit, adsorption analysis, or simulation was performed for this table.

Error magnitudes are reproduced as supplied in `c2sSTD` and `STD`. The author
specifies that binding ratios are ensemble averages and that error bars denote
standard deviations. The table caption follows this supplied definition; no
error magnitudes are recalculated or converted.

The original numerical source data for main-text Figure 1g,h have not been
located on this computer. The author indicated that they may be on another
computer and subsequently requested provisional extraction from the original
figure. `md_binding_figure1_digitized.csv` therefore contains 15 provisional
binding ratios and error-bar magnitudes recovered directly from PDF vector
geometry by `SI/code/digitize_md_binding_figure1.py`. Bar colors identify the
PFAS, axis ranges are 0-100% in panel g and 0-50% in panel h, and bar/error-cap
coordinates determine the values. Axis limits, color identities, and extracted
values were checked against a rendered crop of the original panels. The source
PDF hash and extraction coordinates are retained in the extraction audit under
`output/editorial_review_20261019/md_binding_table_revision_20261021/`.

These values describe the plotted figure, not the original simulation output.
Figure-derived binding ratios and errors are displayed to one decimal place
in the SI; cross-figure duplicates are represented by their numerical sources.
All 43 entries are combined in one table, ordered by PFAS, displayed NaCl
concentration and unrounded surface excess. Figure-derived entries retain a
dash for missing local counterion concentrations; they have no separate group
heading or dagger in the displayed table. Their provenance remains recorded
in the source CSVs and reconciliation file.
Original numerical values remain blank in `md_binding_figure1_pending.csv`,
which now also includes the provisional extracted estimates. The 15 condition
labels are transcribed from the figure and caption only, not independently
verified simulation loadings. Local counterion concentrations are not shown in
Figure 1 and are left blank.
No missing values were substituted from other conditions or from the separate
GROMACS analysis project. Figure-derived error magnitudes remain as drawn;
their original numerical sources are still pending.

At the author's request, repeated means across the two figures are listed only
once. Matching requires the same PFAS, the same NaCl condition and a unique
binding-ratio match within 0.05 percentage points (the rounding precision of the
one-decimal Figure 1 means). Eight Figure 1 entries match available numerical
records: three PFOA, three PFHxA and two PFOS entries, all at 100 mM. The table
retains numerical-source Gamma, binding ratio and errors for these records.
The Figure 1 error bars differ from the numerical-source errors; their values
remain archived, and no SD/SEM conversion is performed.

Five matches have approximately corresponding Gamma labels. The three PFHxA
matches have substantial Gamma discrepancies that cannot be explained by
rounding: Figure 1 labels 0.04, 0.26 and 0.47 micromol/m2 correspond to numerical
source Gamma values 0.0673, 0.202 and 0.359 micromol/m2, respectively. Their means
match within 0.0001 percentage points. These repeated means are consolidated as
requested, but the Figure 1 Gamma labels require checking against the original
simulation records. This table revision does not change the manuscript figure.
All 15 Figure 1 readings are retained in the digitized and pending CSVs.
`md_binding_cross_figure_reconciliation.csv` records every match and label
conflict. The table now contains 43 rows: 36 numerical-source records and seven
additional figure-derived records (six at 20 mM and one PFOS record at 100 mM).

To reproduce the figure extraction, run `SI/code/digitize_md_binding_figure1.py`.
Then run `SI/code/assemble_md_binding_table.py` from any directory to rebuild the
longtable and CSVs. All surface excesses in the table use two decimal places,
and NaCl concentrations use a thousands separator (1,000 mM).
Other numerical-source values use three significant figures, and provisional
digitized binding ratios and error magnitudes use one decimal place. Original numerical CSV values
are unrounded; geometry-derived estimates retain their computed coordinates
without implying equivalent accuracy in the underlying simulation results.
