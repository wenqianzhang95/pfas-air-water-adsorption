# Final manuscript mass-action parameter provenance

Updated 19 October 2026 against the current `supporting_information.tex`.

The CSV records the final SI Section S3 parameter sets and their verified use in Figure 4 source data and current v5 forward/inverse transport audits. Repeated compound rows describe different checked uses of the same parameter set, not separate fits.

- K1 = K_ad, in m3/mol.
- K2 = K_ad * K_bind, in m6/mol2.
- K_bind = K2/K1, in m3/mol; this column is calculated, not independently fitted.
- Gamma_m is in mol/m2; multiply by 1e6 to obtain the micromol/m2 values displayed in SI S3.
- Concentrations for these concentration-based constants are in mol/m3. These are mass-action constants, not the salt-specific Langmuir parameters.
- PFHxA has no transport row because the reported forward/inverse comparison is for PFOA and PFOS.

Historical fitting and MD master-curve records were removed from this final-parameter table. The former table included obsolete values (including PFHxA K1=10.4578 and PFOA K2=0.13284) and stale status labels. Its exact prior version is retained in the local editorial backup, not presented as current source data.

Verification covers numerical agreement with SI S3, all 13 Figure 4 parameter-series rows, and every row of both v5 transport audits. It does not establish that historical calibration/MD plotting scripts or previously exported figures have been regenerated with these final constants. MD master-curve provenance and other figure-specific workflows require their own checks; changing this CSV does not change underlying calculations.

Paths in the CSV are relative to the manuscript project. In the repository preparation layout, Figure 4 data are under `data/fig4/`, and the SI datasets are under `data/si/`. The parameter declaration is traceable to the SI Section S3 table rather than to a refit performed during this update.
