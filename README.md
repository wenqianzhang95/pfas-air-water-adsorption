# PFAS adsorption at the air-water interface: repository preparation

Authors: Wenqian Zhang and Bo Guo, University of Arizona.
Manuscript: Adsorption of ionic PFAS at the air-water interface at low concentrations.

**Status: private repository preparation, 19 October 2026. Hosted at https://github.com/wenqianzhang95/pfas-air-water-adsorption. No public release, DOI or redistribution licence has been assigned. This is a partial archive, not a complete submission archive.**

## Included material

- `data/fig4/`: processed model curves, measurements, fitted parameters, eligible/displayed points and regression results for Figure 4. Its README describes selection and statistics.
- `data/si/`: current comparison metrics, pointwise prediction ratios, dilute-limit results, low-concentration refit summary, current v5 transport and inverse results, and initial-profile-to-SSL ratios. The final manuscript mass-action parameter provenance table and its verification notes are included; obsolete parameter records and legacy error summaries are excluded.
- `data/transport_figure/`: profile and regime data associated with the transport figure. Cross-panel reconciliation remains to be completed.
- `data/transport_v5_checksums.json`: portable checksums for the three files covered by the existing v5 transport provenance record. All three were verified against that record during preparation. These checksums do not authenticate all analyses or establish scientific validity.
- `code/adsorption/`: Figure 4 plotting/model modules and their file-based PFOS inputs. The selected entry point is `plot_selected_pfas_master_curve.py`; the base module contains other workflows and must not be treated as the final standalone figure entry point.
- `code/transport/`: current nonlinear forward/inverse implementation, inputs and selected launch helpers. The original `README_SI_FINAL.md` describes the source package and includes commands requiring additional HPC packaging; it is retained as provenance, not as a tested installation guide for this subset.
- `metadata/manuscript_references.bib`: original literature references. Experimental observations reused from publications retain their original citations; this archive does not grant rights over third-party material.
- `metadata/python_preparation_environment.json`: installed Python library versions on the preparation machine; these are not minimum supported versions or a validated environment specification.
- `file_manifest.csv`: SHA-256, size, archive-relative path and category for each copied file, plus preparation-generated metadata.

## Penalty-tolerance sensitivity

SI Section S9 is supported by `data/si/regime_threshold_sensitivity/`, the matching figure in `figures/supplementary/`, and `code/sensitivity/analyze_regime_threshold_sensitivity.py`. The source-data README documents units, numerical domains, exact boundary definitions, statistical formulas and reproduction commands. This analysis tests tolerances 0.05/0.10/0.20 without refitting. Its complete inclusion does not resolve the other pending MD/calibration archive gaps listed below.

## Figure 4 entry point

With Python, NumPy, pandas, SciPy and Matplotlib installed, run from `code/adsorption`:

```text
python plot_selected_pfas_master_curve.py --output-dir ../../results/fig4
```

The command and its relative input paths have been inspected. A clean-environment regeneration of this staged subset remains pending.

## Transport dependency

The MATLAB implementation uses MRST 2022a. MRST is not included in this subset. Its acquisition route, licence, modules and MATLAB version must be documented and a clean installation tested before release. Keep the supplied XLSX files with the script. Existing source-package HPC launch commands may require files outside this subset.

## Required before release

1. Reconcile every main and supplementary panel with the final calculations, inputs and scripts. MD input files, final trajectories or shareable processed profiles, and their analysis code have not yet been assembled here. Other calibration/regime/potential analyses are not completely covered.
2. Complete portable reproduction instructions and clean-environment checks. Source-data inclusion is not proof that all results can currently be regenerated.
3. Resolve third-party data and bundled-software redistribution terms; select licences for author-owned data and code. No licence is inferred here.
4. Complete and check this private repository, then publish and archive a versioned release; insert the actual version and permanent identifier into manuscript availability statements. A Zenodo archive is planned but has not been created.
5. Supply metadata and formal dataset citation after the repository assigns its identifier. Do not replace pending fields with invented links.
