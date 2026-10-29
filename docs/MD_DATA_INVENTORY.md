# MD data inventory for the archive (status 2026-10-28)

Purpose: record where the molecular dynamics material behind main-text Figures 1 and 2
and Supplementary Section S1 currently lives on the authors' machine, so the authors can
decide what to deposit. Nothing here has been copied into the repository yet; paths are
local to the preparation machine and are not part of the public archive.

## What the manuscript uses

- Engine and force field (SI S1): LAMMPS (27 Oct 2021 build per available logs), SPC/E water,
  OPLS-AA PFAS, Dang-Smith Na+/Cl-, NVT 300 K, PPPM, SHAKE on water, 5 ns + 5 ns.
- 36 distinct MD records (37 plotted points, one PFHxA duplicate) for PFOA, PFHxA, PFOS at
  20-1,000 mM NaCl. The numerical source of the Figure 2d-f points is
  `data/si/md_binding_inputs/*_plotResults_panels_data.json` (exported arrays, not raw output).
- Figure 1g,h values: the original numerical files have not been located; seven table rows
  were digitized from the figure PDF (`data/si/md_binding_figure1_digitized.csv`,
  `data/si/md_binding_figure1_pending.csv`).

## What exists locally

| Location (preparation machine) | Content | Size | Candidate for archive |
|---|---|---|---|
| `E:/PFAS/MD/binding model/low ionic strength (I = 20)/` | LAMMPS input (`test.lmp`), data file, `log.lammps`, `*.profile` (Na, Cl, O, water, PFAS), `rdf_*.rdf`, `out.charge_density`, `out.surfTension.txt`, `trajectory.data` (~337 MB) | 426 MB | inputs, logs, profiles, RDFs, charge density: yes (few MB). Trajectory: Zenodo only or on request |
| `E:/PFAS/MD/binding model/I = 100 midle surface coverage/` | same layout | 349 MB | same |
| `E:/PFAS/MD/binding model/I = 500 middle surface coverage/` | same layout | 430 MB | same |
| `E:/PFAS/MD/binding model/high coverage/` | same layout (trajectory ~61 MB) | 117 MB | same |
| `E:/PFAS/MD/binding model/36+100/` | LAMMPS case | 18 MB | yes |
| `E:/PFAS/MD/binding model/*.mlx`, `MCmethod.mlx`, `bindingRatioPlotting*.mlx`, `RawData.xlsx`, `Data set.xlsb.xlsx` | MATLAB analysis of profiles -> binding ratio, RDF; raw data tables | small | analysis code and tables: yes, after converting `.mlx` to `.m` for readability |
| `E:/PFAS/MD/air-water interface/` | PFOS air-water interface LAMMPS setup and NPT/NVT stages | 320 MB | inputs only |
| `E:/PFAS/MD/PFOA/`, `E:/PFAS/MD/PFOS/` | force-field files (`PFOA.lj`, `PFOA.lmp`), validation | 1.6 GB / 0.2 MB | force-field definitions: yes |
| `E:/PFAS/MD/PBE/` (`PFOA/PFOA16+100`, `PFHxA/PFHxA16+100_20mM`, `PFOS/...`, plus many other PFAS) | **GROMACS** campaign used by the companion EST manuscript (`binding_ratio_workflow_package/`); `density*.xvg`, `ff/`, `input/`, `run/` per case | large | **not** the Nature Water source unless a provenance link is verified (see `MD_S1_methods_audit_20261019.md` in the editorial-review folder). Do not mix the two campaigns. |

Only five LAMMPS case directories were found locally, whereas the manuscript reports 36
records. The remaining cases were presumably run on the HPC; their inputs, logs and
processed profiles need to be retrieved before the MD part of the archive is complete.

## Proposed archive layout (to be filled by the authors)

```
data/md/<compound>/<NaCl_mM>mM_<Gamma>/   in.lmp, system.data, log.lammps,
                                          Na.profile, Cl.profile, PFAS.profile, water.profile,
                                          rdf_*.rdf, out.charge_density, out.surfTension.txt
data/md/README.md                         case table mapping each directory to Table S1 rows
code/md/                                  LAMMPS input template, force-field files,
                                          profile -> binding-ratio and Poisson-potential scripts
```

Trajectories: total ~1.2 GB for the five local cases. Zenodo accepts up to 50 GB per record,
so they can go into the Zenodo deposit (not GitHub) if the authors want full reproducibility;
otherwise state "available from the corresponding author on request" in the Data availability
statement.

## Decisions needed from the authors

1. Which campaign (LAMMPS `binding model` vs GROMACS `PBE`) produced each of the 36 Table S1
   records, and where the missing ~31 LAMMPS cases are stored.
2. Whether to deposit trajectories.
3. Whether the Figure 1g,h cases should be rerun or the digitized rows marked as such in Table S1.
