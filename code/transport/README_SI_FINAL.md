# PFAS transport + inverse SSL — SI concentration / bounded retry v5

This is the only package for the current PFOS/PFOA comparison. It inherits the v4 bounded-retry numerical machinery but does **not** inherit its physical parameters or outputs.

| Case | AWI relation | Parameters |
|---|---|---|
| PFOS Langmuir | matched 10 mM | KL=123.794092031724 m3/mol; Gamma_m=2.49648324768916e-6 mol/m2 |
| PFOS mass action | Table S1 concentration | K1=1553 m3/mol; K2=3.8825 m6/mol2; Gamma_m=3.0695e-6 mol/m2 |
| PFOA Langmuir | matched 10 mM | KL=6.227660735 m3/mol; Gamma_m=2.988939298e-6 mol/m2 |
| PFOA mass action | Table S1 concentration | K1=49.2 m3/mol; K2=0.13 m6/mol2; Gamma_m=3.36e-6 mol/m2 |

All four use NaCl=10 mM, T=298 K, `use_activity_correction=false`, Fs=Faw=1, n=1, PFOS Kd=6.9 L/kg, PFOA Kd=2.2 L/kg, the same measured total-PFAS profile, grid, flow forcing, Aaw, boundary conditions and receptor. The only comparison change is AWI adsorption relation. Preflight aborts on any violation.

## Run and monitor

```bash
tar -xzf PFASlive.tgz
cd PFASlive
./go t
./go s
```

These are foreground runs: the terminal continuously shows case, simulated day/year, dt, retries, mass error, AF/DAF and wall time. `./go s 3` runs only SSL case 3. `./go x` performs the four-case smoke workflow. `./go m` follows the transport log. Slurm scripts remain available but are not the default workflow.

For Slurm: `bash submit_standard.sh YOUR_PI_GROUP transport` or `bash submit_standard.sh YOUR_PI_GROUP ssl`. Windfall equivalents are `bash submit_windfall.sh transport` and `bash submit_windfall.sh ssl`. Monitor with `squeue -u "$USER"`, `tail -f slurm-*.out`, and each case `progress_case_*.csv`.

Transport runs 5 forcing cycles (50 years) with monthly profiles. SSL uses the long horizon configured for each case and rejects an accepted endpoint if `Horizon_Peak_Covered` is false. Every trial is isolated; only same-version completed trials in this package’s output tree can resume.

`Groundwater_DAF=C_leachate/C_receptor` is the groundwater dilution divisor. `Groundwater_AF=C_receptor/C_leachate=1/DAF` is the groundwater mapping factor. The vadose/source transport attenuation is represented separately by the simulated leachate trajectory; it is not called DAF.

Do not use `Results_transport_50yr_v4`, legacy Figure source data, Ryan parameter outputs, or activity-corrected mass-action results for the formal SI Table comparison.
