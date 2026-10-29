# Molecular dynamics: representative case

One complete LAMMPS case is deposited so that the simulation protocol described in
Supplementary Section S1 can be inspected and rerun. The other cases of the campaign use the
same input script with a different initial configuration (`Data.data`: number of PFAS
molecules, number of Na+/Cl- ions, compound-specific force-field block). Trajectories are not
deposited because of their size; they are available from the corresponding author on request.

## `PFOS_NaCl100mM_Gamma0.26/`

Composition of `Data.data` (counted from the atom list): 72 PFOS- (36 per interface),
135 Na+, 63 Cl-, 75,150 SPC/E water molecules; box 151 x 151 x 220 A with periodic
boundaries and a water slab exposing two air-water interfaces. With the 151 x 151 A
cross-section, 36 PFOS per interface corresponds to a surface excess of 0.26 umol m-2.
The directory label follows the authors' case naming ("I = 100, middle surface coverage");
the bulk NaCl concentration reported in Supplementary Table S1 is the value evaluated from
the ion profiles in the bulk region of the film, not the nominal count.

| File | Content |
|---|---|
| `in.lmp` | LAMMPS input script (originally `test.lmp`): OPLS-AA PFOS, SPC/E water with SHAKE, Dang-Smith Na+/Cl-, `lj/cut/coul/long 10.0` with PPPM 1e-4, minimization, NVT 300 K (Nose-Hoover, 100 fs damping), 1 fs step, `run 2000000`; defines the chunk/atom profiles, RDFs, charge-density and surface-tension outputs listed below |
| `Data.data` | Initial configuration and topology (LAMMPS `atom_style full`) |
| `log.lammps` | Run log (LAMMPS 27 Oct 2021 build); thermodynamic output every 10,000 steps |
| `Na.profile`, `Cl.profile`, `PFAS.profile`, `O.profile`, `water.profile` | Number-density profiles along z (1 A bins, 220 chunks), time-averaged over 100,000-step windows (`fix ave/chunk`); `O.profile` is the PFOS sulfonate oxygen |
| `out.charge_density` | Charge-density profile along z (same binning), input to the one-dimensional Poisson solution for the interfacial potential (SI S1.3) |
| `rdf_SNa.rdf`, `rdf_ONa.rdf`, `rdf_SCl.rdf`, `rdf_OCl.rdf` | Radial distribution functions between the sulfonate S or O and Na+/Cl- (200 bins, 10 A) |
| `NaPFAS.txt`, `ClPFAS.txt`, `PFASPFAS.txt`, `waterPFAS.txt` | Time-averaged group/group interaction energies with PFOS (kcal/mol) |
| `out.surfTension.txt` | Instantaneous surface tension from the pressure tensor, 0.5 lz (Pzz - (Pxx + Pyy)/2), written every 50 steps (mN/m, two interfaces) |

Not deposited from this case directory: `trajectory.data` (atom dump, 298 MB), `2ns.nvt`
and `minimized.data` (restart/intermediate configurations), `charge_xy.dat` (2-D charge
map, 12.6 MB).

## Post-processing

The counterion binding ratio is the net countercharge (Na+ minus Cl-) within the Stern
layer divided by the number of adsorbed PFAS, evaluated from the profiles above after
reflecting and averaging the two interfaces; the interfacial potential is obtained from
`out.charge_density` by integrating the one-dimensional Poisson equation (Supplementary
Section S1). The MATLAB live scripts used for this step are not yet included in the
repository (see `docs/MD_DATA_INVENTORY.md`).

## Rerunning

```text
lmp -in in.lmp
```

A run of 2,000,000 steps with 227,736 atoms takes on the order of days on a single node;
the archived log reports the wall time of the original run.
