"""Evaluate the declared manuscript parameters without fitting or editing inputs.

Run from the Nature_Water paper workspace. Surface-tension observations and
the historical rendered curves are read from the existing Fig. 2 source table.
Binding observations are read from the archived numerical MD records, retain
the pre-existing plotting exclusion, and count the exact PFHxA duplicate once.
The current parameter sets come from model_parameter_provenance.csv.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq


PROJECT = Path(__file__).resolve().parents[2]
PAPER = PROJECT.parents[2]
DATA = PROJECT / "SI/source_data"
OUTPUT = DATA / "model_evaluation_20261028"
AUDIT = PROJECT.parent / "model_evaluation_20261028"
SOURCE = PAPER / "NatureWater_Overleaf/figures/fig3_surface_tension_source_data.csv"
MD_SOURCE = DATA / "md_counterion_binding_source_records.csv"
PARAM_SOURCE = DATA / "model_parameter_provenance.csv"
F, R, T, EPSILON = 96485.0, 8.314, 298.0, 8.85e-12 * 78.4
SIGMA0 = 72.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def state(c: float, salt: float, k1: float, k2: float, gm: float) -> dict:
    """Solve the stated MA/GC equations; c and salt are both in mol/m3."""
    total = salt + c  # Includes the counterions accompanying dissolved PFAS.
    diffuse = np.sqrt(8 * EPSILON * R * T * total)

    def populations(phi: float) -> tuple[float, float]:
        activity = k1 * np.exp(phi) + k2 * total
        occupancy = c * activity
        return gm * occupancy / (1 + occupancy), k2 * total / activity

    def charge_residual(phi: float) -> float:
        gamma, beta = populations(phi)
        return F * gamma * (1 - beta) + diffuse * np.sinh(phi / 2)

    phi = brentq(charge_residual, -30, 0, xtol=1e-14, rtol=1e-14)
    gamma, beta = populations(phi)
    pressure = 1000 * (
        -R * T * gm * np.log1p(-gamma / gm)
        + 4 * R * T / F * np.sqrt(2 * EPSILON * R * T * total)
        * (np.cosh(phi / 2) - 1)
    )
    charge_error = abs(charge_residual(phi)) / max(
        F * gamma * (1 - beta), np.finfo(float).tiny
    )
    return {"sigma_predicted_mN_m": SIGMA0 - pressure,
            "Gamma_predicted_mol_m2": gamma, "phi_predicted": phi,
            "beta_predicted": beta, "charge_relative_residual": charge_error}


def summarize(frame: pd.DataFrame, error: str) -> dict:
    residuals = frame[error].to_numpy()
    return {"N": len(frame), "SSE": float(residuals @ residuals),
            "RMSE": float(np.sqrt(np.mean(residuals**2))),
            "mean_signed_error": float(np.mean(residuals)),
            "maximum_absolute_error": float(np.max(np.abs(residuals)))}


def main() -> None:
    OUTPUT.mkdir(exist_ok=True)
    (AUDIT / "before").mkdir(parents=True, exist_ok=True)
    for name in ("main.tex", "supporting_information.tex"):
        backup = AUDIT / "before" / name
        if not backup.exists():
            shutil.copy2(PROJECT / name, backup)
    canonical = pd.read_csv(PARAM_SOURCE)
    canonical = canonical[canonical.Purpose.eq(
        "Final manuscript mass-action parameter set")].set_index("Compound")
    assert set(canonical.index) == {"PFOA", "PFHxA", "PFOS"}
    canonical.to_csv(OUTPUT / "parameter_sets_used.csv")
    source = pd.read_csv(SOURCE)
    assert len(source) == 4951
    observations = source[source.Series.eq("Surface tension experiment")].copy()
    pointwise = []
    for row in observations.itertuples():
        pars = canonical.loc[row.Compound]
        result = state(row.X, row.IonicStrength_mM, pars.K1_m3_mol,
                       pars.K2_m6_mol2, pars.Gamma_m_mol_m2)
        role = ("prediction" if row.IonicStrength_mM == 0
                or (row.Compound == "PFOA" and row.IonicStrength_mM == 10)
                else "calibration-condition comparison")
        # PFOS 50 mM comes from an excluded/comparison file in the original
        # exporter. Do not silently designate it as an independently verified
        # calibration or held-out set; no PFOS aggregate calibration metric is used.
        if row.Compound == "PFOS" and row.IonicStrength_mM == 50:
            role = "PFOS comparison; original exporter reads excluded-data file"
        if row.Compound == "PFOA" and row.IonicStrength_mM == 10:
            reference = "Lyu et al. 2018; as designated in the manuscript"
        elif row.Compound in {"PFOA", "PFHxA"}:
            reference = "Le et al. 2021; attribution follows the manuscript"
        else:
            reference = "Screened PFOS compilation in original figure source"
        pointwise.append({"PFAS": row.Compound, "NaCl_mM": row.IonicStrength_mM,
                          "c_mol_m3": row.X, "sigma_observed_mN_m": row.Y,
                          "role": role, "reference": reference, **result})
    st = pd.DataFrame(pointwise)
    st["residual_mN_m"] = st.sigma_predicted_mN_m - st.sigma_observed_mN_m
    st.to_csv(OUTPUT / "surface_tension_pointwise.csv", index=False)
    summaries = []
    for (name, salt), frame in st.groupby(["PFAS", "NaCl_mM"], sort=False):
        summaries.append({"PFAS": name, "dataset": f"{salt:g} mM NaCl",
                          "role": frame.role.iloc[0],
                          **summarize(frame, "residual_mN_m")})
    for name, salts in (("PFOA", [50, 100, 200, 300, 500]),
                        ("PFHxA", [100, 200, 300, 500])):
        frame = st[st.PFAS.eq(name) & st.NaCl_mM.isin(salts)]
        summaries.append({"PFAS": name, "dataset": "pooled calibration conditions",
                          "role": "calibration-condition comparison",
                          **summarize(frame, "residual_mN_m")})
    frame = st[st.PFAS.eq("PFOA") & st.NaCl_mM.isin([0, 10])]
    summaries.append({"PFAS": "PFOA", "dataset": "pooled 0 and 10 mM predictions",
                      "role": "prediction", **summarize(frame, "residual_mN_m")})
    st_summary = pd.DataFrame(summaries)
    st_summary.to_csv(OUTPUT / "surface_tension_summary.csv", index=False)

    records = pd.read_csv(MD_SOURCE)
    plotted = records[records.status.eq("plotted")].copy()
    assert len(plotted) == 37
    distinct_columns = ["PFAS", "NaCl_mM", "Gamma_umol_m2", "c2s_mM",
                        "c2s_reported_error_mM", "beta_percent",
                        "beta_reported_error_percentage_points"]
    md = plotted.drop_duplicates(distinct_columns).copy()
    assert len(md) == 36
    for index, row in md.iterrows():
        pars = canonical.loc[row.PFAS]
        kb = pars.K2_m6_mol2 / pars.K1_m3_mol
        md.loc[index, "Kbind_m3_mol"] = kb
        md.loc[index, "beta_predicted_percent"] = 100 * kb * row.c2s_mM / (1 + kb * row.c2s_mM)
    md["residual_percentage_points"] = md.beta_predicted_percent - md.beta_percent
    md.to_csv(OUTPUT / "binding_pointwise.csv", index=False)
    binding_summaries = []
    for name, frame in md.groupby("PFAS", sort=False):
        binding_summaries.append({"PFAS": name, "Kbind_m3_mol": frame.Kbind_m3_mol.iloc[0],
                                  **summarize(frame, "residual_percentage_points")})
    binding_summaries.append({"PFAS": "pooled", "Kbind_m3_mol": np.nan,
                              **summarize(md, "residual_percentage_points")})
    binding_summary = pd.DataFrame(binding_summaries)
    binding_summary.to_csv(OUTPUT / "binding_summary.csv", index=False)

    # Check the MA state using the independently parameterized implicit-coverage
    # solver already used in SI Fig. S3. No fitting or plotting is performed.
    spec = importlib.util.spec_from_file_location(
        "slope_analysis", PROJECT / "SI/code/plot_adsorption_slopes.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    gamma_errors = []
    from types import SimpleNamespace
    for row in st.itertuples():
        pars = canonical.loc[row.PFAS]
        alternate, _ = module.state(row.c_mol_m3, SimpleNamespace(
            K1_m3_mol=pars.K1_m3_mol, K2_m6_mol2=pars.K2_m6_mol2,
            Gamma_m_MA_mol_m2=pars.Gamma_m_mol_m2, I_mM=row.NaCl_mM))
        gamma_errors.append(abs(alternate / row.Gamma_predicted_mol_m2 - 1))
    assert max(gamma_errors) < 1e-8
    assert st.charge_relative_residual.max() < 1e-7

    # Diagnose the current exported curves separately. Do not overwrite them,
    # infer new best-fit parameters, or report their interpolated errors as
    # errors of the current declared parameter sets.
    archive_metrics, curve_differences, binding_archive = [], [], []
    for (name, salt), frame in observations.groupby(["Compound", "IonicStrength_mM"]):
        curve = source[source.Compound.eq(name) & source.Series.eq("Surface tension model")
                       & source.IonicStrength_mM.eq(salt)].sort_values("X")
        interpolated = PchipInterpolator(np.log(curve.X), curve.Y)(np.log(frame.X))
        residual = interpolated - frame.Y.to_numpy()
        archive_metrics.append({"PFAS": name, "NaCl_mM": salt, "N": len(frame),
                                "archive_interpolated_RMSE_mN_m": float(np.sqrt(np.mean(residual**2))),
                                "method": "PCHIP in log concentration; diagnostic only"})
        pars = canonical.loc[name]
        predicted = np.array([state(x, salt, pars.K1_m3_mol, pars.K2_m6_mol2,
                                    pars.Gamma_m_mol_m2)["sigma_predicted_mN_m"] for x in curve.X])
        curve_differences.append({"PFAS": name, "NaCl_mM": salt,
                                  "max_abs_current_minus_exported_mN_m": float(np.max(np.abs(predicted - curve.Y)))})
    for name, frame in md.groupby("PFAS", sort=False):
        curve = source[source.Compound.eq(name) & source.Series.eq("Binding master curve")
                       & source.X.gt(0) & source.Y.gt(0) & source.Y.lt(100)]
        fraction = curve.Y / 100
        archived_kb = float(np.median(fraction / (1 - fraction) / curve.X))
        current_kb = float(frame.Kbind_m3_mol.iloc[0])
        predicted = 100 * archived_kb * frame.c2s_mM / (1 + archived_kb * frame.c2s_mM)
        binding_archive.append({"PFAS": name, "exported_curve_Kbind_m3_mol": archived_kb,
                                "current_declared_Kbind_m3_mol": current_kb,
                                "exported_curve_RMSE_percentage_points": float(np.sqrt(np.mean((predicted-frame.beta_percent)**2)))})
    pd.DataFrame(archive_metrics).to_csv(OUTPUT / "exported_surface_curve_diagnostics.csv", index=False)
    pd.DataFrame(curve_differences).to_csv(OUTPUT / "current_vs_exported_surface_curves.csv", index=False)
    pd.DataFrame(binding_archive).to_csv(OUTPUT / "current_vs_exported_binding_curves.csv", index=False)
    verification = {"evaluation_date": "2026-10-28", "parameter_basis": "current declared parameter sets; no refitting",
                    "baseline_surface_tension_mN_m": SIGMA0, "temperature_K": T,
                    "counterion_bulk_convention": "c2,0 = NaCl concentration + dissolved PFAS concentration, in mol/m3",
                    "error_definition": "sqrt(sum((prediction-observation)^2)/N), unweighted",
                    "surface_observations": len(st), "binding_plotted_rows": len(plotted),
                    "binding_distinct_rows": len(md), "binding_preexisting_excluded_rows": int(records.status.ne("plotted").sum()),
                    "max_charge_relative_residual": float(st.charge_relative_residual.max()),
                    "max_independent_solver_Gamma_relative_difference": float(max(gamma_errors)),
                    "source_hashes": {str(p.relative_to(PAPER)): sha256(p) for p in [SOURCE, MD_SOURCE, PARAM_SOURCE]},
                    "figures_modified": False,
                    "caution": "The current rendered Fig. 2 source curves are not all generated with the declared parameter sets. Their diagnostics are archived separately and must not be confused with the reported current-parameter metrics."}
    (OUTPUT / "verification.json").write_text(json.dumps(verification, indent=2), encoding="utf-8")
    print(st_summary.to_string(index=False))
    print(binding_summary.to_string(index=False))
    print(json.dumps(verification, indent=2))


if __name__ == "__main__":
    main()
