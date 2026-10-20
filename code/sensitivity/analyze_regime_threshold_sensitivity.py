"""Quantify operational penalty-tolerance sensitivity without refitting a model.

Uses the final SI parameter table, the fixed-NaCl regime-map equations and
the archived Figure 4 curves. All paths are supplied relative to --project.
Run: python analyze_regime_threshold_sensitivity.py --project <manuscript folder>
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.optimize import brentq

R, T, F = 8.314, 298.0, 96485.0
EPSILON = 8.85e-12 * 78.4
TAUS = (0.05, 0.10, 0.20)
# Molecular weights follow the existing concentration--ionic-strength map.
PARAMETERS = {
    "PFOA": (49.2, 0.13, 3.36e-6, 414.0),
    "PFHxA": (1.69, 0.0046, 2.56e-6, 313.0),
    "PFOS": (1553.0, 3.8825, 3.0695e-6, 499.0),
}
REGIME_COLORS = ("#E4F1FA", "#EBEAF6", "#E7F3E7")
LINE_COLORS = ("#355F8D", "#202020", "#B56C36")
LINE_STYLES = (":", "-", "--")


def setup(name, salt):
    k1, k2, gm, mw = PARAMETERS[name]
    b = k2 / k1 * salt
    q = np.sqrt(8 * EPSILON * R * T * salt)
    return k1, gm, mw, b, q


def theta(phi, gm, b, q):
    return q * np.sinh(phi / 2) * (1 + b * np.exp(phi)) / (F * gm)


def phi_at_theta(value, gm, b, q):
    return brentq(lambda p: theta(p, gm, b, q) - value, 0, 30,
                  xtol=1e-14, rtol=1e-13)


def electrostatic(phi, b):
    return np.log((1 + b) / (np.exp(-phi) + b))


def concentration(phi, k1, gm, mw, b, q):
    th = theta(phi, gm, b, q)
    return th / ((1 - th) * k1 * (np.exp(-phi) + b)) * mw


def classify(le, ls, tau):
    low = (le <= tau) & (ls <= tau)
    return np.where(low, 0, np.where(le >= ls, 1, 2)).astype(np.uint8)


def fields(name, salts, concentrations):
    le = np.empty((len(concentrations), len(salts)))
    ls = np.empty_like(le)
    max_boundary_residual = 0.0
    for j, salt in enumerate(salts):
        k1, gm, mw, b, q = setup(name, salt)
        pmax = phi_at_theta(1 - 1e-10, gm, b, q)
        p = np.linspace(1e-10, pmax, 4000)
        c = concentration(p, k1, gm, mw, b, q)
        pp = np.interp(np.log10(concentrations), np.log10(c), p)
        tt = theta(pp, gm, b, q)
        le[:, j] = electrostatic(pp, b)
        ls[:, j] = -np.log1p(-tt)
        for tau in TAUS:
            onset_phi, _ = onset(name, salt, tau)
            te = theta(onset_phi, gm, b, q)
            residual = abs(max(electrostatic(onset_phi, b), -np.log1p(-te)) - tau)
            max_boundary_residual = max(max_boundary_residual, residual)
    return le, ls, max_boundary_residual


def onset(name, salt, tau):
    k1, gm, mw, b, q = setup(name, salt)
    ps = phi_at_theta(-np.expm1(-tau), gm, b, q)
    term = (1 + b) * np.exp(-tau) - b
    pe = -np.log(term) if term > 0 else np.inf
    return min(ps, pe), "electrostatic" if pe <= ps else "steric"


def equality(name, salt):
    k1, gm, mw, b, q = setup(name, salt)
    # Compare dilute-limit penalty slopes before looking for a nonzero root.
    if 2 * F * gm / (q * (1 + b)**2) <= 1:
        return np.nan, np.nan
    pmax = phi_at_theta(1 - 1e-10, gm, b, q)
    def delta(p):
        return electrostatic(p, b) + np.log1p(-theta(p, gm, b, q))
    p = brentq(delta, 1e-8, pmax, xtol=1e-13)
    return concentration(p, k1, gm, mw, b, q), electrostatic(p, b)


def terminal_salt(name, tau):
    def delta(logsalt):
        k1, gm, mw, b, q = setup(name, 10**logsalt)
        ps = phi_at_theta(-np.expm1(-tau), gm, b, q)
        return electrostatic(ps, b) - tau
    assert delta(0) > 0 and delta(3) < 0
    return 10**brentq(delta, 0, 3, xtol=1e-12)


def transfer_sensitivity(source, out):
    data = pd.read_csv(source / "fig4_a_f_model_curves.csv")
    params = pd.read_csv(source / "fig4_parameters.csv")
    for row in params.itertuples():
        assert np.allclose([row.K1_m3_mol, row.K2_m6_mol2, row.Gamma_m_MA_mol_m2],
                           PARAMETERS[row.PFAS][:3], rtol=1e-12, atol=0)
    data = data.merge(params, on=["PFAS", "I_mM"], validate="many_to_one")
    data["Lambda_st_MA"] = -np.log1p(-data.Gamma_MA_umol_m2 * 1e-6 / data.Gamma_m_MA_mol_m2)
    data["Lambda_st_SL"] = -np.log1p(-data.Gamma_SL_umol_m2 * 1e-6 / data.Gamma_m_SL_mol_m2)
    phi = -F * data.Potential_V / (R * T)
    b = data.K2_m6_mol2 / data.K1_m3_mol * data.I_mM
    le = np.log((1 + b) / (np.exp(-phi) + b))
    lfit = np.log(data.Gamma_m_MA_mol_m2 * (data.K1_m3_mol + data.K2_m6_mol2 * data.I_mM)
                  / (data.Gamma_m_SL_mol_m2 * data.K_SL_m3_mol))
    data["Delta_Lambda_el"] = lfit - le
    data["Reference_exp_neg_delta"] = np.exp(-data.Delta_Lambda_el)
    rows = []
    for tau in TAUS:
        use = (data.Lambda_st_MA <= tau) & (data.Lambda_st_SL <= tau)
        subset = data[use]
        err = subset.GammaSL_over_GammaMA - subset.Reference_exp_neg_delta
        rows.append({"Tolerance": tau, "n_eligible": int(use.sum()),
                     "R_squared_fixed_reference": 1 - np.sum(err**2) /
                     np.sum((subset.GammaSL_over_GammaMA - subset.GammaSL_over_GammaMA.mean())**2),
                     "RMSE_fixed_reference": np.sqrt(np.mean(err**2)),
                     "Maximum_abs_ratio_error": float(np.max(np.abs(err))),
                     "Maximum_relative_reference_error": float(np.max(np.abs(err / subset.Reference_exp_neg_delta)))})
        data[f"Eligible_{tau:.2f}"] = use
    summary = pd.DataFrame(rows)
    baseline = pd.read_csv(source / "selected_pfas_master_curve_linearity.csv").iloc[0]
    actual = summary[summary.Tolerance == .1].iloc[0]
    assert actual.n_eligible == baseline.n_points == 1158
    assert np.isclose(actual.R_squared_fixed_reference, baseline.R_squared, atol=1e-12)
    assert np.isclose(actual.RMSE_fixed_reference, baseline.RMSE, atol=1e-12)
    data.to_csv(out / "threshold_transfer_pointwise.csv", index=False)
    summary.to_csv(out / "threshold_transfer_summary.csv", index=False)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--fig4-source", type=Path,
                        help="Archived Figure 4 CSV directory; defaults to project/figures/source_data/fig4")
    parser.add_argument("--make-figure", action="store_true")
    args = parser.parse_args()
    project = args.project.resolve()
    source = args.fig4_source.resolve() if args.fig4_source else project / "figures/source_data/fig4"
    out = project / "SI/source_data/regime_threshold_sensitivity"
    out.mkdir(parents=True, exist_ok=True)
    salts = np.logspace(0, 3, 600)
    concs = np.logspace(-4, 4, 700)
    boundaries, selected, summaries, plots = [], [], [], []
    max_residual = 0.0
    for name in PARAMETERS:
        le, ls, residual = fields(name, salts, concs)
        max_residual = max(max_residual, residual)
        baseline = classify(le, ls, .1)
        assert np.all((classify(le, ls, .05) == 0) <= (baseline == 0))
        assert np.all((baseline == 0) <= (classify(le, ls, .2) == 0))
        equality_values = [equality(name, salt) for salt in salts]
        ceq = np.array([x[0] for x in equality_values])
        for tau in TAUS:
            codes = classify(le, ls, tau)
            cutoffs = []
            for j, salt in enumerate(salts):
                k1, gm, mw, b, q = setup(name, salt)
                p, limiting = onset(name, salt, tau)
                c = concentration(p, k1, gm, mw, b, q)
                cutoffs.append(c)
                boundaries.append({"PFAS": name, "Tolerance": tau, "NaCl_mM": salt,
                                   "Surface_activity_exit_mg_L": c,
                                   "Limiting_penalty": limiting,
                                   "Equal_penalty_concentration_mg_L": ceq[j],
                                   "Penalty_at_equality": equality_values[j][1]})
            crit = terminal_salt(name, tau)
            summary = {"PFAS": name, "Tolerance": tau,
                       "Single_term_attenuation_percent": 100 * -np.expm1(-tau),
                       "Maximum_combined_attenuation_percent": 100 * -np.expm1(-2 * tau),
                       "Maximum_coverage_in_surface_activity_regime": -np.expm1(-tau),
                       "Reclassified_log_grid_percent_vs_0p10": 100 * np.mean(codes != baseline),
                       "Terminal_salt_for_distinct_electrostatic_regime_mM": crit,
                       "Grid_concentration_count": len(concs), "Grid_salt_count": len(salts)}
            for code, label in enumerate(["Surface_activity", "Electrostatic", "Steric"]):
                summary[f"{label}_log_grid_percent"] = 100 * np.mean(codes == code)
            for salt in (1, 10, 100, 200, 300, 500, 1000):
                k1, gm, mw, b, q = setup(name, salt)
                p, limiting = onset(name, salt, tau)
                selected.append({"PFAS": name, "Tolerance": tau, "NaCl_mM": salt,
                                 "Surface_activity_exit_mg_L": concentration(p, k1, gm, mw, b, q),
                                 "Limiting_penalty": limiting})
            summaries.append(summary)
        plots.append((name, baseline, ceq))
    frame = pd.DataFrame(boundaries)
    frame.to_csv(out / "threshold_regime_boundaries.csv", index=False)
    pd.DataFrame(selected).to_csv(out / "threshold_selected_concentrations.csv", index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(out / "threshold_regime_summary.csv", index=False)
    # Same deterministic grid under all tolerances. No experimental replication.
    transfer = transfer_sensitivity(source, out)
    diagnostics = {"maximum_absolute_onset_penalty_residual": max_residual,
                   "Figure_4g_baseline_reproduced": True,
                   "fixed_parameters": PARAMETERS,
                   "concentration_domain_mg_L": [float(concs[0]), float(concs[-1])],
                   "NaCl_domain_mM": [float(salts[0]), float(salts[-1])],
                   "grid_shape": [len(concs), len(salts)],
                   "equations": "Fixed-background-NaCl regime map; Figure 4 curves reused without refitting",
                   "input_SHA256": {n: hashlib.sha256((source/n).read_bytes()).hexdigest() for n in
                       ["fig4_a_f_model_curves.csv", "fig4_parameters.csv", "selected_pfas_master_curve_linearity.csv"]}}
    assert max_residual < 1e-10
    (out / "verification.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    if args.make_figure:
        make_figure(project, salts, concs, plots, frame)
    print(summary.to_string(index=False))
    print(transfer.to_string(index=False))


def make_figure(project, salts, concs, plots, boundaries):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 7, "axes.labelsize": 8, "axes.titlesize": 9,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": .65})
    fig, axes = plt.subplots(1, 3, figsize=(183/25.4, 98/25.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=.085, right=.985, bottom=.32, top=.91, wspace=.18)
    cmap = ListedColormap(REGIME_COLORS)
    norm = BoundaryNorm(np.arange(-.5, 3.5), 3)
    for index, (ax, (name, codes, ceq)) in enumerate(zip(axes, plots)):
        ax.pcolormesh(salts, concs, codes, cmap=cmap, norm=norm, shading="auto", rasterized=True)
        for tau, color, ls in zip(TAUS, LINE_COLORS, LINE_STYLES):
            part = boundaries[(boundaries.PFAS == name) & (boundaries.Tolerance == tau)]
            ax.plot(salts, part.Surface_activity_exit_mg_L, color=color, linestyle=ls, lw=1.1)
        ax.plot(salts, ceq, color="#777777", linestyle="-.", lw=.95)
        ax.set(xscale="log", yscale="log", xlim=(1, 1000), ylim=(1e-4, 1e4), title=name)
        ax.set_xticks([1,10,100,1000])
        ax.set_yticks([1e-4,1e-2,1,1e2,1e4])
        ax.tick_params(which="minor", length=0)
        ax.tick_params(which="major", direction="out", length=3, width=.65)
        ax.set_xlabel("NaCl concentration (mM)")
        ax.text(-.16 if index == 0 else -.07, 1.05, chr(97+index), transform=ax.transAxes,
                fontsize=10, fontweight="bold")
    axes[0].set_ylabel(r"PFAS concentration (mg L$^{-1}$)")
    handles = [Line2D([], [], color=c, linestyle=s, lw=1.1, label=rf"$\Lambda_c={t:.2f}$")
               for t,c,s in zip(TAUS,LINE_COLORS,LINE_STYLES)]
    handles.append(Line2D([], [], color="#777777", ls="-.", lw=.95,
                          label=r"$\Lambda_{\rm el}=\Lambda_{\rm st}$"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.52,.16), ncol=4,
               frameon=False, handlelength=2.8, columnspacing=1.7)
    patches = [Patch(facecolor=c, label=label) for c,label in zip(REGIME_COLORS,
        ["Surface-activity control", "Electrostatic inhibition", "Steric exclusion"])]
    fig.legend(handles=patches, loc="lower center", bbox_to_anchor=(.52,.065), ncol=3,
               title=r"Background classification at $\Lambda_c=0.10$", title_fontsize=7,
               frameon=False, columnspacing=1.6)
    output = project / "SI/figures/regime_threshold_sensitivity"
    output.parent.mkdir(parents=True,exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(output.with_suffix("."+ext), dpi=400 if ext=="png" else 300)
    plt.close(fig)


if __name__ == "__main__":
    main()
