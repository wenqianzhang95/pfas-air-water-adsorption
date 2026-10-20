"""Compare full-model and Langmuir surface excess for PFOA, PFHxA, and PFOS.

Regime backgrounds follow the adsorption-penalty definition:
1) both electrostatic and steric penalties are <= 0.1;
2) electrostatic inhibition is the larger penalty;
3) steric exclusion is the larger penalty.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq, least_squares
from matplotlib.ticker import FuncFormatter, MaxNLocator


ROOT = Path(__file__).resolve().parent
F = 96485.0
R = 8.314
T = 298.0
EPSILON = 8.85e-12 * 78.4
PENALTY_THRESHOLD = 0.1

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.8,
        "text.color": "#202020",
        "axes.labelcolor": "#202020",
        "axes.edgecolor": "#202020",
        "xtick.color": "#202020",
        "ytick.color": "#202020",
        "axes.spines.top": False,
        "legend.frameon": False,
    }
)

REGIME_COLORS = ("#edf5fa", "#f2f1f7", "#eff7f0")
ION_COLORS = ("#802703", "#fc8c3b", "#fddebe")
EXPERIMENT_COLOR = "#767676"


def format_log_tick(value: float, _: int) -> str:
    """Render log ticks through mathtext so Arial is used consistently."""
    if value <= 0:
        return ""
    exponent = int(round(np.log10(value)))
    if np.isclose(value, 10.0 ** exponent, rtol=1e-10):
        return rf"$10^{{{exponent}}}$"
    return ""


LOG_TICK_FORMATTER = FuncFormatter(format_log_tick)


@dataclass(frozen=True)
class Compound:
    name: str
    molecular_weight: float
    ionic_strengths: tuple[float, ...]
    k1: tuple[float, ...]
    k2: tuple[float, ...]
    gamma_m_full: tuple[float, ...]
    k_langmuir: tuple[float, ...]
    gamma_m_langmuir: tuple[float, ...]
    c_range: tuple[float, float]
    x_limits: tuple[float, float]


PFOA = Compound(
    "PFOA", 413.0, (10.0, 50.0, 100.0),
    (49.2,) * 3, (0.13,) * 3, (3.36e-6,) * 3,
    (6.227663, 20.634773, 25.406345),
    (2.988939e-6, 2.937332e-6, 3.216727e-6),
    (1e-7, 10 ** 1.5), (1e-5, 1e4),
)
PFHXA = Compound(
    "PFHxA", 313.0, (100.0, 200.0, 300.0),
    (1.69,) * 3, (0.0046,) * 3, (2.56e-6,) * 3,
    (0.678870, 1.259270, 1.671099),
    (3.002462e-6, 2.869480e-6, 2.869755e-6),
    (1e-7, 1e2), (1e-3, 1e4),
)

def solve_full_curve(c: np.ndarray, ionic_strength: float, k1: float, k2: float, gamma_m: float):
    """Solve the Kralchevsky/Stern equations in theta, binding ratio, and phi."""
    theta = np.empty_like(c)
    binding_ratio = np.empty_like(c)
    phi = np.empty_like(c)
    previous = np.array([1e-7, 0.01, 1e-5])

    for index, concentration in enumerate(c):
        total_ionic_strength = ionic_strength + concentration

        def residual(values: np.ndarray) -> np.ndarray:
            coverage, ratio, phi_value = values
            adsorption = k1 * concentration * np.exp(-phi_value) + k2 * total_ionic_strength * concentration
            target_ratio = k2 * total_ionic_strength / (k1 * np.exp(-phi_value) + k2 * total_ionic_strength)
            charge_balance = (
                F * gamma_m * coverage * (ratio - 1.0)
                + np.sqrt(8.0 * EPSILON * R * T * total_ionic_strength) * np.sinh(phi_value / 2.0)
            ) / (F * gamma_m)
            return np.array([
                np.log(coverage / (1.0 - coverage)) - np.log(max(adsorption, 1e-300)),
                ratio - target_ratio,
                charge_balance,
            ])

        if index == 0:
            adsorption_guess = k1 * concentration + k2 * total_ionic_strength * concentration
            previous = np.array([
                np.clip(adsorption_guess / (1.0 + adsorption_guess), 1e-12, 1e-3),
                np.clip(k2 * total_ionic_strength / (k1 + k2 * total_ionic_strength), 1e-9, 0.99),
                1e-8,
            ])
        result = least_squares(
            residual,
            previous,
            bounds=([1e-14, 1e-12, 1e-12], [0.999999, 0.999999, 25.0]),
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
            max_nfev=500,
        )
        if not result.success or np.linalg.norm(result.fun) > 1e-6:
            raise RuntimeError(f"Full-model solve failed at I={ionic_strength:g} mM, c={concentration:.3g} mM.")
        previous = result.x
        theta[index], binding_ratio[index], phi[index] = result.x

    return gamma_m * theta, binding_ratio, -R * T * phi / F, theta, phi


def classify_regimes(theta: np.ndarray, phi: np.ndarray, ionic_strength: float, k1: float, k2: float) -> np.ndarray:
    """Classify regimes with the revised two-penalty definition."""
    electrostatic_penalty = calculate_electrostatic_penalty(phi, ionic_strength, k1, k2)
    steric_penalty = -np.log1p(-theta)
    affinity = (electrostatic_penalty <= PENALTY_THRESHOLD) & (steric_penalty <= PENALTY_THRESHOLD)
    codes = np.zeros(theta.size, dtype=int)
    codes[~affinity & (electrostatic_penalty >= steric_penalty)] = 1
    codes[~affinity & (steric_penalty > electrostatic_penalty)] = 2
    return codes


def calculate_electrostatic_penalty(
    phi: np.ndarray, ionic_strength: float, k1: float, k2: float
) -> np.ndarray:
    """Return the dimensionless electrostatic adsorption penalty Lambda_el."""
    binding_term = (k2 / k1) * ionic_strength
    return np.log((1.0 + binding_term) / (np.exp(-phi) + binding_term))


def regime_background(ax: plt.Axes, x: np.ndarray, codes: np.ndarray) -> None:
    """Draw contiguous regime spans, extending them to the visible x limits."""
    edges = np.empty(x.size + 1)
    edges[0], edges[-1] = ax.get_xlim()
    edges[1:-1] = np.sqrt(x[:-1] * x[1:])
    start = 0
    while start < codes.size:
        stop = start + 1
        while stop < codes.size and codes[stop] == codes[start]:
            stop += 1
        left = max(edges[start], ax.get_xlim()[0])
        right = min(edges[stop], ax.get_xlim()[1])
        if right > left:
            ax.axvspan(left, right, color=REGIME_COLORS[codes[start]], ec="none", zorder=0)
            if start > 0 and ax.get_xlim()[0] < left < ax.get_xlim()[1]:
                ax.axvline(left, color="#c9c9c9", lw=0.6, zorder=1)
        start = stop


def format_log_axes(ax: plt.Axes, compound: Compound, show_y_label: bool = True) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*compound.x_limits)
    ax.set_ylim(1e-5, 1e1)
    ax.set_yticks([1e-5, 1e-3, 1e-1, 1e1])
    ax.xaxis.set_major_formatter(LOG_TICK_FORMATTER)
    ax.yaxis.set_major_formatter(LOG_TICK_FORMATTER)
    ax.tick_params(which="major", direction="out", width=0.9, length=3.5, labelsize=12)
    ax.tick_params(which="minor", length=0)
    ax.set_xlabel(rf"{compound.name} concentration (mg L$^{{-1}}$)", fontsize=12)
    if show_y_label:
        ax.set_ylabel("Surface excess\n($\\mu$mol m$^{-2}$)", fontsize=12)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)


def add_panel_label(
    ax: plt.Axes, label: str, x_position: float = -0.20, y_position: float = 1.035
) -> None:
    ax.text(
        x_position, y_position, label, transform=ax.transAxes, ha="left", va="bottom",
        fontsize=13, fontweight="bold",
    )


def longest_monotonic_segment(
    x: np.ndarray, y: np.ndarray, eligible: np.ndarray, increasing: bool = True
) -> np.ndarray:
    """Return the eligible contiguous segment with the largest selected-slope span."""
    candidates = []
    eligible_indices = np.flatnonzero(eligible)
    for block in np.split(eligible_indices, np.where(np.diff(eligible_indices) > 1)[0] + 1):
        if block.size < 2:
            continue
        products = np.diff(x[block]) * np.diff(y[block])
        tolerance = max(float(np.max(np.abs(products))) * 1e-10, 1e-14)
        breaks = products < -tolerance if increasing else products > tolerance
        for segment in np.split(block, np.where(breaks)[0] + 1):
            if segment.size >= 2:
                candidates.append(segment)
    if not candidates:
        raise RuntimeError("No monotonic prediction segment satisfies the selected regime filter.")
    return max(candidates, key=lambda segment: (np.ptp(y[segment]), segment.size))


def plot_isotherm_panel(
    ax: plt.Axes,
    letter: str,
    compound: Compound,
    ionic_strength: float,
    gamma_full: np.ndarray,
    gamma_langmuir: np.ndarray,
    x: np.ndarray,
    regimes: np.ndarray,
    show_legend: bool,
    show_y_label: bool,
    scatter_data: tuple[np.ndarray, np.ndarray, str] | None = None,
) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*compound.x_limits)
    ax.set_ylim(1e-5, 1e1)
    regime_background(ax, x, regimes)
    h_full, = ax.plot(
        x, gamma_full * 1e6, color="black", lw=1.8, zorder=3, label="Mass action model"
    )
    h_lang, = ax.plot(x, gamma_langmuir * 1e6, color="#59ad71", lw=1.8, ls=":", zorder=3, label="Langmuir isotherm")
    handles = [h_full, h_lang]
    if scatter_data is not None:
        h_data = ax.scatter(
            scatter_data[0], scatter_data[1] * 1e6, s=28, color=EXPERIMENT_COLOR, zorder=4, label=scatter_data[2]
        )
        handles.append(h_data)
    format_log_axes(ax, compound, show_y_label=show_y_label)
    ax.set_title(f"I = {ionic_strength:g} mM", fontsize=12, pad=8)
    add_panel_label(ax, letter, x_position=-0.20, y_position=1.13)
    # Shared legends are positioned in the figure header to preserve the 12 pt
    # minimum text size without covering curves or data points.


def plot_ratio_panel(
    ax: plt.Axes,
    letter: str,
    compound: Compound,
    records: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    lambda_fit: tuple[float, float, float],
) -> None:
    all_penalty_shifts = []
    all_log10_ratios = []
    ionic_handles = []
    for index, record in enumerate(records):
        _, gamma_full, gamma_langmuir, _, psi, regimes = record
        phi = -F * psi / (R * T)
        lambda_el = calculate_electrostatic_penalty(
            phi, compound.ionic_strengths[index], compound.k1[index], compound.k2[index]
        )
        penalty_shift = lambda_fit[index] - lambda_el
        log10_ratio = np.log10(gamma_langmuir / gamma_full)
        theta_full = np.clip(gamma_full / compound.gamma_m_full[index], 0.0, 1.0 - 1e-12)
        steric_penalty = -np.log1p(-theta_full)
        eligible = (steric_penalty <= PENALTY_THRESHOLD) & (penalty_shift >= 0.0)
        selected = longest_monotonic_segment(penalty_shift, log10_ratio, eligible, increasing=False)
        all_penalty_shifts.append(penalty_shift[selected])
        all_log10_ratios.append(log10_ratio[selected])
        sample = np.unique(np.r_[selected[::4], selected[-1]])
        if sample.size:
            for regime_code in (0, 1):
                regime_points = sample[regimes[sample] == regime_code]
                if regime_points.size:
                    ax.scatter(
                        penalty_shift[regime_points], log10_ratio[regime_points], s=20,
                        marker="o", facecolors="none",
                        edgecolors=ION_COLORS[index], linewidths=1.2, zorder=2,
                    )
            ionic_handles.append(
                mpl.lines.Line2D(
                    [], [], color=ION_COLORS[index], marker="o", markerfacecolor="none",
                    markeredgecolor=ION_COLORS[index], markeredgewidth=1.2, ls="", markersize=5.5,
                    label=f"I = {compound.ionic_strengths[index]:g} mM",
                )
            )

    shift_max = max(values.max() for values in all_penalty_shifts)
    log10_ratio_min = min(values.min() for values in all_log10_ratios)
    tick_locator = MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10])
    x_ticks = tick_locator.tick_values(0.0, shift_max)
    x_upper = float(x_ticks[-1])
    y_magnitude = max(-log10_ratio_min, x_upper / np.log(10.0))
    y_magnitude_ticks = tick_locator.tick_values(0.0, y_magnitude)
    y_ticks = -y_magnitude_ticks[::-1]
    ax.set_xlim(0.0, x_upper)
    ax.set_ylim(float(y_ticks[0]), 0.0)
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    master_x = np.linspace(0.0, x_upper, 200)
    ax.plot(master_x, -master_x / np.log(10.0), color="#666666", lw=1.0, ls=(0, (3, 2)), zorder=4)
    ax.set_xlabel(r"$\Delta\Lambda_{\mathrm{el}}$", fontsize=12)
    ax.set_ylabel(r"lg($\Gamma_{\mathrm{L}}/\Gamma_{\mathrm{m}}$)", fontsize=14, labelpad=2)
    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
    ax.tick_params(which="major", direction="out", width=0.9, length=3.5, labelsize=12)
    ax.tick_params(which="minor", length=0)
    ax.legend(
        handles=ionic_handles, loc="lower left", fontsize=11, frameon=False,
        handlelength=1.3, handletextpad=0.4, labelspacing=0.25,
        borderpad=0.2, borderaxespad=0.3,
    )
    add_panel_label(ax, letter, x_position=-0.16, y_position=1.07)


def pfoa_experimental_data(ionic_strength: float) -> tuple[np.ndarray, np.ndarray, str]:
    """Return the recorded Buckley PFOA surface-excess measurements."""
    datasets = {
        10.0: np.array([
            [8.309558516304949e-7, 1.1982644441374595e-10],
            [2.6072527543749535e-6, 2.085119479272611e-10],
            [5.235059376744198e-6, 3.2488986712388876e-10],
            [7.614389482996116e-6, 1.3208155685780405e-9],
            [1.8843209918922845e-5, 2.499050488724877e-9],
            [3.2928393745285405e-5, 2.5318768529687156e-9],
            [3.433115486031564e-5, 2.7902783461980338e-9],
            [6.253474540786146e-5, 3.6298208124593774e-9],
            [7.492330415796452e-5, 5.66015901205234e-9],
            [1.730136777597597e-4, 7.360560542480119e-9],
        ]),
        100.0: np.array([
            [1.2110883820290888e-7, 1.2704188950706232e-10],
            [2.3018827031543657e-7, 1.0594852066878877e-10],
            [7.412722870049782e-7, 6.526430987054734e-10],
            [1.0657253437853875e-6, 5.676832964384601e-10],
            [1.8087647911371963e-6, 1.3421990278306638e-9],
            [2.109531587161151e-6, 1.0894778930807533e-9],
            [2.7079292080660584e-6, 2.8799014132938425e-9],
            [4.23310864134842e-6, 2.5397229048151797e-9],
            [5.216909618968103e-6, 3.2152411727962222e-9],
            [6.610210978641479e-6, 4.5487080288425235e-9],
        ]),
    }
    values = datasets[float(ionic_strength)]
    return values[:, 0] * PFOA.molecular_weight, values[:, 1], "Experimental data"


def pfos_state(
    concentration: np.ndarray, ionic_strength: float, k1: float, k2: float, gamma_m: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate the state variables from the PFOS mass-action model used for the exported curves."""
    gamma = np.empty_like(concentration)
    binding_ratio = np.empty_like(concentration)
    psi = np.empty_like(concentration)
    theta = np.empty_like(concentration)
    phi = np.empty_like(concentration)

    for index, c_raw in enumerate(concentration):
        c_total = ionic_strength + c_raw
        c_surfactant = c_raw
        c_electrolyte = c_total

        def charge_balance(potential: float) -> float:
            exp_term = np.exp(F * potential / (R * T))
            adsorption = k1 * c_surfactant * exp_term + k2 * c_electrolyte * c_surfactant
            coverage = adsorption / (1.0 + adsorption)
            gamma_value = gamma_m * coverage
            counterion_excess = gamma_value * k2 * c_electrolyte / (k1 * exp_term + k2 * c_electrolyte)
            return F * (-gamma_value + counterion_excess) - np.sqrt(8.0 * EPSILON * R * T * c_electrolyte) * np.sinh(
                F * potential / (2.0 * R * T)
            )

        potential = brentq(charge_balance, -0.3, 0.0, xtol=1e-13, rtol=1e-13)
        exp_term = np.exp(F * potential / (R * T))
        adsorption = k1 * c_surfactant * exp_term + k2 * c_electrolyte * c_surfactant
        coverage = adsorption / (1.0 + adsorption)
        gamma[index] = gamma_m * coverage
        binding_ratio[index] = k2 * c_electrolyte / (k1 * exp_term + k2 * c_electrolyte)
        psi[index] = potential
        theta[index] = coverage
        phi[index] = -F * potential / (R * T)

    return gamma, binding_ratio, psi, theta, phi


def pfos_curves() -> tuple[
    Compound,
    list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    list[tuple[np.ndarray, np.ndarray, str] | None],
]:
    """Evaluate updated PFOS mass action curves and load recorded Langmuir fits and measurements."""
    curve_path = ROOT / "PFAS_plotResults_package" / "output_matlab_PFOS_Langmuir_vs_Kralchevsky" / "PFOS_Gamma_model_comparison_curves.csv"
    fit_path = ROOT / "PFAS_plotResults_package" / "output_matlab_PFOS_Langmuir_vs_Kralchevsky" / "PFOS_Langmuir_fit_summary.csv"
    point_path = ROOT / "PFAS_plotResults_package" / "output_matlab_PFOS_surface_excess_models" / "PFOS_surface_excess_experimental_selection.csv"
    curve_data = pd.read_csv(curve_path)
    fit_data = pd.read_csv(fit_path)
    point_data = pd.read_csv(point_path)

    k1 = 1553.0
    k2 = 3.8825
    gamma_m = 3.0695e-6
    ionic_strengths = (1.0, 10.0, 50.0)
    fit_rows = [fit_data.loc[fit_data["NaCl_mM"].eq(value)].iloc[0] for value in ionic_strengths]
    k_langmuir_values = tuple(float(row["K_Langmuir_m3_mol"]) for row in fit_rows)
    gamma_m_langmuir_values = tuple(float(row["Gamma_m_Langmuir_mol_m2"]) for row in fit_rows)
    compound = Compound(
        "PFOS", 499.0, ionic_strengths,
        (k1,) * 3, (k2,) * 3, (gamma_m,) * 3,
        k_langmuir_values, gamma_m_langmuir_values,
        (1e-8, 1e3), (1e-5, 1e6),
    )

    curves = []
    points = []
    for ionic_strength in compound.ionic_strengths:
        rows = curve_data["NaCl_mM"].eq(ionic_strength)
        c_reference = curve_data.loc[rows, "Concentration_mol_m3"].to_numpy(dtype=float)
        gamma_langmuir_reference = curve_data.loc[rows, "Gamma_Langmuir_umol_m2"].to_numpy(dtype=float) * 1e-6
        fit_row = fit_data[fit_data["NaCl_mM"].eq(ionic_strength)]
        if len(fit_row) != 1:
            raise RuntimeError(f"Missing PFOS Langmuir fit for {ionic_strength:g} mM.")
        k_langmuir = float(fit_row["K_Langmuir_m3_mol"].iloc[0])
        gamma_m_langmuir = float(fit_row["Gamma_m_Langmuir_mol_m2"].iloc[0])
        reproduced_langmuir = gamma_m_langmuir * k_langmuir * c_reference / (1.0 + k_langmuir * c_reference)
        if not np.allclose(reproduced_langmuir, gamma_langmuir_reference, rtol=2e-6, atol=1e-13):
            raise RuntimeError(f"PFOS Langmuir fit does not reproduce the exported {ionic_strength:g} mM curve.")

        # Extend only the deterministic model predictions beyond the supplied curve table;
        # direct measurements remain limited to the recorded 1 and 10 mM points below.
        c = np.logspace(-8, 3, 480)
        gamma_full, binding_ratio, psi, theta, phi = pfos_state(c, ionic_strength, k1, k2, gamma_m)
        gamma_langmuir = gamma_m_langmuir * k_langmuir * c / (1.0 + k_langmuir * c)
        codes = classify_regimes(theta, phi, ionic_strength, k1, k2)
        curves.append((c * compound.molecular_weight, gamma_full, gamma_langmuir, binding_ratio, psi, codes))

        direct = point_data[
            point_data["Include_in_primary_comparison"].eq(1) & point_data["NaCl_mM"].eq(ionic_strength)
        ]
        if direct.empty:
            points.append(None)
        else:
            points.append((
                direct["PFOS_concentration_mol_m3"].to_numpy(dtype=float) * compound.molecular_weight,
                direct["Surface_excess_umol_m2"].to_numpy(dtype=float) * 1e-6,
                "Direct measurement",
            ))
    return compound, curves, points


def calculate_fit_implied_penalties(
    compound: Compound,
) -> tuple[tuple[float, float, float], list[dict[str, float | str]]]:
    """Infer the constant electrostatic penalty from the two fitted dilute-limit slopes."""
    penalties = []
    summary = []
    for index, ionic_strength in enumerate(compound.ionic_strengths):
        binding_term = (compound.k2[index] / compound.k1[index]) * ionic_strength
        full_slope = compound.gamma_m_full[index] * compound.k1[index] * (1.0 + binding_term)
        langmuir_slope = compound.gamma_m_langmuir[index] * compound.k_langmuir[index]
        lambda_fit = float(np.log(full_slope / langmuir_slope))
        penalties.append(lambda_fit)
        summary.append({
            "Compound": compound.name,
            "I_mM": ionic_strength,
            "full_model_dilute_slope": full_slope,
            "langmuir_dilute_slope": langmuir_slope,
            "Lambda_el_fit": lambda_fit,
            "definition": "ln(full-model dilute slope / Langmuir dilute slope)",
        })
    return tuple(penalties), summary


def main() -> None:
    pfoa_data = []
    for index, ionic_strength in enumerate(PFOA.ionic_strengths):
        c = np.logspace(np.log10(PFOA.c_range[0]), np.log10(PFOA.c_range[1]), 220)
        gamma_full, binding_ratio, psi, theta, phi = solve_full_curve(
            c, ionic_strength, PFOA.k1[index], PFOA.k2[index], PFOA.gamma_m_full[index]
        )
        gamma_langmuir = PFOA.gamma_m_langmuir[index] * PFOA.k_langmuir[index] * c / (1.0 + PFOA.k_langmuir[index] * c)
        pfoa_data.append((c * PFOA.molecular_weight, gamma_full, gamma_langmuir, binding_ratio, psi,
                          classify_regimes(theta, phi, ionic_strength, PFOA.k1[index], PFOA.k2[index])))

    pfhxa_data = []
    for index, ionic_strength in enumerate(PFHXA.ionic_strengths):
        c = np.logspace(np.log10(PFHXA.c_range[0]), np.log10(PFHXA.c_range[1]), 220)
        gamma_full, binding_ratio, psi, theta, phi = solve_full_curve(
            c, ionic_strength, PFHXA.k1[index], PFHXA.k2[index], PFHXA.gamma_m_full[index]
        )
        gamma_langmuir = PFHXA.gamma_m_langmuir[index] * PFHXA.k_langmuir[index] * c / (1.0 + PFHXA.k_langmuir[index] * c)
        pfhxa_data.append((c * PFHXA.molecular_weight, gamma_full, gamma_langmuir, binding_ratio, psi,
                           classify_regimes(theta, phi, ionic_strength, PFHXA.k1[index], PFHXA.k2[index])))

    pfos, pfos_data, pfos_points = pfos_curves()

    fit_penalties = {}
    fit_penalty_summary = []
    for compound in (PFOA, PFHXA, pfos):
        values, summary_rows = calculate_fit_implied_penalties(compound)
        fit_penalties[compound.name] = values
        fit_penalty_summary.extend(summary_rows)
    pd.DataFrame(fit_penalty_summary).to_csv(ROOT / "Lambda_el_fit_summary.csv", index=False)

    # Preserve A4 landscape width while adding height for clearer separation between rows.
    fig, axes = plt.subplots(3, 4, figsize=(11.69, 9.25))
    fig.subplots_adjust(left=0.08, right=0.95, bottom=0.20, top=0.91, wspace=0.50, hspace=0.80)
    for axis in axes[:, 3]:
        position = axis.get_position()
        axis.set_position([position.x0 + 0.020, position.y0, position.width, position.height])
    rows = ((PFOA, pfoa_data), (PFHXA, pfhxa_data), (pfos, pfos_data))
    pfoa_points = {
        0: pfoa_experimental_data(10.0),
    }
    for row_index, (compound, records) in enumerate(rows):
        for column_index, record in enumerate(records):
            x, gamma_full, gamma_langmuir, binding_ratio, psi, codes = record
            points = None
            if compound.name == "PFOA" and column_index in pfoa_points:
                points = pfoa_points[column_index]
            elif compound.name == "PFOS":
                points = pfos_points[column_index]
            plot_isotherm_panel(
                axes[row_index, column_index], f"({chr(97 + row_index)}{column_index + 1})",
                compound, compound.ionic_strengths[column_index],
                gamma_full, gamma_langmuir, x, codes, show_legend=False,
                show_y_label=(column_index == 0), scatter_data=points,
            )
        plot_ratio_panel(
            axes[row_index, 3], f"({chr(97 + row_index)}4)", compound, records, fit_penalties[compound.name]
        )

    for axis in axes.flat:
        axis.set_facecolor("white")

    model_legend = fig.legend(
        handles=[
            mpl.lines.Line2D([], [], color="black", lw=1.8, label="Mass action model"),
            mpl.lines.Line2D([], [], color="#59ad71", lw=1.8, ls=":", label="Langmuir isotherm"),
            mpl.lines.Line2D([], [], color=EXPERIMENT_COLOR, marker="o", ls="", markersize=6, label="Experimental data"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.40, 0.070),
        ncol=3,
        fontsize=12,
        frameon=False,
        handlelength=1.8,
        columnspacing=0.7,
        handletextpad=0.3,
    )
    fig.add_artist(model_legend)
    regime_legend = fig.legend(
        handles=[
            mpl.patches.Patch(
                facecolor=REGIME_COLORS[0], edgecolor="none", label="adsorption-affinity-controlled"
            ),
            mpl.patches.Patch(
                facecolor=REGIME_COLORS[1], edgecolor="none", label="electrostatic-inhibition"
            ),
            mpl.patches.Patch(
                facecolor=REGIME_COLORS[2], edgecolor="none", label="steric-exclusion"
            ),
        ],
        loc="lower center",
        bbox_to_anchor=(0.40, 0.015),
        ncol=3,
        fontsize=12,
        frameon=False,
        handlelength=1.0,
        columnspacing=0.9,
        handletextpad=0.3,
    )
    fig.add_artist(regime_legend)

    output_bases = (
        ROOT / "Full_vs_Langmuir_Isotherms_PFOA_PFHxA_PFOS",
        ROOT / "Full_vs_Langmuir_Isotherms",
    )
    publication_pdf = ROOT / "NatureWater_Overleaf" / "figures" / "fig5_full_vs_langmuir.pdf"
    publication_pdf.parent.mkdir(parents=True, exist_ok=True)
    for output_base in output_bases:
        fig.savefig(output_base.with_suffix(".svg"))
        fig.savefig(output_base.with_suffix(".png"), dpi=600)
    fig.savefig(publication_pdf)
    plt.close(fig)
    print(
        "Exported SVG/PNG previews and publication PDF: "
        + str(publication_pdf.relative_to(ROOT))
    )


if __name__ == "__main__":
    main()
