"""Create an eight-panel comparison for selected PFAS and ionic strengths."""

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogFormatterMathtext, LogLocator, MaxNLocator, NullFormatter
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd

import plot_full_vs_langmuir_three_pfas as base


ROOT = Path(__file__).resolve().parent
MASS_ACTION_COLOR = "#557D9A"
LANGMUIR_COLOR = "#C66F61"
RATIO_COLOR = "#7C7097"
PFAS_MARKERS = {
    "PFOA": "o",
    "PFHxA": "s",
    "PFOS": "^",
}
PFAS_COLORS = {
    "PFOA": "#557D9A",
    "PFHxA": "#D0A25D",
    "PFOS": "#6E9D88",
}
PFOA_ALL_IONIC_STRENGTHS = base.Compound(
    "PFOA", 413.0, (10.0, 50.0, 100.0, 200.0, 300.0, 500.0),
    (49.2,) * 6, (0.13,) * 6, (3.36e-6,) * 6,
    (6.227663, 20.634773, 25.406345, 44.640446, 43.536832, 75.203134),
    (2.988939e-6, 2.937332e-6, 3.216727e-6, 3.297605e-6, 3.441808e-6, 3.629248e-6),
    base.PFOA.c_range, base.PFOA.x_limits,
)
PFHXA_ALL_IONIC_STRENGTHS = base.Compound(
    "PFHxA", 313.0, (100.0, 200.0, 300.0, 500.0),
    (1.69,) * 4, (0.0046,) * 4, (2.56e-6,) * 4,
    (0.678870, 1.259270, 1.671099, 2.440557),
    (3.002462e-6, 2.869480e-6, 2.869755e-6, 2.955345e-6),
    base.PFHXA.c_range, base.PFHXA.x_limits,
)
ION_COLORS = {
    # Broad aurora spectrum: indigo, blue, teal, green, and yellow-green.
    1.0: "#482475",
    10.0: "#355F8D",
    50.0: "#2A788E",
    100.0: "#53B66D",
    200.0: "#9FD63B",
}
PANEL_LABELS = (("a", "b", "c"), ("d", "e", "f"))
REGIME_LABELS = {
    0: "surface-\nactivity-\ncontrolled",
    1: "electrostatic-\ninhibition",
    2: "steric-\nexclusion",
}
REGIME_TEXT_COLORS = {0: "#4f809f", 1: "#666666", 2: "#508f6c"}
REGIME_TEXT_Y_POSITIONS = {0: 0.90, 1: 0.23, 2: 0.24}
PANEL_X_LIMITS = (1e-5, 1e6)


def format_alternate_decade_tick(value: float, position: int) -> str:
    """Show every second decade label while retaining all decade tick marks."""
    if value <= 0.0:
        return ""
    exponent = int(np.round(np.log10(value)))
    if not np.isclose(value, 10.0 ** exponent) or exponent % 2 == 0:
        return ""
    return base.LOG_TICK_FORMATTER(value, position)


ALTERNATE_LOG_TICK_FORMATTER = FuncFormatter(format_alternate_decade_tick)


def subdued_point_edge(color: str, text_fraction: float = 0.35) -> tuple[float, ...]:
    """Create a low-chroma, same-family marker edge from a palette color."""
    color_rgb = np.asarray(mpl.colors.to_rgb(color))
    text_rgb = np.asarray(mpl.colors.to_rgb("#282D32"))
    return tuple((1.0 - text_fraction) * color_rgb + text_fraction * text_rgb)


def calculate_records(compound: base.Compound, indices: tuple[int, ...]) -> list[tuple[np.ndarray, ...]]:
    records = []
    for index in indices:
        ionic_strength = compound.ionic_strengths[index]
        concentration = np.logspace(
            np.log10(PANEL_X_LIMITS[0] / compound.molecular_weight),
            np.log10(PANEL_X_LIMITS[1] / compound.molecular_weight),
            320,
        )
        gamma_full, binding_ratio, psi, theta, phi = base.solve_full_curve(
            concentration,
            ionic_strength,
            compound.k1[index],
            compound.k2[index],
            compound.gamma_m_full[index],
        )
        gamma_langmuir = (
            compound.gamma_m_langmuir[index]
            * compound.k_langmuir[index]
            * concentration
            / (1.0 + compound.k_langmuir[index] * concentration)
        )
        regimes = base.classify_regimes(
            theta, phi, ionic_strength, compound.k1[index], compound.k2[index]
        )
        records.append(
            (
                concentration * compound.molecular_weight,
                gamma_full,
                gamma_langmuir,
                binding_ratio,
                psi,
                regimes,
            )
        )
    return records


def annotate_regime_names(ax: plt.Axes, x: np.ndarray, regimes: np.ndarray) -> None:
    """Place compact regime names above their visible background bands."""
    x_min, x_max = ax.get_xlim()
    for code, label in REGIME_LABELS.items():
        indices = np.flatnonzero((regimes == code) & (x >= x_min) & (x <= x_max))
        if indices.size == 0:
            continue
        left = max(float(x[indices[0]]), x_min)
        right = min(float(x[indices[-1]]), x_max)
        center = np.sqrt(left * right)
        center_fraction = (
            np.log10(center) - np.log10(x_min)
        ) / (np.log10(x_max) - np.log10(x_min))
        y_position = REGIME_TEXT_Y_POSITIONS[code]
        ax.text(
            center_fraction,
            y_position,
            label,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=11,
            color=REGIME_TEXT_COLORS[code],
            linespacing=0.95,
            zorder=8,
            bbox={
                "facecolor": base.REGIME_COLORS[code],
                "edgecolor": "none",
                "alpha": 0.92,
                "pad": 0.25,
            },
        )


def plot_isotherm_ratio_panel(
    ax: plt.Axes,
    letter: str,
    compound: base.Compound,
    ionic_strength: float,
    record: tuple[np.ndarray, ...],
    scatter_data: tuple[np.ndarray, np.ndarray, str] | None,
    annotate_regimes: bool,
    show_left_label: bool,
    show_right_label: bool,
    panel_label_x_position: float,
) -> plt.Axes:
    x, gamma_full, gamma_langmuir, _, _, regimes = record
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*PANEL_X_LIMITS)
    ax.set_ylim(1e-5, 1e1)
    base.regime_background(ax, x, regimes)
    ax.plot(
        x,
        gamma_full * 1e6,
        color=MASS_ACTION_COLOR,
        lw=1.8,
        zorder=3,
    )
    ax.plot(
        x,
        gamma_langmuir * 1e6,
        color=LANGMUIR_COLOR,
        lw=1.7,
        ls=(0, (3.4, 2.2)),
        zorder=3,
    )
    if scatter_data is not None:
        ax.scatter(
            scatter_data[0],
            scatter_data[1] * 1e6,
            s=56,
            facecolors=mpl.colors.to_rgba(base.EXPERIMENT_COLOR, alpha=0.5),
            edgecolors=base.EXPERIMENT_COLOR,
            linewidths=1.0,
            zorder=6,
        )
    y_ticks = 10.0 ** np.arange(-5, 2)
    x_decades = np.arange(
        np.ceil(np.log10(PANEL_X_LIMITS[0])),
        np.floor(np.log10(PANEL_X_LIMITS[1])) + 1,
    )
    ax.set_yticks(y_ticks)
    ax.set_xticks(10.0 ** x_decades)
    ax.xaxis.set_major_formatter(ALTERNATE_LOG_TICK_FORMATTER)
    ax.yaxis.set_major_formatter(ALTERNATE_LOG_TICK_FORMATTER)
    ax.tick_params(which="major", direction="out", width=0.8, length=3.2, labelsize=11)
    ax.tick_params(which="minor", length=0)
    ax.set_title(f"{compound.name}, I = {ionic_strength:g} mM", fontsize=12, pad=6)
    if show_left_label:
        ax.set_ylabel("Surface excess\n($\\mu$mol m$^{-2}$)", fontsize=12)
    else:
        ax.tick_params(axis="y", labelleft=False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
    base.add_panel_label(
        ax, letter, x_position=panel_label_x_position, y_position=1.04
    )

    ratio_ax = ax.twinx()
    ratio = gamma_langmuir / gamma_full
    ratio_ax.plot(x, ratio, color=RATIO_COLOR, lw=1.45, zorder=5)
    ratio_ax.set_ylim(0.0, 1.7)
    ratio_ax.set_yticks([0.0, 0.5, 1.0, 1.5])
    ratio_ax.tick_params(
        which="major",
        direction="out",
        width=0.8,
        length=3.2,
        labelsize=11,
        pad=2,
        colors=RATIO_COLOR,
    )
    ratio_ax.tick_params(which="minor", length=0)
    ratio_ax.patch.set_visible(False)
    for spine in ratio_ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
    ratio_ax.spines["right"].set_color(RATIO_COLOR)
    if show_right_label:
        ratio_ax.set_ylabel(
            r"$\Gamma_{\mathrm{SL}}/\Gamma_{\mathrm{MA}}$",
            fontsize=12,
            rotation=90,
            labelpad=16,
            color=RATIO_COLOR,
        )

    if annotate_regimes:
        annotate_regime_names(ax, x, regimes)
    return ratio_ax


def plot_integrated_panel(
    ax: plt.Axes,
    selected_rows: list[tuple[base.Compound, list[tuple[np.ndarray, ...]]]],
    fit_penalties: dict[str, tuple[float, float, float]],
    panel_label: str,
    surface_activity_only: bool = False,
    show_legends: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    plotted_x = []
    plotted_y = []
    source_rows = []

    plot_order = ("PFOS", "PFOA", "PFHxA")
    rows_by_name = {compound.name: (compound, records) for compound, records in selected_rows}
    point_zorders = {"PFOS": 3, "PFOA": 4, "PFHxA": 5}

    for compound_name in plot_order:
        compound, records = rows_by_name[compound_name]
        for index, record in enumerate(records):
            concentration_mg_l, gamma_full, gamma_langmuir, _, psi, regimes = record
            ionic_strength = float(compound.ionic_strengths[index])
            phi = -base.F * psi / (base.R * base.T)
            lambda_ele = base.calculate_electrostatic_penalty(
                phi, ionic_strength, compound.k1[index], compound.k2[index]
            )
            delta_lambda_ele = fit_penalties[compound.name][index] - lambda_ele
            ratio = gamma_langmuir / gamma_full
            theta_full = np.clip(
                gamma_full / compound.gamma_m_full[index], 0.0, 1.0 - 1e-12
            )
            steric_penalty = -np.log1p(-theta_full)
            theta_langmuir = np.clip(
                gamma_langmuir / compound.gamma_m_langmuir[index], 0.0, 1.0 - 1e-12
            )
            steric_penalty_langmuir = -np.log1p(-theta_langmuir)
            eligible = (
                (steric_penalty <= base.PENALTY_THRESHOLD)
                & (steric_penalty_langmuir <= base.PENALTY_THRESHOLD)
            )
            if surface_activity_only:
                eligible &= regimes == 0
            selected = np.flatnonzero(eligible)
            if selected.size == 0:
                raise RuntimeError(
                    f"No points satisfy both low-coverage conditions for {compound.name}, {ionic_strength:g} mM."
                )
            sampled = np.unique(np.r_[selected[::15], selected[-1]])
            displayed = np.zeros(eligible.size, dtype=bool)
            displayed[sampled] = True

            color = ION_COLORS[ionic_strength]
            marker = PFAS_MARKERS[compound.name]
            ax.scatter(
                delta_lambda_ele[sampled],
                ratio[sampled],
                s=92,
                marker=marker,
                facecolors=mpl.colors.to_rgba(color, alpha=0.58),
                edgecolors=mpl.colors.to_rgba(subdued_point_edge(color), alpha=0.78),
                linewidths=0.8,
                zorder=point_zorders[compound.name],
            )
            plotted_x.append(delta_lambda_ele[selected])
            plotted_y.append(ratio[selected])

            for point_index in selected:
                source_rows.append(
                    {
                        "PFAS": compound.name,
                        "I_mM": ionic_strength,
                        "Panel": panel_label.strip("()"),
                        "Delta_Lambda_ele": float(delta_lambda_ele[point_index]),
                        "GammaSL_over_GammaMA": float(ratio[point_index]),
                        "Concentration_mg_L": float(concentration_mg_l[point_index]),
                        "theta_MA": float(theta_full[point_index]),
                        "theta_SL": float(theta_langmuir[point_index]),
                        "Lambda_st_MA": float(steric_penalty[point_index]),
                        "Lambda_st_SL": float(steric_penalty_langmuir[point_index]),
                        "Reference_exp_neg_delta": float(np.exp(-delta_lambda_ele[point_index])),
                        "Displayed": bool(displayed[point_index]),
                        "marker": marker,
                        "color": color,
                    }
                )

    x_min = min(values.min() for values in plotted_x)
    x_max = max(values.max() for values in plotted_x)
    locator = MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10])
    x_ticks = locator.tick_values(x_min, x_max)
    x_lower = float(x_ticks[0])
    x_upper = float(x_ticks[-1])
    all_x = np.concatenate(plotted_x)
    all_y = np.concatenate(plotted_y)
    y_upper = max(1.05, float(np.ceil(np.max(all_y) * 4.0) / 4.0))

    theoretical_y = np.exp(-all_x)
    residual_sum_squares = np.sum((all_y - theoretical_y) ** 2)
    total_sum_squares = np.sum((all_y - np.mean(all_y)) ** 2)
    r_squared = 1.0 - residual_sum_squares / total_sum_squares
    rmse = np.sqrt(np.mean((all_y - theoretical_y) ** 2))

    ax.set_xlim(x_lower, x_upper)
    ax.set_ylim(0.0, y_upper)
    ax.set_xticks(x_ticks)
    ax.set_yticks(np.arange(0.0, y_upper + 1e-9, 0.25))
    master_x = np.linspace(x_lower, x_upper, 300)
    ax.plot(
        master_x,
        np.exp(-master_x),
        color="#707070",
        lw=1.25,
        ls=(0, (3, 2)),
        zorder=7,
    )
    annotation_x = 0.62
    annotation_y = 0.20
    ax.text(
        annotation_x,
        annotation_y + 0.18,
        r"$\Gamma_{\mathrm{SL}}/\Gamma_{\mathrm{MA}}"
        r" = \exp(-\Delta\Lambda_{\mathrm{ele}})$",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        color="#333333",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.90, "pad": 1.0},
    )
    ax.text(
        annotation_x,
        annotation_y + 0.09,
        f"$R^2$ = {r_squared:.5f}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        color="#333333",
    )
    ax.text(
        annotation_x,
        annotation_y,
        f"RMSE = {rmse:.4f}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        color="#333333",
    )
    ax.set_xlabel(r"$\Delta\Lambda_{\mathrm{ele}}$", fontsize=12)
    ax.set_ylabel(r"$\Gamma_{\mathrm{SL}}/\Gamma_{\mathrm{MA}}$", fontsize=12)
    ax.tick_params(which="major", direction="out", width=0.8, length=3.2, labelsize=11)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
    base.add_panel_label(ax, panel_label, x_position=-0.12, y_position=1.025)

    if surface_activity_only:
        ax.set_title("Surface-activity-controlled regime", fontsize=11, pad=5)

    pfas_handles = [
        mpl.lines.Line2D(
            [],
            [],
            color="#555555",
            marker=PFAS_MARKERS[name],
            markerfacecolor=mpl.colors.to_rgba("#555555", alpha=0.5),
            markeredgewidth=1.2,
            ls="",
            markersize=6.0,
            label=name,
        )
        for name in ("PFOA", "PFHxA", "PFOS")
    ]
    if show_legends:
        pfas_legend = ax.legend(
            handles=pfas_handles,
            loc="upper right",
            fontsize=11,
            frameon=False,
            handlelength=1.2,
            handletextpad=0.5,
            labelspacing=0.3,
        )
        ax.add_artist(pfas_legend)

    ion_handles = [
        mpl.lines.Line2D(
            [],
            [],
            color=ION_COLORS[ionic_strength],
            lw=2.2,
            label=f"{ionic_strength:g}",
        )
        for ionic_strength in sorted(ION_COLORS)
    ]
    if show_legends:
        ax.legend(
            handles=ion_handles,
            title="I (mM)",
            loc="lower left",
            ncol=2,
            fontsize=11,
            title_fontsize=11,
            frameon=False,
            handlelength=1.2,
            handletextpad=0.5,
            labelspacing=0.3,
        )
    fit_summary = pd.DataFrame(
        [
            {
                "n_points": len(all_x),
                "Panel": panel_label.strip("()"),
                "selection": (
                    "both steric penalties <= 0.1, surface-activity-controlled"
                    if surface_activity_only
                    else "both steric penalties <= 0.1, all eligible points"
                ),
                "model": "fixed theoretical exponential",
                "amplitude": 1.0,
                "decay_constant": 1.0,
                "R_squared": r_squared,
                "RMSE": rmse,
                "theoretical_amplitude": 1.0,
                "theoretical_decay_constant": 1.0,
            }
        ]
    )
    return pd.DataFrame(source_rows), fit_summary


def add_shared_legend(fig: plt.Figure, center_x: float, y_position: float) -> None:
    model_legend = fig.legend(
        handles=[
            mpl.lines.Line2D(
                [], [], color=MASS_ACTION_COLOR, lw=1.8, label="Mass action model"
            ),
            mpl.lines.Line2D(
                [], [], color=LANGMUIR_COLOR, lw=1.7, ls=(0, (3.4, 2.2)), label="Langmuir isotherm"
            ),
            mpl.lines.Line2D(
                [],
                [],
                color=RATIO_COLOR,
                lw=1.45,
                label=r"$\Gamma_{\mathrm{SL}}/\Gamma_{\mathrm{MA}}$",
            ),
            mpl.lines.Line2D(
                [],
                [],
                markerfacecolor=mpl.colors.to_rgba(base.EXPERIMENT_COLOR, alpha=0.5),
                markeredgecolor=base.EXPERIMENT_COLOR,
                markeredgewidth=1.0,
                marker="o",
                ls="",
                markersize=7,
                label="Experimental data",
            ),
        ],
        loc="lower center",
        bbox_to_anchor=(center_x, y_position),
        ncol=4,
        fontsize=11,
        frameon=False,
        handlelength=1.5,
        columnspacing=0.9,
        handletextpad=0.3,
    )
    fig.add_artist(model_legend)


def calculate_minimum_gamma_ratio(
    compound: base.Compound, index: int
) -> tuple[float, float, str]:
    """Return the minimum Langmuir-to-mass-action ratio for one fitted series."""
    concentration = np.logspace(
        np.log10(compound.c_range[0]), np.log10(compound.c_range[1]), 640
    )
    ionic_strength = compound.ionic_strengths[index]
    if compound.name == "PFOS":
        gamma_full, *_ = base.pfos_state(
            concentration,
            ionic_strength,
            compound.k1[index],
            compound.k2[index],
            compound.gamma_m_full[index],
        )
    else:
        gamma_full, *_ = base.solve_full_curve(
            concentration,
            ionic_strength,
            compound.k1[index],
            compound.k2[index],
            compound.gamma_m_full[index],
        )
    gamma_langmuir = (
        compound.gamma_m_langmuir[index]
        * compound.k_langmuir[index]
        * concentration
        / (1.0 + compound.k_langmuir[index] * concentration)
    )
    ratio = gamma_langmuir / gamma_full
    sampled_index = int(np.argmin(ratio))
    sampled_minimum = float(ratio[sampled_index])
    binding_term = (compound.k2[index] / compound.k1[index]) * ionic_strength
    dilute_limit_ratio = float(
        compound.gamma_m_langmuir[index]
        * compound.k_langmuir[index]
        / (compound.gamma_m_full[index] * compound.k1[index] * (1.0 + binding_term))
    )
    if dilute_limit_ratio < sampled_minimum:
        return dilute_limit_ratio, 0.0, "dilute limit"
    return (
        sampled_minimum,
        float(concentration[sampled_index] * compound.molecular_weight),
        "model curve",
    )


def plot_minimum_gamma_ratio_panel(
    ax: plt.Axes,
    compounds: tuple[base.Compound, ...],
) -> pd.DataFrame:
    """Show the minimum Langmuir-to-mass-action ratio for each salt series."""
    source_rows = []
    for compound in compounds:
        ionic_strengths = np.asarray(compound.ionic_strengths, dtype=float)
        minima = [
            calculate_minimum_gamma_ratio(compound, index)
            for index in range(len(ionic_strengths))
        ]
        minimum_ratios = [minimum[0] for minimum in minima]
        color = PFAS_COLORS[compound.name]
        marker = PFAS_MARKERS[compound.name]
        ax.scatter(
            ionic_strengths,
            minimum_ratios,
            s=92,
            marker=marker,
            facecolors=mpl.colors.to_rgba(color, alpha=0.65),
            edgecolors=mpl.colors.to_rgba(subdued_point_edge(color), alpha=0.80),
            linewidths=0.8,
            label=compound.name,
            zorder=4,
        )
        for ionic_strength, minimum in zip(ionic_strengths, minima):
            minimum_ratio, concentration_at_minimum, minimum_source = minimum
            source_rows.append(
                {
                    "PFAS": compound.name,
                    "I_mM": ionic_strength,
                    "GammaSL_over_GammaMA_min": minimum_ratio,
                    "Concentration_at_min_mg_L": concentration_at_minimum,
                    "minimum_source": minimum_source,
                }
            )

    source_data = pd.DataFrame(source_rows)
    log_ionic_strength = np.log10(source_data["I_mM"].to_numpy(dtype=float))
    log_minimum_ratio = np.log10(
        source_data["GammaSL_over_GammaMA_min"].to_numpy(dtype=float)
    )
    slope, intercept = np.polyfit(log_ionic_strength, log_minimum_ratio, 1)
    predicted_log_ratio = slope * log_ionic_strength + intercept
    residual_sum_squares = np.sum((log_minimum_ratio - predicted_log_ratio) ** 2)
    total_sum_squares = np.sum((log_minimum_ratio - np.mean(log_minimum_ratio)) ** 2)
    r_squared = 1.0 - residual_sum_squares / total_sum_squares
    rmse_log10 = np.sqrt(np.mean((log_minimum_ratio - predicted_log_ratio) ** 2))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.7, 700.0)
    ax.set_ylim(5e-3, 1.0)
    ax.set_xticks([1.0, 10.0, 100.0, 500.0])
    ax.set_xticklabels(["1", "10", "100", "500"])
    ax.yaxis.set_major_locator(LogLocator(base=10.0))
    ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("Ionic strength (mM)", fontsize=12)
    ax.set_ylabel(
        r"$\left(\Gamma_{\mathrm{SL}}/\Gamma_{\mathrm{MA}}\right)_{\min}$",
        fontsize=12,
    )
    ax.tick_params(which="major", direction="out", width=0.8, length=3.2, labelsize=11)
    ax.tick_params(which="minor", direction="out", width=0.6, length=2.0)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.55, alpha=0.65)
    fit_ionic_strength = np.logspace(0.0, np.log10(500.0), 240)
    ax.plot(
        fit_ionic_strength,
        10.0 ** (intercept + slope * np.log10(fit_ionic_strength)),
        color="#707070",
        lw=1.25,
        ls=(0, (3, 2)),
        label="Pooled empirical fit",
        zorder=2,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
    ax.legend(
        loc="upper left",
        ncol=2,
        fontsize=11,
        frameon=False,
        handlelength=1.4,
        handletextpad=0.5,
        labelspacing=0.3,
    )
    ax.text(
        0.58,
        0.10,
        rf"$y_{{\min}} = 10^{{{intercept:.2f}}} I^{{{slope:.2f}}}$"
        + f"\n$R^2$ = {r_squared:.2f}"
        + "\n"
        + rf"$\mathrm{{RMSE}}_{{\log_{{10}}}}$ = {rmse_log10:.2f}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        multialignment="left",
        fontsize=11,
        color="#333333",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.90, "pad": 1.0},
    )
    base.add_panel_label(ax, "h", x_position=-0.12, y_position=1.025)
    return source_data


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path,
        help="Export the complete figure and its source data to this directory.",
    )
    args = parser.parse_args()
    pfoa_records = calculate_records(base.PFOA, (0, 1))
    pfhxa_records = calculate_records(base.PFHXA, (0, 1))
    pfos, pfos_records_all, pfos_points_all = base.pfos_curves()
    pfos_records = pfos_records_all[:2]
    pfos_points = pfos_points_all[:2]

    selected_rows = [
        (base.PFOA, pfoa_records),
        (base.PFHXA, pfhxa_records),
        (pfos, pfos_records),
    ]
    fit_penalties = {
        compound.name: base.calculate_fit_implied_penalties(compound)[0]
        for compound, _ in selected_rows
    }

    fig = plt.figure(figsize=(11.69, 11.3))
    top_grid_bottom = 0.46
    top_grid = fig.add_gridspec(
        2,
        3,
        left=0.105,
        right=0.910,
        bottom=top_grid_bottom,
        top=0.94,
        wspace=0.20,
        hspace=0.32,
    )
    top_axes = np.empty((2, 3), dtype=object)
    ratio_axes = np.empty((2, 3), dtype=object)
    for row_index in range(2):
        for column_index in range(3):
            share_x_axis = top_axes[0, column_index] if row_index else None
            share_y_axis = top_axes[0, 0] if (row_index or column_index) else None
            top_axes[row_index, column_index] = fig.add_subplot(
                top_grid[row_index, column_index],
                sharex=share_x_axis,
                sharey=share_y_axis,
            )

    pfoa_points = {0: base.pfoa_experimental_data(10.0)}
    point_maps = (pfoa_points, {}, {0: pfos_points[0], 1: pfos_points[1]})

    panel_index = 0
    for compound_index, (compound, records) in enumerate(selected_rows):
        points_by_column = point_maps[compound_index]
        for record_index, record in enumerate(records):
            row_index, column_index = divmod(panel_index, 3)
            ratio_axes[row_index, column_index] = plot_isotherm_ratio_panel(
                top_axes[row_index, column_index],
                PANEL_LABELS[row_index][column_index],
                compound,
                compound.ionic_strengths[record_index],
                record,
                points_by_column.get(record_index),
                annotate_regimes=(panel_index == 1),
                show_left_label=(column_index == 0),
                show_right_label=(column_index == 2),
                panel_label_x_position=(-0.14 if column_index == 0 else -0.08),
            )
            if row_index < 1:
                top_axes[row_index, column_index].tick_params(
                    axis="x", which="both", labelbottom=False
                )
            panel_index += 1

    top_grid_left = top_axes[0, 0].get_position().x0
    top_grid_right = top_axes[0, 2].get_position().x1
    for column_index in range(3):
        column_axis = top_axes[1, column_index]
        column_position = column_axis.get_position()
        fig.text(
            column_position.x0 + column_position.width / 2,
            top_grid_bottom - 0.037,
            r"Concentration (mg L$^{-1}$)",
            ha="center",
            va="center",
            fontsize=12,
        )
    add_shared_legend(
        fig,
        (top_grid_left + top_grid_right) / 2,
        top_grid_bottom - 0.105,
    )

    master_grid = fig.add_gridspec(
        1,
        2,
        left=top_grid_left,
        right=top_grid_right,
        bottom=0.055,
        top=0.335,
        wspace=0.25,
    )
    integrated_axis = fig.add_subplot(master_grid[0, 0])
    surface_activity_axis = fig.add_subplot(master_grid[0, 1])
    source_data_g, fit_summary_g = plot_integrated_panel(
        integrated_axis,
        selected_rows,
        fit_penalties,
        panel_label="g",
        surface_activity_only=False,
        show_legends=True,
    )
    minimum_ratio_data = plot_minimum_gamma_ratio_panel(
        surface_activity_axis,
        (PFOA_ALL_IONIC_STRENGTHS, PFHXA_ALL_IONIC_STRENGTHS, pfos),
    )
    source_data = source_data_g
    fit_summary = fit_summary_g

    if args.output_dir is None:
        output_base = ROOT / "Selected_PFAS_Full_vs_Langmuir_Master"
        publication_pdf = ROOT / "NatureWater_Overleaf" / "figures" / "fig5_full_vs_langmuir.pdf"
        source_dir = ROOT
    else:
        output_base = args.output_dir.resolve() / "fig5_full_vs_langmuir"
        publication_pdf = output_base.with_suffix(".pdf")
        source_dir = output_base.parent / "source_data" / "fig4"
    publication_pdf.parent.mkdir(parents=True, exist_ok=True)
    source_dir.mkdir(parents=True, exist_ok=True)
    source_data.to_csv(source_dir / "selected_pfas_master_curve_data.csv", index=False)
    fit_summary.to_csv(source_dir / "selected_pfas_master_curve_linearity.csv", index=False)
    minimum_ratio_data.to_csv(source_dir / "selected_pfas_minimum_ratio.csv", index=False)
    curve_tables = []
    observation_tables = []
    parameter_rows = []
    panel_index = 0
    for compound_index, (compound, records) in enumerate(selected_rows):
        for index, record in enumerate(records):
            concentration, gamma_ma, gamma_sl, beta, psi, regimes = record
            panel = chr(ord("a") + panel_index)
            ionic_strength = compound.ionic_strengths[index]
            curve_tables.append(pd.DataFrame({
                "Panel": panel, "PFAS": compound.name, "I_mM": ionic_strength,
                "Concentration_mg_L": concentration,
                "Gamma_MA_umol_m2": gamma_ma * 1e6,
                "Gamma_SL_umol_m2": gamma_sl * 1e6,
                "GammaSL_over_GammaMA": gamma_sl / gamma_ma,
                "Counterion_binding_fraction": beta,
                "Potential_V": psi, "Regime_code": regimes,
            }))
            observed = point_maps[compound_index].get(index)
            if observed is not None:
                observation_tables.append(pd.DataFrame({
                    "Panel": panel, "PFAS": compound.name, "I_mM": ionic_strength,
                    "Concentration_mg_L": observed[0],
                    "Gamma_observed_umol_m2": observed[1] * 1e6,
                    "Reference_key": "buckley2023using" if compound.name == "PFOA" else "schaefer2019uptake",
                }))
            panel_index += 1
    for compound in (PFOA_ALL_IONIC_STRENGTHS, PFHXA_ALL_IONIC_STRENGTHS, pfos):
        for index, ionic_strength in enumerate(compound.ionic_strengths):
            parameter_rows.append({
                "PFAS": compound.name, "I_mM": ionic_strength,
                "K1_m3_mol": compound.k1[index], "K2_m6_mol2": compound.k2[index],
                "Gamma_m_MA_mol_m2": compound.gamma_m_full[index],
                "K_SL_m3_mol": compound.k_langmuir[index],
                "Gamma_m_SL_mol_m2": compound.gamma_m_langmuir[index],
            })
    pd.concat(curve_tables, ignore_index=True).to_csv(source_dir / "fig4_a_f_model_curves.csv", index=False)
    pd.concat(observation_tables, ignore_index=True).to_csv(source_dir / "fig4_a_f_measurements.csv", index=False)
    pd.DataFrame(parameter_rows).to_csv(source_dir / "fig4_parameters.csv", index=False)
    figure_width, figure_height = fig.get_size_inches()
    vertical_crop = Bbox.from_extents(
        0.0,
        0.0,
        figure_width,
        0.97 * figure_height,
    )
    fig.savefig(output_base.with_suffix(".svg"), bbox_inches=vertical_crop)
    fig.savefig(publication_pdf, bbox_inches=vertical_crop)
    fig.savefig(output_base.with_suffix(".png"), dpi=600, bbox_inches=vertical_crop)
    plt.close(fig)
    # Output may be outside the script directory (e.g. --output-dir elsewhere),
    # so report the path as given instead of relative to ROOT.
    print(f"Exported new figure: {publication_pdf}")


if __name__ == "__main__":
    main()
