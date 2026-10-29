"""Plot recorded surface-tension calibrations without refitting parameters."""

from pathlib import Path
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatterMathtext


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    directory = args.project / "SI/source_data/szyszkowski_langmuir_calibration"
    parameters = pd.read_csv(directory / "parameters.csv")
    observations = pd.read_csv(directory / "observations.csv")
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 8, "axes.labelsize": 8,
        "axes.titlesize": 9, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "axes.linewidth": 0.7, "lines.linewidth": 1.15,
        "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False,
        "mathtext.fontset": "dejavusans",
    })
    palette = json.loads((directory / "palette.json").read_text(encoding="utf-8"))
    salts = palette["NaCl_mM"]
    colors = {salt: palette["colors"][index] for salt, index in
              zip(salts, palette["palette_indices"])}
    markers = dict(zip(salts, ["o", "s", "^", "D", "v", "P", "X"]))
    limits = {"PFOA": (1e-2, 1e4), "PFHxA": (1e-1, 1e5), "PFOS": (1e-2, 1e3)}
    ticks = {"PFOA": [1e-2, 1, 1e2, 1e4], "PFHxA": [1e-1, 1e1, 1e3, 1e5],
             "PFOS": [1e-2, 1, 1e2, 1e3]}
    fig, axes = plt.subplots(1, 3, figsize=(183/25.4, 90/25.4), sharey=True)
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.29, top=0.88, wspace=0.18)
    curves = []
    for index, (compound, ax) in enumerate(zip(("PFOA", "PFHxA", "PFOS"), axes)):
        subset = parameters[parameters.PFAS.eq(compound)].sort_values("NaCl_mM")
        for _, fit in subset.iterrows():
            points = observations[observations.PFAS.eq(compound) &
                                  observations.NaCl_mM.eq(fit.NaCl_mM)]
            assert len(points) == fit.n
            if compound == "PFOS":
                # Show each fitted isotherm across the same displayed concentration range.
                concentration = np.geomspace(
                    limits[compound][0]/fit.MW_g_mol,
                    limits[compound][1]/fit.MW_g_mol, 400)
            else:
                concentration = np.geomspace(points.c_mol_m3.min(), points.c_mol_m3.max(), 250)
            predicted = fit.sigma0_mN_m - 1000*8.314*298*fit.Gamma_m_SL_mol_m2 * np.log1p(
                fit.K_SL_m3_mol * concentration)
            color = colors[int(fit.NaCl_mM)]
            marker = markers[int(fit.NaCl_mM)]
            ax.plot(concentration*fit.MW_g_mol, predicted, color=color)
            ax.scatter(points.c_mg_L, points.sigma_mN_m, s=15, marker=marker,
                       facecolors="none", edgecolors=color, linewidths=0.85, zorder=3)
            for c, sigma in zip(concentration, predicted):
                curves.append({"PFAS": compound, "NaCl_mM": fit.NaCl_mM,
                               "c_mol_m3": c, "c_mg_L": c*fit.MW_g_mol,
                               "sigma_fitted_mN_m": sigma})
        ax.set(xscale="log", xlim=limits[compound], ylim=(20, 80),
               xlabel=f"{compound} concentration (mg L$^{{-1}}$)", title=compound)
        ax.set_xticks(ticks[compound])
        ax.xaxis.set_major_formatter(LogFormatterMathtext())
        ax.set_yticks([20, 30, 40, 50, 60, 70, 80])
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(which="major", direction="out", length=3, width=0.7)
        ax.tick_params(which="minor", length=1.8, width=0.5)
        ax.text(-0.08, 1.04, "abc"[index], transform=ax.transAxes, fontweight="bold", fontsize=10)
    axes[0].set_ylabel("Surface tension (mN m$^{-1}$)")
    handles = [Line2D([], [], color=colors[salt], marker=markers[salt],
                      markerfacecolor="none", markersize=3.4, label=f"{salt:g} mM")
               for salt in salts]
    fig.legend(handles=handles, title="NaCl concentration", loc="lower center",
               bbox_to_anchor=(0.535, 0.015), ncol=7, fontsize=7.5, title_fontsize=8,
               handlelength=1.5, columnspacing=1.0, handletextpad=0.5)
    pd.DataFrame(curves).to_csv(directory / "curves.csv", index=False)
    output = args.project / "SI/figures/surface_tension_langmuir_subset"
    pending = []
    for extension in ("pdf", "svg", "png"):
        temporary = output.with_suffix(".pending."+extension)
        fig.savefig(temporary, dpi=300, facecolor="white")
        pending.append((temporary, output.with_suffix("."+extension)))
    for temporary, destination in pending:
        temporary.replace(destination)
    plt.close(fig)
    print(f"Exported three-panel calibration figure: {len(observations)} observations, {len(parameters)} series.")


if __name__ == "__main__":
    main()
