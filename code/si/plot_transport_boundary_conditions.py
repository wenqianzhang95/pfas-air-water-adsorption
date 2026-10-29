"""Draw the transport domain, boundaries, and prescribed initial profiles."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch


def main():
    project = Path(__file__).resolve().parents[2]
    plt.rcParams.update({"font.family": "Arial", "font.size": 9,
                         "pdf.fonttype": 42, "svg.fonttype": "none",
                         "mathtext.fontset": "dejavusans"})
    fig = plt.figure(figsize=(183/25.4, 183/25.4))
    ax = fig.add_axes([0, .475, 1, .525])
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    blue, green, ink = "#355F8D", "#2A9D8F", "#282D32"
    def text(x, y, label, **kwargs):
        kwargs.setdefault("color", ink)
        ax.text(x, y, label, va="top", **kwargs)
    def arrow(start, end, color=blue, **kwargs):
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
                                   mutation_scale=12, linewidth=1.5,
                                   color=color, **kwargs))
    # Soil geometry is schematic; no unmeasured profile or water table is drawn.
    left, right, bottom, top = .38, .61, .23, .75
    ax.add_patch(Rectangle((left, bottom), right-left, top-bottom,
                           facecolor="#F3F1E9", edgecolor=ink, linewidth=1.1))
    ax.plot([left-.025, right+.025], [top, top], color=ink, lw=1.3)
    ax.plot([left-.025, right+.025], [bottom, bottom], color=ink, lw=1.1)
    text((left+right)/2, .66, "Unsaturated\nsoil profile", ha="center",
         fontsize=10, linespacing=1.4)
    arrow((.495, .53), (.495, .39), color=green)
    text(.495, .35, "PFAS transport", ha="center", fontsize=8.5)
    ax.annotate("", xy=(.345, bottom), xytext=(.345, top),
                arrowprops={"arrowstyle": "<->", "lw": .9, "color": ink})
    ax.text(.32, (bottom+top)/2, "1.5 m", rotation=90, ha="center", va="center", color=ink)
    text(left+.015, top-.015, "$z=0$", ha="left", fontsize=8)
    text(left+.015, bottom+.04, "$z=1.5$ m", ha="left", fontsize=8)
    for x in (.435, .495, .555):
        arrow((x, .87), (x, top+.01))
    text(.495, .975, "Rainfall-driven infiltration", ha="center", color=blue,
         fontsize=10, fontweight="bold")
    text(.495, .925, "$q_{\\mathrm{in}}(t)=0.47P(t)$", ha="center", fontsize=10)
    text(.035, .78, "Upper boundary", fontweight="bold", fontsize=10)
    text(.035, .715, "Water input: 0.47 × daily rainfall\nNo separate evapotranspiration term\nNo added AFFF water", fontsize=8.5, linespacing=1.65)
    text(.035, .54, "No further PFAS input\n$J_{\\mathrm{PFAS}}=0$", fontsize=9,
         linespacing=1.6)
    text(.035, .39, "If ponding develops:\nswitch to a prescribed\npressure-head boundary", fontsize=8.2,
         linespacing=1.45)
    text(.67, .76, "Initial state", fontweight="bold", fontsize=10)
    text(.67, .69, "Water pressure head from\nthe site water-content profile", fontsize=8.5,
         linespacing=1.5)
    text(.67, .57, "Measured total soil PFAS profile\nInitial aqueous concentration\nfrom local adsorption equilibrium", fontsize=8.5, linespacing=1.5)
    text(.67, .405, "Lower boundary", fontweight="bold", fontsize=10)
    text(.67, .34, "Free drainage: $\\partial h/\\partial z=0$\nPFAS: $\\partial c/\\partial z=0$\nAdvective PFAS outflow: $J_{\\mathrm{PFAS}}=qc$",
         fontsize=8.5, linespacing=1.7)
    arrow((.46, bottom-.005), (.46, .09), color=blue)
    arrow((.54, bottom-.005), (.54, .09), color=green)
    text(.495, .065, "Water and dissolved PFAS leave the model", ha="center", fontsize=8.5)
    fig.text(.025, .99, "a", fontweight="bold", fontsize=11, va="top", color=ink)

    # The water values are prescribed inputs used to reconstruct pressure head.
    # They are not the isotherm-dependent water contents after initialization.
    with (project / "SI/source_data/transport_initial_water_contents.csv").open(newline="") as handle:
        water = list(csv.DictReader(handle))
    with (project / "SI/source_data/transport_initial_soil_pfas_profiles.csv").open(newline="") as handle:
        soil = list(csv.DictReader(handle))
    ax_b = fig.add_axes([.12, .085, .30, .33])
    ax_c = fig.add_axes([.59, .085, .36, .33])
    for axis in (ax_b, ax_c):
        axis.set_ylim(1.5, 0)
        axis.set_yticks([0, .5, 1, 1.5])
        axis.set_ylabel("Depth (m)", fontsize=9)
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(labelsize=8, direction="out", width=.7, length=3)
        axis.spines[["bottom", "left"]].set_linewidth(.8)
        axis.set_axisbelow(True)
        axis.grid(axis="y", color="#E5E7E9", linewidth=.5)
    ax_b.plot([float(row["water_content_full_precision"]) for row in water],
              [float(row["depth_cm"])/100 for row in water],
              linestyle="none", marker="o", markersize=5, color=blue, clip_on=False)
    ax_b.set_xlim(.09, .15)
    ax_b.set_xticks([.10, .12, .14])
    ax_b.set_xlabel(r"Water content, $\theta_w$ (m$^3$ m$^{-3}$)", fontsize=8.5)
    ax_b.set_title("Initial water-content inputs", fontsize=9, pad=9)
    for pfas, color, style in (("PFOA", "#6E9D88", "-"),
                               ("PFOS", "#D0A25D", "-")):
        rows = sorted((row for row in soil if row["PFAS"] == pfas),
                      key=lambda row: float(row["Depth_cm"]))
        assert len(rows) == 150, f"Incomplete {pfas} initial profile"
        ax_c.plot([float(row["TotalSoil_ug_kg"]) for row in rows],
                  [float(row["Depth_cm"])/100 for row in rows],
                  color=color, linestyle=style, linewidth=1.6, label=pfas)
    ax_c.set_xscale("log")
    ax_c.set_xlim(1e-2, 2e4)
    ax_c.set_xticks([1e-2, 1, 1e2, 1e4])
    ax_c.set_xlabel(r"Total soil PFAS ($\mu$g kg$^{-1}$)", fontsize=8.5)
    ax_c.set_title("Initial PFAS concentration", fontsize=9, pad=9)
    ax_c.legend(loc="upper left", frameon=False, fontsize=8, handlelength=2)
    fig.text(.025, .455, "b", fontweight="bold", fontsize=11, va="top", color=ink)
    fig.text(.50, .455, "c", fontweight="bold", fontsize=11, va="top", color=ink)
    output = project / "SI/figures/transport_boundary_conditions"
    for suffix in ("pdf", "svg", "png"):
        temp = output.with_suffix(".pending."+suffix)
        fig.savefig(temp, dpi=300, facecolor="white")
        temp.replace(output.with_suffix("."+suffix))
    plt.close(fig)


if __name__ == "__main__":
    main()
