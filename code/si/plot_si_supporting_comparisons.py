"""Restyle supplementary comparisons from recorded data, without refitting."""
from pathlib import Path
import json
import importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter

PROJECT = Path(__file__).resolve().parents[2]
DATA = PROJECT / "SI/source_data"
INPUT = DATA / "style_inputs"
STYLE = json.loads((INPUT / "figure_style.json").read_text())
COL = STYLE["roles"]
INK = STYLE["ink"]
AUDIT = {}


def theme():
    plt.rcParams.update({"font.family": "Arial", "font.size": 8.5,
        "axes.labelsize": 8.5, "axes.titlesize": 9, "axes.labelcolor": INK,
        "text.color": INK, "xtick.color": INK, "ytick.color": INK,
        "axes.edgecolor": INK, "axes.linewidth": .75,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "legend.frameon": False,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "mathtext.fontset": "dejavusans", "lines.linewidth": 1.5})


def axis_style(ax, letter=None):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(which="major", direction="out", length=3, width=.75)
    ax.tick_params(which="minor", length=1.7, width=.5)
    if letter:
        ax.text(-.14, 1.08, letter, transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")


def export(fig, name):
    output = PROJECT / "SI/figures" / name
    pending = []
    for ext in ["pdf", "svg", "png"]:
        temp = output.with_suffix(".pending."+ext)
        fig.savefig(temp, dpi=300, facecolor="white")
        pending.append((temp, output.with_suffix("."+ext)))
    for temp, target in pending:
        temp.replace(target)
    plt.close(fig)


def surface_tension_validation():
    model = pd.read_csv(INPUT / "s2_model.csv")
    obs = pd.read_csv(INPUT / "s2_observations.csv")
    fig, ax = plt.subplots(figsize=(150/25.4, 110/25.4))
    fig.subplots_adjust(left=.14, right=.97, bottom=.17, top=.87)
    ax.plot(model.c_mg_L, model.sigma_mN_m, color=COL["mass_action"], label="Mass-action model")
    ax.scatter(obs.c_mg_L, obs.sigma_mN_m, marker="s", s=20,
               facecolors="white", edgecolors=COL["experiment"], linewidths=.8,
               label="Measurements", zorder=3)
    ax.set(xscale="log", xlim=(1e-2,1e4), ylim=(22,80),
           xlabel=r"PFOA concentration (mg L$^{-1}$)",
           ylabel=r"Surface tension (mN m$^{-1}$)", title="PFOA, 10 mM NaCl")
    ax.set_xticks([1e-2,1e-1,1,10,100,1000,10000])
    ax.legend(loc="lower left", handlelength=2)
    axis_style(ax)
    AUDIT["S2"]={"observations":len(obs),"model_points":len(model)}
    export(fig,"PFOA_10mM_ST")


def direct_comparisons():
    data = pd.read_csv(DATA / "surface_excess_error_pointwise_current.csv")
    group = data[data.PFAS.eq("PFOA") & data.NaCl_mM.eq(100)].sort_values("concentration_mol_m3")
    x=group.concentration_mol_m3.to_numpy()*414
    fig, ax=plt.subplots(figsize=(150/25.4,110/25.4))
    fig.subplots_adjust(left=.15,right=.97,bottom=.18,top=.91)
    ax.scatter(x,group.Gamma_observed_umol_m2,s=23,facecolors="white",
               edgecolors=INK,linewidths=.8,label="Measurements",zorder=3)
    ax.plot(x,group.Gamma_MA_umol_m2,color=COL["mass_action"],label="Mass-action model")
    ax.plot(x,group.Gamma_SL_umol_m2,color=COL["szyszkowski_langmuir"],ls="--",label="Szyszkowski–Langmuir isotherm")
    ax.set(xscale="log",yscale="log",xlim=(x.min()*.9,x.max()*1.1),ylim=(1e-5,1e-2),
        xlabel=r"PFOA concentration (mg L$^{-1}$)",ylabel=r"Surface excess ($\mu$mol m$^{-2}$)")
    ax.legend(loc="upper left",handlelength=2,fontsize=7.8)
    axis_style(ax)
    export(fig,"PFOA_Fit_I100")
    fig, axes=plt.subplots(1,2,figsize=(183/25.4,88/25.4))
    fig.subplots_adjust(left=.09,right=.98,bottom=.30,top=.85,wspace=.29)
    for i,(ax,salt) in enumerate(zip(axes,[10,100])):
        sub=data[data.PFAS.eq("PFOA") & data.NaCl_mM.eq(salt)].sort_values("concentration_mol_m3")
        x=sub.concentration_mol_m3.to_numpy()*414
        for code,color,ls,marker in [("MA",COL["mass_action"],"-","o"),("SL",COL["szyszkowski_langmuir"],"--","s")]:
            y=sub[f"squared_log10_error_{code}"]
            expected=np.log10(sub[f"Gamma_pred_over_obs_{code}"])**2
            np.testing.assert_allclose(y,expected,rtol=1e-12,atol=1e-14)
            ax.plot(x,y,color=color,ls=ls,marker=marker,ms=3,mew=.6)
        ax.set(xscale="log",xlim=(1e-4,1e-1) if salt==10 else (5e-7,5e-3),
               ylim=(0,1) if salt==10 else (0,1.4),title=f"{salt} mM NaCl",
               xlabel=r"PFOA concentration (mg L$^{-1}$)")
        ax.xaxis.set_major_locator(LogLocator(base=10,numticks=5))
        axis_style(ax,"ab"[i])
        AUDIT[f"S5_{salt}"]={"points_per_model":len(sub)}
    axes[0].set_ylabel("Squared logarithmic error")
    handles=[Line2D([],[],color=COL["mass_action"],marker="o",ms=3,label="Mass-action model"),
             Line2D([],[],color=COL["szyszkowski_langmuir"],ls="--",marker="s",ms=3,label="Szyszkowski–Langmuir isotherm")]
    fig.legend(handles=handles,ncol=2,loc="lower center",bbox_to_anchor=(.54,.04),columnspacing=2)
    export(fig,"PFOA_Error_Comparison")
    AUDIT["S4"]={"points_per_series":len(group),"predictions":"unchanged recorded pointwise values"}


def dilute_refits():
    params=pd.read_csv(DATA/"low_concentration_langmuir_refit_summary.csv")
    obs=pd.read_csv(DATA/"szyszkowski_langmuir_calibration/observations.csv")
    fig,axes=plt.subplots(1,2,figsize=(183/25.4,91/25.4),sharey=True)
    fig.subplots_adjust(left=.09,right=.98,bottom=.32,top=.85,wspace=.18)
    for i,(ax,(_,row)) in enumerate(zip(axes,params.iterrows())):
        mw=414 if row.PFAS=="PFOA" else 314
        x=np.logspace(-8 if i==0 else -7,2,72)
        points=obs[obs.PFAS.eq(row.PFAS)&obs.NaCl_mM.eq(row.NaCl_mM)]
        keep=points.c_mol_m3>=1/(10*row.K_L_full_m3_mol)
        assert int((~keep).sum())==row.n_removed
        for version,color,ls in [("full",COL["full_data_fit"],"-"),("trim",COL["trimmed_data_fit"],"--")]:
            k=row[f"K_L_{version}_m3_mol"];gm=row[f"Gamma_m_{version}_umol_m2"]*1e-6
            sigma=71-1000*8.314*298*gm*np.log1p(k*x)
            ax.plot(x*mw,sigma,color=color,ls=ls,lw=1.6)
        for retain in [True,False]:
            sub=points[keep if retain else ~keep]
            ax.scatter(sub.c_mg_L,sub.sigma_mN_m,s=20,marker="s" if i==0 else "o",
                facecolors=COL["full_data_fit"] if retain else "white",
                edgecolors=COL["full_data_fit"],linewidths=.8,zorder=3,clip_on=False)
        ax.set(xscale="log",xlim=(1e-2,1e5) if i==0 else (1e-1,1e5),ylim=(30,80),
            xlabel=f"{row.PFAS} concentration (mg L$^{{-1}}$)",title=f"{row.PFAS}, {row.NaCl_mM:g} mM NaCl")
        ax.set_xticks([1e-2,1,1e2,1e4] if i==0 else [1e-1,1e1,1e3,1e5])
        axis_style(ax,"ab"[i])
    axes[0].set_ylabel(r"Surface tension (mN m$^{-1}$)")
    handles=[Line2D([],[],ls="none",marker="o",color=COL["full_data_fit"],ms=4,label="Retained observations"),
        Line2D([],[],ls="none",marker="o",mfc="white",color=COL["full_data_fit"],ms=4,label="Excluded observations"),
        Line2D([],[],color=COL["full_data_fit"],label="Full-data fit"),
        Line2D([],[],color=COL["trimmed_data_fit"],ls="--",label="Fit with dilute points excluded")]
    fig.legend(handles=handles,ncol=2,loc="lower center",bbox_to_anchor=(.54,.025),columnspacing=2)
    AUDIT["S6"]={"refitting":False,"parameter_source":"recorded full and trimmed fits"}
    export(fig,"surface_tension_langmuir_subset_trimmed")


def tension_potential():
    observations=pd.read_csv(DATA/"szyszkowski_langmuir_calibration/observations.csv")
    fig,axes=plt.subplots(1,2,figsize=(183/25.4,94/25.4))
    fig.subplots_adjust(left=.09,right=.90,bottom=.32,top=.86,wspace=.58)
    for i,(ax,salt,key,k,gm,limits) in enumerate(zip(axes,[10,100],["data","data100"],
            [6.2277,25.4063],[2.9889e-6,3.2167e-6],[(40,72),(20,80)])):
        cache=pd.read_csv(INPUT/f"s7_{key}.csv")
        x=cache.C_R_mM
        ax.plot(x,71-1000*8.314*298*gm*np.log1p(k*x),color=COL["szyszkowski_langmuir"])
        obs=observations[observations.PFAS.eq("PFOA")&observations.NaCl_mM.eq(salt)]
        ax.scatter(obs.c_mol_m3,obs.sigma_mN_m,s=16,facecolors="white",
            edgecolors=COL["szyszkowski_langmuir"],linewidths=.8,zorder=3)
        ax.axvline(1/(10*k),color="#81909A",ls="--",lw=.9)
        ax.set(xscale="log",xlim=(1e-6,1e2),ylim=limits,title=f"{salt} mM NaCl",
            xlabel="PFOA concentration (mM)")
        ax.set_xticks([1e-6,1e-4,1e-2,1,100])
        ax.set_ylabel(r"Surface tension (mN m$^{-1}$)",color=COL["szyszkowski_langmuir"])
        ax.tick_params(axis="y",colors=COL["szyszkowski_langmuir"])
        axis_style(ax,"ab"[i])
        right=ax.twinx()
        phi=96485*cache.surface_potential_V/(8.314*298)
        right.plot(x,phi,color=COL["potential"],lw=1.5)
        right.set_ylim((-5,0) if i==0 else (-3,0))
        right.spines["top"].set_visible(False)
        right.spines["right"].set_color(COL["potential"])
        right.tick_params(axis="y",colors=COL["potential"],length=3,width=.7)
        right.set_ylabel(r"Dimensionless potential, $\phi$",color=COL["potential"],labelpad=6)
        AUDIT[f"S7_{salt}"]={"potential_points":len(cache),"potential":"unchanged saved arrays"}
    handles=[Line2D([],[],color=COL["szyszkowski_langmuir"],label="Surface tension: Szyszkowski–Langmuir isotherm"),
        Line2D([],[],color=COL["szyszkowski_langmuir"],ls="none",marker="o",mfc="white",ms=3,label="Surface tension: measurements"),
        Line2D([],[],color=COL["potential"],label="Interfacial potential"),
        Line2D([],[],color="#81909A",ls="--",label=r"$c_{1,0}=1/(10K_{\mathrm{SL}})$")]
    fig.legend(handles=handles,ncol=2,loc="lower center",bbox_to_anchor=(.51,.025),fontsize=7.5,columnspacing=1.6)
    export(fig,"surface_tension_and_dimensionless_potential_comparison_10mM_100mM")


def threshold_map():
    file=Path(__file__).with_name("analyze_regime_threshold_sensitivity.py")
    spec=importlib.util.spec_from_file_location("threshold",file)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    boundaries=pd.read_csv(DATA/"regime_threshold_sensitivity/threshold_regime_boundaries.csv")
    salts=np.sort(boundaries.NaCl_mM.unique());concs=np.logspace(-4,4,700)
    plots=[]
    for name in module.PARAMETERS:
        le,ls,residual=module.fields(name,salts,concs)
        ceq=boundaries[boundaries.PFAS.eq(name)&boundaries.Tolerance.eq(.1)].sort_values("NaCl_mM").Equal_penalty_concentration_mg_L.to_numpy()
        plots.append((name,module.classify(le,ls,.1),ceq))
    module.make_figure(PROJECT,salts,concs,plots,boundaries)
    AUDIT["S9"]={"boundaries":"unchanged archived values","classification":"unchanged equations and grid"}


if __name__=="__main__":
    theme()
    surface_tension_validation()
    direct_comparisons()
    dilute_refits()
    tension_potential()
    threshold_map()
    (INPUT/"style_verification.json").write_text(json.dumps(AUDIT,indent=2))
    print(json.dumps(AUDIT,indent=2))
