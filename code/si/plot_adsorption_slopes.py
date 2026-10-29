"""Plot local adsorption slopes and their analytical dilute limits.

Contract: six panels match the Figure 4a-f parameter series; blue solid MA
and red dashed SL slopes approach distinct low-concentration plateaus.
Grey bands start at the second-lowest surface-tension concentration, assigning
adjacent-point slopes to their upper concentration. Top-edge ticks retain all
observations used for each SL fit, including the first point.
Use the archived parameter sets without refitting. Export editable PDF/SVG,
a 300-dpi preview, numerical curves, all 13 dilute limits and verification.
Concentration derivatives use mol/m3, so the displayed slope has units m.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import LogLocator, NullFormatter

PROJECT = Path(__file__).resolve().parents[2]
DATA = PROJECT / 'SI/source_data'
FIG = PROJECT / 'SI/figures'
PARAMS = PROJECT / 'figures/source_data/fig4/fig4_parameters.csv'
R, T, F = 8.314, 298.0, 96485.0
EPSILON = 8.85e-12 * 78.4
MW = {'PFOA': 413., 'PFHxA': 313., 'PFOS': 499.}
CASES = [('PFOA', 10), ('PFOA', 50), ('PFHxA', 100),
         ('PFHxA', 200), ('PFOS', 1), ('PFOS', 10)]
XMAX = {'PFOA': 1e4, 'PFHxA': 1e5, 'PFOS': 1e3}


def state(c, row):
    """Return Gamma and dGamma/dc from the implicit MA/GC solution."""
    k = row.K1_m3_mol
    gm = row.Gamma_m_MA_mol_m2
    total = row.I_mM + c
    bc = row.K2_m6_mol2 / k
    b = bc * total
    a = np.sqrt(8 * EPSILON * R * T * total) / (F * gm)

    def coverage(p):
        return a * np.sinh(p / 2) * (1 + b * np.exp(p))

    # p = -phi >= 0. The upper bracket stays strictly below saturation.
    pmax = brentq(lambda p: coverage(p) - (1 - 1e-12), 0, 30)

    def residual(p):
        th = coverage(p)
        return th - k * c * (np.exp(-p) + b) * (1 - th)

    p = brentq(residual, 0, pmax, xtol=1e-28, rtol=1e-14)
    th = coverage(p)
    tp = a * (.5 * np.cosh(p/2) * (1+b*np.exp(p))
              + np.sinh(p/2)*b*np.exp(p))
    u = np.exp(-p) + b
    # Figure 4 includes the counterions accompanying dissolved PFAS: c2=I+c.
    # Differentiate both the adsorption equation and GC charge balance.
    tc = th/(2*total) + a*np.sinh(p/2)*bc*np.exp(p)
    slope = gm*(1-th)*(tp*k*(u+c*bc)+tc*k*c*np.exp(-p)) / (
        tp*(1+k*c*u)+k*c*np.exp(-p)*(1-th))
    return gm * th, slope


def main():
    params = pd.read_csv(PARAMS)
    observations = pd.read_csv(DATA / 'szyszkowski_langmuir_calibration/observations.csv')
    limits = params[['PFAS', 'I_mM']].copy()
    limits['k_aw_MA_m'] = params.Gamma_m_MA_mol_m2 * (
        params.K1_m3_mol + params.K2_m6_mol2 * params.I_mM)
    limits['k_aw_SL_m'] = params.Gamma_m_SL_mol_m2 * params.K_SL_m3_mol
    limits['k_aw_SL_over_MA'] = limits.k_aw_SL_m / limits.k_aw_MA_m
    reference = pd.read_csv(DATA / 'empirical_dilute_limit_ratios.csv')
    merged = limits.merge(reference, on=['PFAS', 'I_mM'], validate='one_to_one')
    np.testing.assert_allclose(merged.k_aw_SL_over_MA,
        merged.GammaSL_over_GammaMA_min, rtol=1e-12)
    limits.to_csv(DATA / 'adsorption_slope_dilute_limits.csv', index=False)

    style = json.loads((DATA / 'style_inputs/figure_style.json').read_text())
    colors = [style['roles']['mass_action'], style['roles']['szyszkowski_langmuir']]
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8,
        'axes.labelsize': 8, 'axes.titlesize': 8.5, 'axes.linewidth': .7,
        'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
        'pdf.fonttype': 42, 'svg.fonttype': 'none',
        'mathtext.fontset': 'dejavusans', 'legend.frameon': False})
    fig, axes = plt.subplots(3, 2, figsize=(183/25.4, 174/25.4))
    fig.subplots_adjust(left=.13, right=.98, bottom=.09, top=.88,
                        hspace=.53, wspace=.26)
    curves, low_errors, fd_errors = [], [], []
    ranges, calibration_points = [], []
    archive = pd.read_csv(PROJECT / 'figures/source_data/fig4/fig4_a_f_model_curves.csv')
    archive_errors = []
    for index, (ax, (name, salt)) in enumerate(zip(axes.flat, CASES)):
        row = params[params.PFAS.eq(name) & params.I_mM.eq(salt)].iloc[0]
        obs = observations[observations.PFAS.eq(name) & observations.NaCl_mM.eq(salt)].copy()
        assert len(obs) > 0
        # Convert molar observations using this figure's existing mass-axis
        # convention; do not mix acid/salt molecular-weight conventions.
        obs['c_figure_mg_L'] = obs.c_mol_m3 * MW[name]
        concentrations = np.unique(obs.c_figure_mg_L)
        assert len(concentrations) >= 3
        lower, upper = concentrations[1], concentrations[-1]
        assert 1e-8 < lower < upper < XMAX[name]
        ranges.append({'PFAS': name, 'I_mM': salt, 'n_observations': len(obs),
            'c_min_mol_m3': obs.c_mol_m3.min(), 'c_max_mol_m3': obs.c_mol_m3.max(),
            'c_min_figure_mg_L': concentrations[0], 'c_max_figure_mg_L': upper,
            'shading_min_figure_mg_L': lower, 'shading_max_figure_mg_L': upper,
            'shading_min_mol_m3': lower / MW[name],
            'shading_rule': 'second-lowest distinct concentration to highest concentration',
            'figure_molecular_weight_g_mol': MW[name],
            'source_datasets': '; '.join(obs.dataset_id.unique()),
            'source_references': '; '.join(obs.reference.unique())})
        x = np.logspace(-8, np.log10(XMAX[name]), 651)
        c = x / MW[name]
        gamma, ma = np.array([state(v, row) for v in c]).T
        sl = row.Gamma_m_SL_mol_m2 * row.K_SL_m3_mol / (1 + row.K_SL_m3_mol*c)**2
        kma = row.Gamma_m_MA_mol_m2 * (row.K1_m3_mol + row.K2_m6_mol2*salt)
        ksl = row.Gamma_m_SL_mol_m2 * row.K_SL_m3_mol
        low_errors.extend([abs(ma[0]/kma - 1), abs(sl[0]/ksl - 1)])
        assert max(low_errors) < 1e-4
        for j in (0, 130, 325, 520, 650):
            h = c[j]*1e-4
            fd = (state(c[j]+h,row)[0]-state(c[j]-h,row)[0])/(2*h)
            fd_errors.append(abs(fd/ma[j]-1))
        ax.axvspan(lower, upper, color='#E7E7E7', lw=0, zorder=0)
        ax.plot(obs.c_figure_mg_L, np.full(len(obs), .985), '|',
                transform=ax.get_xaxis_transform(), color='#737373',
                markersize=4, markeredgewidth=.7, zorder=2)
        ax.plot(x, ma, color=colors[0], lw=1.5, zorder=3)
        ax.plot(x, sl, color=colors[1], lw=1.5, ls='--', zorder=3)
        for k, color, suffix in [(kma, colors[0], 'MA'), (ksl, colors[1], 'SL')]:
            value = f'{k:.2e}'.split('e')
            annotation = (rf'$k_{{\mathrm{{aw,{suffix}}}}} = '
                          rf'{value[0]}\times10^{{{int(value[1])}}}$ m')
            ax.text(.035, .32 if suffix == 'MA' else .20, annotation,
                    transform=ax.transAxes, color=color, fontsize=7.5, va='bottom')
        ax.set(xscale='log', yscale='log', xlim=(1e-8, XMAX[name]),
               title=f'{name}, {salt:g} mM NaCl')
        ax.set_xticks([1e-8, 1e-4, 1, 1e3 if name == 'PFOS' else 1e4])
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=5))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(which='major', direction='out', length=3, width=.7)
        ax.tick_params(which='minor', length=1.5, width=.5)
        ax.text(-.13, 1.06, 'abcdef'[index], transform=ax.transAxes,
                weight='bold', fontsize=10)
        if index % 2 == 0:
            ax.set_ylabel(r'$\mathrm{d}\Gamma/\mathrm{d}c_{1,0}$ (m)')
        if index >= 4:
            ax.set_xlabel(r'PFAS concentration (mg L$^{-1}$)')
        curves.append(pd.DataFrame({'PFAS': name, 'I_mM': salt,
            'c_mg_L': x, 'c_mol_m3': c, 'Gamma_MA_mol_m2': gamma,
            'slope_MA_m': ma, 'slope_SL_m': sl}))
        obs_gamma, obs_ma = np.array([state(v, row) for v in obs.c_mol_m3]).T
        obs['Gamma_MA_mol_m2'] = obs_gamma
        obs['Gamma_SL_mol_m2'] = row.Gamma_m_SL_mol_m2 * row.K_SL_m3_mol * obs.c_mol_m3 / (1 + row.K_SL_m3_mol * obs.c_mol_m3)
        obs['slope_MA_m'] = obs_ma
        obs['slope_SL_m'] = ksl / (1 + row.K_SL_m3_mol * obs.c_mol_m3)**2
        obs['Gamma_SL_over_MA'] = obs.Gamma_SL_mol_m2 / obs.Gamma_MA_mol_m2
        obs['slope_SL_over_MA'] = obs.slope_SL_m / obs.slope_MA_m
        calibration_points.append(obs)
        # Independent numerical agreement with recorded Figure 4 surface excesses.
        sub = archive[archive.PFAS.eq(name) & archive.I_mM.eq(salt)].iloc[::80]
        for point in sub.itertuples():
            value = state(point.Concentration_mg_L/MW[name], row)[0]
            archive_errors.append(abs(value/(point.Gamma_MA_umol_m2*1e-6)-1))
    assert max(fd_errors) < 1e-5
    assert max(archive_errors) < 1e-4
    for left, right in axes:
        lo = min(left.get_ylim()[0], right.get_ylim()[0])
        hi = max(left.get_ylim()[1], right.get_ylim()[1])
        left.set_ylim(lo, hi*2)
        right.set_ylim(lo, hi*2)
    handles = [Line2D([], [], color=colors[0], lw=1.5, label='Mass-action model'),
        Line2D([], [], color=colors[1], ls='--', lw=1.5,
               label='Szyszkowski–Langmuir isotherm')]
    fig.legend(handles=handles, ncol=2, loc='upper center',
               bbox_to_anchor=(.55, .995), fontsize=8.5)
    fig.legend(handles=[Patch(facecolor='#E7E7E7', edgecolor='none',
               label='Surface-tension data range')],
               loc='upper center', bbox_to_anchor=(.55, .952), fontsize=8)
    for ext in ('pdf','svg','png'):
        fig.savefig(FIG / f'adsorption_slopes.{ext}', dpi=300, facecolor='white')
    plt.close(fig)
    pd.concat(curves).to_csv(DATA / 'adsorption_slope_curves.csv', index=False)
    pd.DataFrame(ranges).to_csv(DATA / 'adsorption_slope_calibration_ranges.csv', index=False)
    pd.concat(calibration_points).to_csv(DATA / 'adsorption_slope_calibration_observations.csv', index=False)
    audit = {'parameters_refitted': False, 'temperature_K': T,
        'maximum_low_concentration_relative_error': max(low_errors),
        'maximum_finite_difference_relative_error': max(fd_errors),
        'maximum_archived_Gamma_relative_difference': max(archive_errors),
        'panels': 6, 'curves': 12, 'points_per_curve': 651,
        'calibration_observations_shown': sum(r['n_observations'] for r in ranges),
        'shading_source': 'SI/source_data/szyszkowski_langmuir_calibration/observations.csv',
        'shading_definition': 'second-lowest to highest distinct retained surface-tension concentration; adjacent-point slopes assigned to upper concentration',
        'first_observation_retained_in_fit_and_top_ticks': True,
        'mass_axis_conversion': 'molar observations times the existing anion molecular weight',
        'x_limits_by_compound_mg_L': {name: [1e-8, xmax] for name, xmax in XMAX.items()},
        'dilute_limit_series': len(limits),
        'source_parameter_file': 'figures/source_data/fig4/fig4_parameters.csv'}
    (DATA / 'adsorption_slopes_verification.json').write_text(json.dumps(audit,indent=2))
    print(json.dumps(audit,indent=2))
    print(limits[limits.I_mM.isin([1,10,50,100,200])].to_string(index=False))


if __name__ == '__main__':
    main()
