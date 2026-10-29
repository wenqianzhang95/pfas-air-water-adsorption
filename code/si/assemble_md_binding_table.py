"""Tabulate archived MD plot inputs without recalculating or fitting results."""
from pathlib import Path
import csv
import hashlib
import json
import math

PROJECT = Path(__file__).resolve().parents[2]
DATA = PROJECT / 'SI' / 'source_data'
INPUT = DATA / 'md_binding_inputs'
FIELDS = ['isFiltered', 'gamma', 'c2s', 'brMD', 'c2sSTD', 'STD']
COMPOUNDS = ['PFOA', 'PFHxA', 'PFOS']
NACL_DISPLAY = {499: 500, 997: 1000}


def close(a, b):
    return math.isclose(float(a), float(b), rel_tol=1e-10, abs_tol=1e-9)


def display(value):
    if value == 0:
        return '0'
    places = 2 - math.floor(math.log10(abs(value)))
    return f'{round(value, places):.{max(0, places)}f}'


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    with (INPUT / 'main_binding_plot_points.csv').open(newline='', encoding='utf-8-sig') as handle:
        plotted = list(csv.DictReader(handle))
    records = []
    verified = 0
    for compound in COMPOUNDS:
        source = INPUT / f'{compound}_plotResults_panels_data.json'
        data = json.loads(source.read_text(encoding='utf-8'))
        assert len({len(data[key]) for key in FIELDS}) == 1
        expected = [row for row in plotted if row['Compound'] == compound]
        plotted_index = 0
        for index, values in enumerate(zip(*(data[key] for key in FIELDS)), start=1):
            salt, gamma, c2s, beta, c2err, betaerr = values
            excluded = compound == 'PFOS' and 600 < c2s < 900 and beta < 40 and betaerr > 40
            row = {
                'PFAS': compound, 'NaCl_mM': salt,
                'Gamma_umol_m2': gamma * 1e6,
                'c2s_mM': c2s, 'c2s_reported_error_mM': c2err,
                'beta_percent': beta, 'beta_reported_error_percentage_points': betaerr,
                'source_point_index': index,
                'main_figure_2_point_index': '',
                'status': 'excluded_from_current_figure' if excluded else 'plotted',
                'source_file': source.name,
            }
            if not excluded:
                item = expected[plotted_index]
                for actual, key in [(salt, 'IonicStrength_mM'), (c2s, 'X'),
                                    (beta, 'Y'), (c2err, 'X_error'), (betaerr, 'Y_error')]:
                    assert close(actual, item[key]), (compound, index, key)
                plotted_index += 1
                row['main_figure_2_point_index'] = plotted_index
                verified += 1
            records.append(row)
        assert plotted_index == len(expected)
    assert verified == len(plotted) == 37
    write_csv(DATA / 'md_counterion_binding_source_records.csv', records)

    unique = {}
    for row in records:
        if row['status'] != 'plotted':
            continue
        key = tuple(row[name] for name in ['PFAS', 'NaCl_mM', 'Gamma_umol_m2',
                                         'c2s_mM', 'c2s_reported_error_mM',
                                         'beta_percent', 'beta_reported_error_percentage_points'])
        if key in unique:
            old = unique[key]
            old['source_point_index'] = f"{old['source_point_index']};{row['source_point_index']}"
            old['main_figure_2_point_index'] = f"{old['main_figure_2_point_index']};{row['main_figure_2_point_index']}"
        else:
            unique[key] = dict(row)
    rows = sorted(unique.values(), key=lambda r: (COMPOUNDS.index(r['PFAS']), r['NaCl_mM'], r['Gamma_umol_m2']))
    assert len(rows) == 36
    for row in rows:
        row['NaCl_display_mM'] = NACL_DISPLAY.get(row['NaCl_mM'], row['NaCl_mM'])

    pending = []
    figure1_file = DATA / 'md_binding_figure1_digitized.csv'
    with figure1_file.open(newline='', encoding='utf-8') as handle:
        figure1 = list(csv.DictReader(handle))
    figure1_lookup = {(r['PFAS'], int(r['NaCl_label_mM']), float(r['Gamma_label_umol_m2'])): r
                      for r in figure1}
    assert len(figure1_lookup) == 15
    for compound in COMPOUNDS:
        for salt, panel, gammas in [(100, '1g', [0.04, 0.26, 0.47]),
                                    (20, '1h', [0.03, 0.09])]:
            for gamma in gammas:
                digitized = figure1_lookup[(compound, salt, gamma)]
                pending.append({'PFAS': compound, 'NaCl_label_mM': salt,
                                'Gamma_label_umol_m2': gamma, 'main_figure_panel': panel,
                                'beta_percent': '', 'beta_reported_error_percentage_points': '',
                                'digitized_beta_percent': digitized['beta_vector_estimate_percent'],
                                'digitized_error_bar_percentage_points': digitized['beta_error_bar_estimate_percentage_points'],
                                'status': 'figure_digitized_original_data_pending'})
    # Reconcile by compound, electrolyte condition and the plotted mean, not
    # by the nominal surface-excess labels alone. Figure 1 rounds some means
    # to one decimal place. Keep the available numerical inputs in the table;
    # preserve all figure readings, including differing errors, in the CSVs.
    additional = []
    reconciliation = []
    for row in pending:
        matches = [source for source in rows
                   if source['PFAS'] == row['PFAS']
                   and source['NaCl_display_mM'] == row['NaCl_label_mM']
                   and abs(source['beta_percent'] - float(row['digitized_beta_percent'])) < 0.05]
        assert len(matches) <= 1, ('ambiguous match', row)
        source = matches[0] if matches else None
        if source is None:
            status = 'additional_figure1_record'
            additional.append(row)
        elif abs(source['Gamma_umol_m2'] - row['Gamma_label_umol_m2']) <= 0.01:
            status = 'duplicate_approximate_gamma_label'
        else:
            status = 'duplicate_mean_gamma_label_conflict'
        row['table_status'] = status
        row['matched_source_point_index'] = source['source_point_index'] if source else ''
        reconciliation.append({
            'PFAS': row['PFAS'], 'NaCl_label_mM': row['NaCl_label_mM'],
            'main_figure_panel': row['main_figure_panel'],
            'figure1_Gamma_label_umol_m2': row['Gamma_label_umol_m2'],
            'figure1_beta_percent': row['digitized_beta_percent'],
            'figure1_error_bar_percentage_points': row['digitized_error_bar_percentage_points'],
            'table_status': status,
            'retained_source_point_index': source['source_point_index'] if source else '',
            'retained_Gamma_umol_m2': source['Gamma_umol_m2'] if source else '',
            'retained_beta_percent': source['beta_percent'] if source else '',
            'retained_error_bar_percentage_points': source['beta_reported_error_percentage_points'] if source else '',
        })
    assert len(additional) == 7
    assert sum(r['table_status'] == 'duplicate_approximate_gamma_label' for r in pending) == 5
    assert sum(r['table_status'] == 'duplicate_mean_gamma_label_conflict' for r in pending) == 3
    write_csv(DATA / 'md_counterion_binding_data.csv', rows)
    write_csv(DATA / 'md_binding_figure1_pending.csv', pending)
    write_csv(DATA / 'md_binding_cross_figure_reconciliation.csv', reconciliation)

    header = r'''PFAS & NaCl & $\Gamma$ & $c_{2,\mathrm{s}}$ & $\beta$ \\
 & (mM) & ($\mu$mol\,m$^{-2}$) & (mM) & (\%) \\'''
    text = r'''% Generated by SI/code/assemble_md_binding_table.py.
\begingroup
\small
\setstretch{1.0}
\begin{longtable}{llrrr}
\caption{MD simulation cases and calculated counterion binding ratios. $\Gamma$ is the PFAS surface excess, $c_{2,\mathrm{s}}$ is the local counterion concentration adjacent to the interface, and $\beta$ is the ensemble-averaged counterion binding ratio. Values after $\pm$ denote standard deviations. A dash indicates missing data.}
\label{tab:md_binding_data_si}\\
\toprule
''' + header + r'''
\midrule
\endfirsthead
\multicolumn{5}{l}{\textbf{Table~\thetable\ continued}}\\
\toprule
''' + header + r'''
\midrule
\endhead
\midrule
\multicolumn{5}{r}{Continued on next page}\\
\endfoot
\bottomrule
\endlastfoot
'''
    table_rows = []
    for row in rows:
        key = (COMPOUNDS.index(row['PFAS']), row['NaCl_display_mM'], row['Gamma_umol_m2'])
        rendered = (f"{row['PFAS']} & {row['NaCl_display_mM']:,.0f} & ${row['Gamma_umol_m2']:.2f}$ & "
                    f"${display(row['c2s_mM'])} \\pm {display(row['c2s_reported_error_mM'])}$ & "
                    f"${display(row['beta_percent'])} \\pm {display(row['beta_reported_error_percentage_points'])}$ "
                    + r'\\' + '\n')
        table_rows.append((key, rendered))
    for row in additional:
        digitized = figure1_lookup[(row['PFAS'], row['NaCl_label_mM'], row['Gamma_label_umol_m2'])]
        key = (COMPOUNDS.index(row['PFAS']), row['NaCl_label_mM'], row['Gamma_label_umol_m2'])
        rendered = (f"{row['PFAS']} & {row['NaCl_label_mM']:,} & ${row['Gamma_label_umol_m2']:.2f}$ & --- & "
                    + f"${digitized['beta_display_percent']} \\pm {digitized['beta_error_bar_display_percentage_points']}$ "
                    + r'\\' + '\n')
        table_rows.append((key, rendered))
    assert len(table_rows) == 43
    previous = None
    for key, rendered in sorted(table_rows):
        if previous is not None and previous != key[0]:
            text += '\\midrule\n'
        previous = key[0]
        text += rendered
    text += '\\end{longtable}\n\\endgroup\n'
    tables = PROJECT / 'SI' / 'tables'
    tables.mkdir(exist_ok=True)
    (tables / 'md_counterion_binding_data.tex').write_text(text, encoding='utf-8')
    report = {
        'source_records': len(records), 'plotted_records_verified': verified,
        'distinct_plotted_records': len(rows), 'duplicate_plot_entries': verified - len(rows),
        'pending_figure1_records': len(pending),
        'figure1_digitized_records': len(figure1),
        'cross_figure_duplicate_records': len(pending) - len(additional),
        'cross_figure_gamma_label_conflicts': 3,
        'additional_figure1_table_records': len(additional),
        'total_table_records': len(rows) + len(additional),
        'figure1_digitized_csv_sha256': hashlib.sha256(figure1_file.read_bytes()).hexdigest(),
        'source_hashes': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted(INPUT.iterdir())
                          if p.name.endswith('_plotResults_panels_data.json')
                          or p.name == 'main_binding_plot_points.csv'},
        'error_convention': 'Author specifies ensemble-averaged binding ratios and standard deviations. Input and figure-extracted error magnitudes retained without conversion.',
    }
    (DATA / 'md_binding_inputs' / 'verification.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'source_hashes'}))


if __name__ == '__main__':
    main()
