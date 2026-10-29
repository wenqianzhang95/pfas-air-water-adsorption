"""Recover provisional plotted values from Figure 1 vector geometry.

This recovers what is displayed in the PDF, not original simulation output.
"""
from pathlib import Path
import csv
import hashlib
import json
import math

import pymupdf

PROJECT = Path(__file__).resolve().parents[2]
PDF = PROJECT / 'figures' / 'Fig1_surface_excess_A4.pdf'
DATA = PROJECT / 'SI' / 'source_data'
QA = PROJECT.parent / 'md_binding_table_revision_20261021'
COLORS = {'PFOA': (0.3333, 0.4902, 0.6039),
          'PFHxA': (0.4314, 0.6157, 0.5333),
          'PFOS': (0.8157, 0.6353, 0.3647)}


def same_color(actual, target):
    return actual is not None and max(abs(a - b) for a, b in zip(actual, target)) < 0.0002


def main():
    document = pymupdf.open(PDF)
    page = document[0]
    drawings = page.get_drawings()
    axes = sorted([d['rect'] for d in drawings
                   if same_color(d['fill'], (1, 1, 1))
                   and 0.29 * page.rect.width < d['rect'].width < 0.32 * page.rect.width
                   and 100 < d['rect'].height < 140
                   and d['rect'].x0 > 0.6 * page.rect.width], key=lambda r: r.y0)
    assert len(axes) == 4, axes
    source_hash = hashlib.sha256(PDF.read_bytes()).hexdigest()
    rows = []
    audit = []
    for axis, panel, salt, maximum, gammas in [
            (axes[0], '1g', 100, 100, [0.04, 0.26, 0.47]),
            (axes[1], '1h', 20, 50, [0.03, 0.09])]:
        scale = maximum / axis.height
        for compound, color in COLORS.items():
            bars = sorted([(i, d) for i, d in enumerate(drawings)
                           if same_color(d['fill'], color)
                           and axis.contains(d['rect'])
                           and abs(d['rect'].y1 - axis.y1) < 0.002
                           and len(d['items']) == 4], key=lambda item: item[1]['rect'].x0)
            assert len(bars) == len(gammas), (panel, compound)
            for gamma, (drawing_index, bar) in zip(gammas, bars):
                box = bar['rect']
                center = (box.x0 + box.x1) / 2
                vertical = [d for d in drawings
                            if same_color(d['color'], (0.702, 0.702, 0.702))
                            and d['fill'] is None and d['rect'].width < 0.001
                            and abs(d['rect'].x0 - center) < 0.001
                            and axis.contains(d['rect'])]
                assert len(vertical) == 2, (panel, compound, gamma, vertical)
                upper_y = min(d['rect'].y0 for d in vertical)
                lower_y = max(d['rect'].y1 for d in vertical)
                upper_error = (box.y0 - upper_y) * scale
                lower_error = (lower_y - box.y0) * scale
                # PDF path coordinates have small rounding differences at joins.
                assert abs(upper_error - lower_error) < 0.002
                caps = [d for d in drawings
                        if same_color(d['color'], (0.702, 0.702, 0.702))
                        and d['fill'] is None and d['rect'].height < 0.001
                        and abs((d['rect'].x0 + d['rect'].x1) / 2 - center) < 0.001
                        and (abs(d['rect'].y0 - upper_y) < 0.001
                             or abs(d['rect'].y0 - lower_y) < 0.001)]
                assert len(caps) == 2, (panel, compound, gamma)
                beta = (axis.y1 - box.y0) * scale
                error = (upper_error + lower_error) / 2
                assert 0 <= beta <= maximum and error >= 0
                rows.append({
                    'PFAS': compound, 'NaCl_label_mM': salt,
                    'Gamma_label_umol_m2': gamma, 'main_figure_panel': panel,
                    'beta_vector_estimate_percent': beta,
                    'beta_error_bar_estimate_percentage_points': error,
                    'beta_display_percent': f'{beta:.1f}',
                    'beta_error_bar_display_percentage_points': f'{error:.1f}',
                    'status': 'figure_digitized_original_data_pending',
                    'source_pdf': PDF.name, 'source_pdf_sha256': source_hash,
                    'source_pdf_page': 1,
                })
                audit.append({'panel': panel, 'PFAS': compound, 'gamma_label': gamma,
                              'axis_rect_pt': list(axis), 'axis_y_max_percent': maximum,
                              'bar_drawing_index': drawing_index, 'bar_rect_pt': list(box),
                              'error_upper_y_pt': upper_y, 'error_lower_y_pt': lower_y,
                              'upper_lower_error_difference_percentage_points': abs(upper_error - lower_error)})
    assert len(rows) == 15
    with (DATA / 'md_binding_figure1_digitized.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    QA.mkdir(exist_ok=True)
    (QA / 'figure1_vector_extraction_audit.json').write_text(
        json.dumps({'source_pdf_sha256': source_hash,
                    'note': 'Axis ranges and legend identities checked against the original PDF rendering.',
                    'records': audit}, indent=2), encoding='utf-8')
    # Inspect the original panel geometry directly; no source figure is changed.
    clip = pymupdf.Rect(365, 0, page.rect.width, axes[1].y1 + 35)
    page.get_pixmap(matrix=pymupdf.Matrix(3, 3), clip=clip, alpha=False).save(
        QA / 'figure1_binding_panels_original.png')
    print(*[f"{r['main_figure_panel']} {r['PFAS']} Gamma={r['Gamma_label_umol_m2']} "
            f"beta={r['beta_display_percent']} +/- {r['beta_error_bar_display_percentage_points']}"
            for r in rows], sep='\n')


if __name__ == '__main__':
    main()
