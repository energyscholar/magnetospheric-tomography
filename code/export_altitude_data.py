#!/usr/bin/env python3
"""Export altitude-resolved Moran's I data as structured JSON for Paper 1.

Uses the existing relational_diagnostics pipeline (two-layer Chapman model)
to compute I(ground), I(120km), I(300km) per epoch, then exports structured
JSON for figure generation and statistical analysis.

Usage:
    python3 export_altitude_data.py --event "Halloween 2003"
    python3 export_altitude_data.py --all
    python3 export_altitude_data.py --all --resolution 30min
"""

import sys
import json
import argparse
import numpy as np
from pathlib import Path
from scipy import stats as sp_stats

import loader_supermag as sm
import loader_tec as tec
import relational_diagnostics as rd
import event_catalog as catalog

BASE_DIR = Path(__file__).parent.parent
RESULTS_DIR = BASE_DIR / 'results' / 'paper1_data'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def safe_float(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return None
    return float(v)


def export_event(event, resolution='2hr'):
    print(f"\n{'='*60}")
    print(f"Exporting: {event.name} ({event.date_str})")
    print(f"{'='*60}")

    sm_data, sm_coords = sm.load_event_day(
        event.year, event.month, event.day, min_valid_frac=0.3)
    if not sm_data:
        print("  SKIP: No SuperMAG data")
        return None

    tec_maps = tec.load_day(event.year, event.month, event.day)
    print(f"  SuperMAG: {len(sm_data)} stations")
    print(f"  TEC: {len(tec_maps) if tec_maps else 0} epochs")

    table = rd.multi_epoch_table(sm_data, sm_coords, tec_maps, resolution)

    epochs = []
    for ut_hours, m in table:
        if m is None:
            continue

        alt_I = m.get('altitude_I', {})
        record = {
            'ut_hours': ut_hours,
            'I_N': safe_float(m.get('I_N')),
            'I_E': safe_float(m.get('I_E')),
            'I_Z': safe_float(m.get('I_Z')),
            'I_120': safe_float(alt_I.get(120)),
            'I_300': safe_float(alt_I.get(300)),
            'delta_120_300': safe_float(m.get('delta_120_300')),
            'rms_nT': safe_float(m.get('rms_N')),
            'n_stations': m.get('n_stations'),
            'n_tec_stations': m.get('n_tec_stations'),
            'coupling_pct': safe_float(m.get('coupling_pct')),
            'vario_ratio': safe_float(m.get('vario_ratio')),
            'f_auroral': safe_float(m.get('f_auroral')),
        }
        if alt_I:
            record['altitude_profile'] = {
                str(k): safe_float(v) for k, v in sorted(alt_I.items())
            }
        epochs.append(record)

    if not epochs:
        print("  SKIP: No valid epochs")
        return None

    event_data = {
        'event_name': event.name,
        'date': event.date_str,
        'category': event.category,
        'dst_min': event.dst_min,
        'resolution': resolution,
        'n_epochs': len(epochs),
        'epochs': epochs,
    }

    slug = event.slug
    data_path = RESULTS_DIR / f'{slug}_altitude_data.json'
    with open(data_path, 'w') as f:
        json.dump(event_data, f, indent=2)
    print(f"  Data: {data_path} ({data_path.stat().st_size / 1024:.1f} KB)")

    summary = compute_summary(event_data)
    summary_path = RESULTS_DIR / f'{slug}_summary.json'
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  Summary: {summary_path}")

    rd.print_summary(table)
    return event_data


def compute_summary(event_data):
    epochs = event_data['epochs']

    def extract(key):
        return np.array([e[key] for e in epochs if e.get(key) is not None])

    I_N = extract('I_N')
    I_120 = extract('I_120')
    I_300 = extract('I_300')
    delta = extract('delta_120_300')
    rms = extract('rms_nT')

    summary = {
        'event_name': event_data['event_name'],
        'date': event_data['date'],
        'n_epochs': len(epochs),
        'n_with_altitude': int(len(I_120)),
    }

    for name, arr in [('I_N', I_N), ('I_120', I_120), ('I_300', I_300),
                      ('delta_120_300', delta), ('rms_nT', rms)]:
        if len(arr) >= 2:
            summary[f'{name}_mean'] = safe_float(np.mean(arr))
            summary[f'{name}_std'] = safe_float(np.std(arr))
            summary[f'{name}_min'] = safe_float(np.min(arr))
            summary[f'{name}_max'] = safe_float(np.max(arr))

    pairs = [('I_N', 'I_120', I_N, I_120),
             ('I_N', 'I_300', I_N, I_300),
             ('I_120', 'I_300', I_120, I_300)]

    correlations = {}
    for name_a, name_b, a, b in pairs:
        n = min(len(a), len(b))
        if n >= 4:
            a_use, b_use = a[:n], b[:n]
            mask = np.isfinite(a_use) & np.isfinite(b_use)
            if np.sum(mask) >= 4:
                r, p = sp_stats.pearsonr(a_use[mask], b_use[mask])
                correlations[f'{name_a}_vs_{name_b}'] = {
                    'r': safe_float(r),
                    'R2': safe_float(r ** 2),
                    'p': safe_float(p),
                    'n': int(np.sum(mask)),
                }
    summary['correlations'] = correlations

    return summary


def main():
    parser = argparse.ArgumentParser(
        description='Export altitude-resolved Moran\'s I data for Paper 1')
    parser.add_argument('--all', action='store_true', help='Process all events')
    parser.add_argument('--event', type=str, help='Event name (partial match)')
    parser.add_argument('--resolution', default='2hr',
                        choices=['2hr', '30min'], help='Temporal resolution')
    args = parser.parse_args()

    if args.all:
        catalog.verify_availability()
        events = [e for e in catalog.EVENTS if e.has_supermag]
        print(f"Exporting {len(events)} events at {args.resolution} resolution")

        all_summaries = []
        for event in events:
            data = export_event(event, args.resolution)
            if data:
                summary = compute_summary(data)
                all_summaries.append(summary)

        cross_path = RESULTS_DIR / 'cross_event_summary.json'
        with open(cross_path, 'w') as f:
            json.dump(all_summaries, f, indent=2)
        print(f"\nCross-event summary: {cross_path}")

        print(f"\n{'='*60}")
        print("CROSS-EVENT SUMMARY")
        print(f"{'='*60}")
        print(f"\n{'Event':<25} {'I_N':>6} {'I_120':>6} {'I_300':>6} "
              f"{'Delta':>6} {'r(N,120)':>8} {'r(N,300)':>8} {'r(120,300)':>9}")
        print('-' * 85)
        for s in all_summaries:
            def g(k):
                v = s.get(k)
                return f'{v:+.3f}' if v is not None else '  N/A'
            def gc(k):
                c = s.get('correlations', {}).get(k, {})
                v = c.get('r')
                return f'{v:+.3f}' if v is not None else '  N/A'
            print(f"{s['event_name']:<25} {g('I_N_mean'):>6} {g('I_120_mean'):>6} "
                  f"{g('I_300_mean'):>6} {g('delta_120_300_mean'):>6} "
                  f"{gc('I_N_vs_I_120'):>8} {gc('I_N_vs_I_300'):>8} "
                  f"{gc('I_120_vs_I_300'):>9}")

    elif args.event:
        matches = [e for e in catalog.EVENTS
                   if args.event.lower() in e.name.lower()]
        if not matches:
            print(f"No event matching '{args.event}'")
            sys.exit(1)
        export_event(matches[0], args.resolution)

    else:
        export_event(catalog.EVENTS[0], args.resolution)


if __name__ == '__main__':
    main()
