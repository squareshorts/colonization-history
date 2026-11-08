"""Observed time to sustained occupancy increase in archived tandem flow-cell images.

Inputs are area_px counts previously extracted from Lee et al. (2018) binary
image stacks, supplied under empirical/observables/ in the public repository.
The endpoint is a within-record fold increase over early observed occupancy;
it is NOT an exponential intercept or an estimate of individual-cell memory.

Recompute:
  python observed_establishment.py --repo /path/to/colonization-history \
      --out /path/to/observed_establishment.csv \
      --sensitivity /path/to/robustness.csv
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

DATASETS = [
    ('FC1-1', 'WT YFP-CFP (T=29.5h;FC1)'),
    ('FC2-1', 'WT YFP-CFP (T=29.5h;FC2)'),
    ('FC1-2', 'WT YFP-CFP (T=26.75h;FC1)'),
    ('FC2-2', 'WT YFP-CFP (T=26.75h;FC2)'),
    ('FC1-3', 'WT YFP-CFP (T1=32.5h,T2=14.75h;FC1)'),
    ('FC2-3', 'WT YFP-CFP (T1=32.5h,T2=14.75h;FC2)'),
    ('FC3-3', 'WT YFP-CFP (T1=32.5h,T2=14.75h;FC3)'),
    ('Tr14', 'WT YFP-CFP (Tr=14h)'),
    ('Tr37', 'WT YFP-CFP (Tr=37h)'),
]


def moving_average(values, width=20):
    """Centered moving average normalized at both record boundaries."""
    values = np.asarray(values, dtype=float)
    if width <= 0:
        raise ValueError('width must be positive')
    k = np.ones(int(width), dtype=float)
    return np.convolve(values, k, mode='same') / np.convolve(np.ones_like(values), k, mode='same')


def time_to_sustained_fold(times, occupancy, *, fold=5., width=20,
                           persistence_minutes=15., minimum_baseline=1.,
                           baseline_hours=2.):
    """First observed time of a maintained fold increase, or NaN if censored.

    Each record's early baseline is the median of smoothed observed occupancy
    for t < 2 h, with the stated optional one-cell-equivalent floor. A
    crossing must begin at/after that reference period and persist through
    the specified real-time horizon. No extrapolation is performed.
    """
    t, y = np.asarray(times, float), np.asarray(occupancy, float)
    if len(t) != len(y) or len(t) < 4 or np.any(np.diff(t) <= 0):
        raise ValueError('times must be increasing and match observations')
    if not np.isfinite(y).all() or np.min(y) < 0 or fold <= 1:
        raise ValueError('occupancy must be finite/nonnegative and fold > 1')
    smoothed = moving_average(y, width)
    early = smoothed[t < baseline_hours]
    if len(early) == 0:
        raise ValueError('no early baseline samples')
    baseline = max(float(np.median(early)), float(minimum_baseline))
    target = fold * baseline
    persist_hours = float(persistence_minutes)/60.
    for i in np.flatnonzero((t >= baseline_hours) & (smoothed >= target)):
        if t[-1] < t[i] + persist_hours:
            continue
        end = int(np.searchsorted(t, t[i] + persist_hours, side='left'))
        if np.all(smoothed[i:end+1] >= target):
            return float(t[i]), baseline, target
    return float('nan'), baseline, target


def read_record(repo, name):
    loc = Path(repo) / 'empirical' / 'observables' / name
    with (loc/'counts.csv').open(newline='') as f:
        rows = list(csv.DictReader(f))
    meta = json.loads((loc/'meta.json').read_text())
    time = np.array([float(r['time_h']) for r in rows])
    area = np.array([float(r['area_px']) for r in rows])
    return time, area, float(meta['median_component_area_first2h'])


def paired_summary(records):
    by = {r['record']: r for r in records}
    vals = np.array([
        float(by[f'FC2-{j}']['time_h']) - float(by[f'FC1-{j}']['time_h'])
        for j in (1, 2, 3)], float)
    if not np.all(np.isfinite(vals)):
        raise ValueError('paired endpoint is censored in one or more records')
    mean = float(np.mean(vals)); se = float(np.std(vals, ddof=1)/np.sqrt(3))
    # Descriptive paired Student-t CI, df=2 (three matched experimental series).
    t975_df2 = 4.302652729911275
    return {'paired_differences_h': vals.tolist(), 'paired_mean_h': mean,
            'descriptive_95pct_interval_h': [mean-t975_df2*se, mean+t975_df2*se],
            'n_independent_series': 3}


def calculate(repo, fold=5., width=20, persistence_minutes=15.,
              minimum_baseline=1.):
    controls = [read_record(repo, name)[2] for id, name in DATASETS if id.startswith('FC1')]
    shared = float(np.median(controls))
    records = []
    for id, name in DATASETS:
        t, area, local = read_record(repo, name)
        val, base, target = time_to_sustained_fold(
            t, area/local, fold=fold, width=width,
            persistence_minutes=persistence_minutes,
            minimum_baseline=minimum_baseline)
        common, _, _ = time_to_sustained_fold(
            t, area/shared, fold=fold, width=width,
            persistence_minutes=persistence_minutes,
            minimum_baseline=minimum_baseline)
        records.append({'record': id, 'source_record': name, 'time_h': val,
                        'shared_calibration_time_h': common,
                        'early_baseline_cell_equivalents': base,
                        'threshold_cell_equivalents': target,
                        'record_end_h': float(t[-1]),
                        'local_cell_area_px': local, 'shared_cell_area_px': shared,
                        'censored': not np.isfinite(val)})
    return records, paired_summary(records)


def write_csv(path, records):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0])); w.writeheader(); w.writerows(records)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', type=Path, required=True)
    ap.add_argument('--out', type=Path, default=Path('observed_establishment.csv'))
    ap.add_argument('--sensitivity', type=Path, default=Path('observed_establishment_sensitivity.csv'))
    ap.add_argument('--full-sensitivity', type=Path, default=None)
    args = ap.parse_args()
    primary, summ = calculate(args.repo)
    no_floor,_ = calculate(args.repo, minimum_baseline=0.)
    records=[]
    for r,raw in zip(primary,no_floor):
        records.append({
            'record':r['record'],'source_record':r['source_record'],
            'time_h':round(r['time_h'],7),
            'shared_calibration_time_h':round(r['shared_calibration_time_h'],7),
            'no_floor_time_h':round(raw['time_h'],7),
            'last_record_h':round(r['record_end_h'],7)})
    write_csv(args.out, records)
    rows=[]; full=[]
    for fold in (3., 5.):
        for width in (10, 20, 40):
            for dwell in (0., 15., 30.):
                _,stats=calculate(args.repo,fold=fold,width=width,persistence_minutes=dwell)
                item={'fold':int(fold),'width':width,'persistence_min':int(dwell),
                      **{f'pair_{i}_h':round(v,7) for i,v in enumerate(stats['paired_differences_h'],1)}}
                full.append(item)
                if dwell==15.:rows.append(item)
    write_csv(args.sensitivity, rows)
    if args.full_sensitivity is not None:write_csv(args.full_sensitivity,full)
    print(json.dumps({'primary':summ,'robustness_settings':len(full)},indent=2))


if __name__ == '__main__':
    main()
