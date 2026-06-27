"""
find_long_runs.py
Find LRG fibers with consecutive failure runs of length >= MIN_RUN.

Usage:
  python find_long_runs.py              # default MIN_RUN=3
  python find_long_runs.py --min 4
"""

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--min', type=int, default=3, dest='min_run',
                    help='Minimum consecutive failure run length (default: 3)')
parser.add_argument('--npz', default='/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/'
                    'all-fibers-vs-time/LRG_fiber_time.npz',
                    help='Path to LRG_fiber_time.npz')
args = parser.parse_args()

print(f'Loading {args.npz} ...')
d = np.load(args.npz)
mjd_all  = d['mjd']
zsuc_all = d['zsuc']
idx      = d['index']

MIN_RUN = args.min_run
print(f'Scanning 5000 fibers for runs of length >= {MIN_RUN} ...\n')

results = []

for fib in range(5000):
    s, e = idx[fib]
    if s < 0:
        continue
    mjds = mjd_all[s:e]
    zsuc = zsuc_all[s:e]
    if len(mjds) < MIN_RUN:
        continue

    # Sort by MJD
    order = np.argsort(mjds)
    zsuc_sorted = zsuc[order]
    mjds_sorted = mjds[order]

    # Find consecutive failure runs
    max_run = 0
    runs = []
    i, n = 0, len(zsuc_sorted)
    while i < n:
        if not zsuc_sorted[i]:
            j = i
            while j < n and not zsuc_sorted[j]:
                j += 1
            run_len = j - i
            if run_len >= MIN_RUN:
                runs.append((run_len, float(mjds_sorted[i]), float(mjds_sorted[j-1])))
            max_run = max(max_run, run_len)
            i = j
        else:
            i += 1

    if runs:
        petal = fib // 500
        results.append((fib, petal, max_run, len(mjds), runs))

# Sort by max run length descending
results.sort(key=lambda x: -x[2])

from datetime import datetime, timedelta
def mjd2date(m):
    return str((datetime(1858, 11, 17) + timedelta(days=float(m))).date())

print(f'Found {len(results)} fibers with at least one run >= {MIN_RUN}\n')
print(f'{"fiber":>6}  {"petal":>5}  {"n_obs":>6}  {"max_run":>8}  runs (len, start, end)')
print('-' * 75)
for fib, petal, max_run, n_obs, runs in results:
    run_strs = [f'len={r[0]} ({mjd2date(r[1])}–{mjd2date(r[2])})' for r in runs]
    print(f'{fib:>6}  {petal:>5}  {n_obs:>6}  {max_run:>8}  {", ".join(run_strs)}')

print(f'\nTotal: {len(results)} fibers')
