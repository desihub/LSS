"""
build_fiber_scores.py — Build fiber_scores.json (5000 entries, one per DESI fiber).

Fields:
  fiber, petal, n_lrg, plog_score, plog_run_lengths, plog_windows,
  plog_sequence, lrg_observations, krolewski_nsig, krolewski_pval

Change from v1: MJD_TOL raised from 0.05 to 0.6 to capture early-survey
observations where the NPZ MJD differs from exposures-matterhorn.fits by
up to ~0.4 days (MEAN_MJD vs individual exposure MJD convention).
"""

import json, os
import numpy as np
from astropy.io import fits

# ── paths ─────────────────────────────────────────────────────────────────────
DATADIR  = os.path.dirname(os.path.abspath(__file__))
NPZ_PATH = os.path.join(DATADIR, 'LRG_fiber_time.npz')
EXP_FILE = '/global/cfs/cdirs/desi/spectro/redux/matterhorn/exposures-matterhorn.fits'
NSIG_FILE = ('/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/'
             'krolewski-pipeline/summaryscale/LRGfibersim_matterhorn-v2.txt')
OUT_PATH = os.path.join(DATADIR, 'fiber_scores.json')

PLOG_THRESH = 2.0   # minimum plog_score to populate detail fields
LOG_INV_P   = np.log10(1.0 / 0.0177)   # ~1.752
MJD_TOL     = 0.6   # days — raised from 0.05 to capture early-survey MJD offset

# ── load data ─────────────────────────────────────────────────────────────────
print(f'Loading {NPZ_PATH}...')
d = np.load(NPZ_PATH)
mjds_all = d['mjd']
zsuc_all = d['zsuc']
index    = d['index']

print(f'Loading {EXP_FILE}...')
exp_data  = fits.open(EXP_FILE)[1].data
exp_mjd   = np.array(exp_data['MJD'], dtype=float)
exp_tile  = np.array(exp_data['TILEID'], dtype=int)
exp_night = np.array(exp_data['NIGHT'], dtype=int)
sort_exp  = np.argsort(exp_mjd)
exp_mjd_s = exp_mjd[sort_exp]
exp_tile_s  = exp_tile[sort_exp]
exp_night_s = exp_night[sort_exp]

print(f'Loading Krolewski scores from {NSIG_FILE}...')
# Columns: FIBER, p_value, var, nsig  (nsig from simulation, not recomputed from p_value)
# p_value can be exactly 0 for very bad fibers; use simulation nsig directly to avoid ±inf
kdata = np.loadtxt(NSIG_FILE, usecols=(0, 1, 3))   # FIBER, p_value, nsig_sim
krow_by_fiber = {int(row[0]): (round(float(row[1]), 7), round(float(row[2]), 5))
                 for row in kdata}

# ── per-fiber computation ─────────────────────────────────────────────────────
def compute_plog(fib):
    """Returns (score, n_obs, windows, run_lengths, sequence) or (0,n,None,None,None)."""
    s, e = int(index[fib, 0]), int(index[fib, 1])
    n_obs = e - s
    if n_obs == 0:
        return 0.0, 0, None, None, None

    order = np.argsort(mjds_all[s:e])
    mj = mjds_all[s:e][order]
    zs = zsuc_all[s:e][order]

    score, windows, run_lengths = 0.0, [], []
    i = 0
    while i < len(zs):
        if not zs[i]:
            j = i
            while j < len(zs) and not zs[j]:
                j += 1
            k = j - i
            if k >= 2:
                contrib = k * LOG_INV_P - np.log10(n_obs)
                score += contrib
                windows.append({
                    'start': round(float(mj[i]), 2),
                    'end':   round(float(mj[j-1]), 2),
                    'k': k,
                    'contribution': round(contrib, 4)
                })
                run_lengths.append(k)
            i = j
        else:
            i += 1

    if score <= PLOG_THRESH:
        return score, n_obs, None, None, None

    # plog_sequence: '1'/'0' with ' ' for gaps > 3 days
    chars = ['1' if z else '0' for z in zs]
    seq_parts = [chars[0]]
    for idx in range(1, len(chars)):
        if float(mj[idx]) - float(mj[idx-1]) > 3.0:
            seq_parts.append(' ')
        seq_parts.append(chars[idx])
    seq = ''.join(seq_parts)

    return score, n_obs, windows, run_lengths, seq


def build_lrg_observations(fib):
    """Return list of {tileid, night, fail} for every LRG observation of this fiber."""
    s, e = int(index[fib, 0]), int(index[fib, 1])
    if e == s:
        return []

    order = np.argsort(mjds_all[s:e])
    mj = mjds_all[s:e][order]
    zs = zsuc_all[s:e][order]

    obs = []
    for m, z in zip(mj, zs):
        m = float(m)
        # Binary search for closest exposure by MJD
        idx_ins = np.searchsorted(exp_mjd_s, m)
        best_dist = MJD_TOL + 1
        best_tile, best_night = None, None
        for candidate in [idx_ins - 1, idx_ins, idx_ins + 1]:
            if 0 <= candidate < len(exp_mjd_s):
                dist = abs(exp_mjd_s[candidate] - m)
                if dist < best_dist:
                    best_dist = dist
                    best_tile  = int(exp_tile_s[candidate])
                    best_night = int(exp_night_s[candidate])
        if best_dist > MJD_TOL:
            best_tile  = None
            best_night = None
        obs.append({'tileid': best_tile, 'night': best_night, 'fail': bool(not z)})

    return obs


# ── main loop ─────────────────────────────────────────────────────────────────
print('Computing per-fiber scores...')
records = []
n_plog2 = 0
for fib in range(5000):
    petal = fib // 500
    score, n_lrg, windows, run_lengths, seq = compute_plog(fib)

    lrg_obs = None
    if score > PLOG_THRESH:
        n_plog2 += 1
        lrg_obs = build_lrg_observations(fib)

    kpval, knsig = krow_by_fiber.get(fib, (None, None))

    records.append({
        'fiber':            fib,
        'petal':            petal,
        'n_lrg':            n_lrg,
        'plog_score':       round(score, 5),
        'plog_run_lengths': run_lengths,
        'plog_windows':     windows,
        'plog_sequence':    seq,
        'lrg_observations': lrg_obs,
        'krolewski_nsig':   knsig,
        'krolewski_pval':   kpval,
    })

    if fib % 500 == 499:
        print(f'  petal {petal} done')

print(f'Fibers with plog > {PLOG_THRESH}: {n_plog2}')

# ── count unmatched ───────────────────────────────────────────────────────────
n_none = sum(
    1 for r in records if r['lrg_observations']
    for o in r['lrg_observations'] if o['tileid'] is None
)
n_total = sum(
    len(r['lrg_observations']) for r in records if r['lrg_observations']
)
print(f'lrg_observations unmatched (tileid=None): {n_none} / {n_total} '
      f'({100*n_none/max(n_total,1):.1f}%)')

# ── write output ──────────────────────────────────────────────────────────────
print(f'Writing {OUT_PATH}...')
with open(OUT_PATH, 'w') as f:
    json.dump(records, f, separators=(',', ':'))

mb = os.path.getsize(OUT_PATH) / 1e6
print(f'Done: {mb:.1f} MB, {len(records)} fibers')
