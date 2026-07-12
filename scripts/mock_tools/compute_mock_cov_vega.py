#!/usr/bin/env python
# coding: utf-8

"""
Build mock-ensemble mean and covariance from xirunpc xipoles outputs.

The primary ASCII output is a flat N×N covariance matrix in the same
ASCII layout used by Vega's ``build_global_cov.py`` (compatible with
RascalC-style files, but computed from mock scatter):

    [xi_0(s_1..s_n), xi_2(s_1..s_n), ..., xi_L(s_1..s_n)]

Examples::

    python compute_mock_cov_vega.py \\
        --xipoles-glob '$SCRATCH/QF3x2/mocks/mock-*/smu/xipoles_QSO_GCcomb_1.77_3.8_default_custom1_njack0_nran4_split20.txt' \\
        --outdir $SCRATCH/QF3x2/mocks/QSO_auto_cov \\
        --s-min 20 --s-max 200 \\
        --z-min 1.77 --z-max 3.8

    python $HOME/Repos/vega/bin/build_global_cov.py \\
        --lya-cov /path/to/lya-qso-cross-covariance.fits \\
        --qso-cov $SCRATCH/QF3x2/mocks/QSO_auto_cov/xi024_QSO_GCcomb_z1.77-3.8_mock_s20-200_cov.txt \\
        --output /path/to/combined_global_cov.fits
"""

import argparse
import glob
import os
import re

import numpy as np


def load_xipoles(path, ells=(0, 2, 4)):
    """
    Read xipoles_*.txt from xirunpc.

    Columns: s_mid, s_avg, xi_0, xi_2, [xi_4, ...]
    """
    raw = np.loadtxt(path, comments='#')
    s_mid = raw[:, 0]
    s_avg = raw[:, 1]
    xi = np.column_stack([raw[:, 2 + ell // 2] for ell in ells])
    return s_mid, s_avg, xi


def apply_s_cuts(s_mid, s_avg, xi, s_min, s_max):
    """Keep bins with s_mid in [s_min, s_max)."""
    mask = (s_mid >= s_min) & (s_mid < s_max)
    if not np.any(mask):
        raise ValueError(
            f'no s-bins in [{s_min}, {s_max}); available range '
            f'{s_mid.min():.1f} – {s_mid.max():.1f}')
    return s_mid[mask], s_avg[mask], xi[mask]


def xi_to_data_vector(xi):
    """Flatten (n_s, n_ells) array to RascalC / Vega data-vector order."""
    return np.concatenate([xi[:, i] for i in range(xi.shape[1])])


def save_ascii_cov(cov, path, comment=''):
    """Write a square covariance matrix as a flat ASCII table (one row per line)."""
    cov = np.asarray(cov, dtype=np.float64)
    if cov.ndim != 2 or cov.shape[0] != cov.shape[1]:
        raise ValueError(f'covariance must be square, got {cov.shape}')
    with open(path, 'w') as f:
        if comment:
            f.write(f'# {comment}\n')
        for row in cov:
            f.write(' '.join(f'{x:.18e}' for x in row) + '\n')


def save_vega_multipoles(path, s_mid, s_avg, xi_mean, cov, ells):
    """Write mock mean + mock-cov std columns in Vega multipole ASCII format."""
    n_ells = len(ells)
    n_s = len(s_mid)
    std = np.sqrt(np.diag(cov)).reshape(n_ells, n_s)
    cols = (['s_mid', 's_avg']
            + [f'xi_{ell}' for ell in ells]
            + [f'std_{ell}' for ell in ells])
    table = np.column_stack(
        [s_mid, s_avg]
        + [xi_mean[:, i] for i in range(n_ells)]
        + [std[i] for i in range(n_ells)])
    header = ' '.join(cols)
    np.savetxt(path, table, header=header, comments='# ')


def default_cov_name(outdir, z_min, z_max, s_min, s_max, ells):
    """Default output filename for the mock covariance matrix."""
    ell_label = ''.join(str(ell) for ell in ells)
    ztag = f'z{z_min:g}-{z_max:g}'
    return os.path.join(
        outdir,
        f'xi{ell_label}_QSO_GCcomb_{ztag}_mock_s{int(s_min)}-{int(s_max)}_cov.txt')


def default_mean_name(cov_path):
    """Default output filename for the mock mean multipole table."""
    if cov_path.endswith('_cov.txt'):
        return cov_path[:-8] + '_mean.txt'
    root, _ = os.path.splitext(cov_path)
    return root + '_mean.txt'


def parse_z_from_xipoles_path(path):
    """Try to read z limits from an xipoles filename."""
    match = re.search(r'_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)_', os.path.basename(path))
    if match:
        return float(match.group(1)), float(match.group(2))
    return None, None


def _natural_sort_key(path):
    """Sort mock-2 before mock-10 when paths contain numeric tokens."""
    parts = re.split(r'(\d+)', path)
    key = []
    for part in parts:
        if part.isdigit():
            key.append((0, int(part)))
        else:
            key.append((1, part))
    return key


def expand_paths(patterns):
    """Expand shell-style globs; keep literal paths that exist."""
    paths = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern), key=_natural_sort_key)
        if matches:
            paths.extend(matches)
        elif os.path.exists(pattern):
            paths.append(pattern)
        else:
            raise FileNotFoundError(f'no matches for pattern: {pattern}')
    seen = set()
    unique = []
    for path in paths:
        if path not in seen:
            seen.add(path)
            unique.append(path)
    return unique


def resolve_xipoles_files(mock_dirs=None, xi_subdir='smu', xipoles_pattern=None,
                           xipoles_glob=None):
    """Return sorted list of xipoles file paths."""
    if xipoles_glob is not None:
        files = expand_paths([xipoles_glob])
        files = [f for f in files if os.path.isfile(f)]
        if not files:
            raise FileNotFoundError(f'no xipoles files matched: {xipoles_glob}')
        return files

    if not mock_dirs:
        raise ValueError('provide --mock-dirs and/or --xipoles-glob')

    mock_dirs = expand_paths(mock_dirs)
    files = []
    for mock_dir in mock_dirs:
        if not os.path.isdir(mock_dir):
            raise NotADirectoryError(mock_dir)
        xi_dir = os.path.join(mock_dir, xi_subdir)
        if xipoles_pattern is not None:
            path = os.path.join(xi_dir, xipoles_pattern)
            if not os.path.isfile(path):
                raise FileNotFoundError(path)
            files.append(path)
        else:
            matches = sorted(glob.glob(os.path.join(xi_dir, 'xipoles_QSO_GCcomb_*.txt')))
            if len(matches) != 1:
                raise RuntimeError(
                    f'expected 1 xipoles file in {xi_dir}, found {len(matches)}')
            files.append(matches[0])
    return files


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '--mock-dirs', nargs='+', default=None,
        help='mock output directories; globs allowed, e.g. .../mock-*')
    parser.add_argument(
        '--xipoles-glob', default=None,
        help='glob of xipoles files directly, e.g. .../mock-*/smu/xipoles_....txt')
    parser.add_argument('--xi-subdir', default='smu',
                        help='subdirectory under each mock dir (default: smu)')
    parser.add_argument(
        '--xipoles-pattern', default=None,
        help='filename inside xi-subdir when using --mock-dirs')
    parser.add_argument('--ells', type=int, nargs='+', default=[0, 2, 4])
    parser.add_argument('--s-min', type=float, default=20.0,
                        help='lower s cut (Mpc/h, inclusive; matches qsoxqso.ini)')
    parser.add_argument('--s-max', type=float, default=200.0,
                        help='upper s cut (Mpc/h, exclusive; matches qsoxqso.ini)')
    parser.add_argument('--z-min', type=float, default=None,
                        help='redshift lower limit for output filename (parsed from '
                             'xipoles filename when omitted)')
    parser.add_argument('--z-max', type=float, default=None,
                        help='redshift upper limit for output filename')
    parser.add_argument('--outdir', required=True,
                        help='where to write mean/cov/output files')
    parser.add_argument(
        '--cov-out', dest='cov_out', default=None,
        help='ASCII covariance output path (default: auto name in --outdir)')
    parser.add_argument(
        '--rascalc-out', dest='cov_out', default=None,
        help=argparse.SUPPRESS)
    parser.add_argument(
        '--mean-out', default=None,
        help='optional Vega multipole ASCII file for the mock mean '
             '(default: replace _cov.txt with _mean.txt in --cov-out)')
    args = parser.parse_args()

    if args.mock_dirs is None and args.xipoles_glob is None:
        parser.error('provide --mock-dirs and/or --xipoles-glob')

    os.makedirs(args.outdir, exist_ok=True)
    ells = tuple(args.ells)
    n_ells = len(ells)

    used_files = resolve_xipoles_files(
        mock_dirs=args.mock_dirs,
        xi_subdir=args.xi_subdir,
        xipoles_pattern=args.xipoles_pattern,
        xipoles_glob=args.xipoles_glob,
    )

    z_min, z_max = args.z_min, args.z_max
    if z_min is None or z_max is None:
        pz_min, pz_max = parse_z_from_xipoles_path(used_files[0])
        z_min = z_min if z_min is not None else pz_min
        z_max = z_max if z_max is not None else pz_max

    s_mid_ref = None
    s_avg_sum = None
    vectors = []
    for fn in used_files:
        s_mid, s_avg, xi = load_xipoles(fn, ells=ells)
        s_mid, s_avg, xi = apply_s_cuts(s_mid, s_avg, xi, args.s_min, args.s_max)
        if s_mid_ref is None:
            s_mid_ref = s_mid
            s_avg_sum = np.zeros_like(s_avg)
        elif len(s_mid) != len(s_mid_ref) or not np.allclose(s_mid, s_mid_ref):
            raise ValueError(
                f's-bin mismatch in {fn}: expected the same s_mid grid as the first mock')
        s_avg_sum += s_avg
        vectors.append(xi_to_data_vector(xi))
        print(f'loaded {fn}  (N={len(vectors[-1])})')
    s_avg_ref = s_avg_sum / len(used_files)

    X = np.array(vectors)
    nmock, ndata = X.shape
    xi_mean = X.mean(axis=0)
    cov = np.cov(X, rowvar=False, ddof=1)

    n_s = ndata // n_ells
    if n_s * n_ells != ndata:
        raise ValueError(f'data vector length {ndata} is not divisible by n_ells={n_ells}')
    xi_mean_2d = np.column_stack(
        [xi_mean[i * n_s:(i + 1) * n_s] for i in range(n_ells)])

    ds = (args.s_max - args.s_min) / n_s
    s_cov_centers = args.s_min + (np.arange(n_s) + 0.5) * ds
    max_s_mismatch = np.max(np.abs(s_cov_centers - s_avg_ref))
    if max_s_mismatch > ds:
        print(f'WARNING: inferred covariance s-grid differs from xipoles s_avg by up to '
              f'{max_s_mismatch:.2f} Mpc/h (bin width {ds:.2f}). '
              f'Check --s-min/--s-max against the xipoles binning.')

    cov_out = args.cov_out
    if cov_out is None:
        if z_min is None or z_max is None:
            cov_out = os.path.join(
                args.outdir,
                f'xi{"".join(str(e) for e in ells)}_QSO_GCcomb_mock_s{int(args.s_min)}-{int(args.s_max)}_cov.txt')
        else:
            cov_out = default_cov_name(
                args.outdir, z_min, z_max, args.s_min, args.s_max, ells)

    comment = (f'mock-ensemble covariance from {nmock} realizations; '
               f'ells={ells}; s in [{args.s_min}, {args.s_max}); '
               f'order=[xi_0(s bins), xi_2(s bins), ...]')
    save_ascii_cov(cov, cov_out, comment=comment)

    mean_out = args.mean_out
    if mean_out is None:
        mean_out = default_mean_name(cov_out)
    save_vega_multipoles(mean_out, s_mid_ref, s_avg_ref, xi_mean_2d, cov, ells)

    np.save(os.path.join(args.outdir, 'xi_mean.npy'), xi_mean)
    np.save(os.path.join(args.outdir, 'cov.npy'), cov)
    np.save(os.path.join(args.outdir, 's_mid.npy'), s_mid_ref)
    np.save(os.path.join(args.outdir, 's_avg.npy'), s_avg_ref)
    np.savetxt(os.path.join(args.outdir, 'mock_list.txt'), used_files, fmt='%s')

    print(f'\nNmock = {nmock}, Ndata = {ndata}  ({n_s} s-bins × {n_ells} ells)')
    print(f' wrote {cov_out}')
    print(f' wrote {mean_out}')
    print(f' wrote {args.outdir}/cov.npy')

    eigvals = np.linalg.eigvalsh(cov)
    min_eig = eigvals.min()
    if min_eig <= 0:
        print(f' WARNING: covariance has non-positive eigenvalue {min_eig:.3e}')
    else:
        print(f' covariance is positive definite (min eigenvalue = {min_eig:.3e})')

    try:
        icov = np.linalg.inv(cov)
        chi2 = [float((x - xi_mean) @ icov @ (x - xi_mean)) for x in X]
        print(f' chi2/mock: min={min(chi2):.1f}, mean={np.mean(chi2):.1f}, max={max(chi2):.1f}')
        print(f' (if cov is good, mean chi2 ~ {ndata}; with few mocks this will not hold yet)')
    except np.linalg.LinAlgError:
        print(' cov is singular — need more mocks and/or fewer bins/multipoles')


if __name__ == '__main__':
    main()
