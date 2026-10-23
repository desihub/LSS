# compare_dr3_matterhorn_stats — DR3 daily vs matterhorn-v2 column statistics

*James Rohlf and Claude Sonnet 4.6 — 2026-05-21*

---

## Purpose

Compare per-column statistics (mean, sigma, [min, max] range) between the
DR3 daily catalog and the matterhorn-v2 reprocessing, for all available tracers:
LRG, ELG_LOPnotqso, ELGnotqso, QSO, BGS_BRIGHT, and BGS_FAINT.  The script
attempts all tracers and skips gracefully if a file is missing.  Both catalogs
use the `full.dat.fits` (pre-clustering-cut) files from their respective
`test/` directories.  Statistics are computed after applying the standard
z-cut for each tracer.

## Files

| File | Description |
|------|-------------|
| `compare_dr3_matterhorn_stats.py` | Analysis script |
| `compare_dr3_matterhorn_stats.txt` | Saved output table |

## Data paths

| Catalog | Path |
|---------|------|
| DR3 daily | `/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/daily/LSScats/test` |
| matterhorn-v2 | `/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/matterhorn-v2/LSScats/test` |

## Tracers and z-cuts

| Tracer | z range | Status |
|--------|---------|--------|
| LRG | 0.4 – 1.1 | both catalogs |
| ELG_LOPnotqso | 0.8 – 1.6 | both catalogs |
| ELGnotqso | 0.8 – 1.6 | both catalogs |
| QSO | 0.8 – 2.1 | both catalogs |
| BGS_BRIGHT | 0.1 – 0.4 | both catalogs |
| BGS_FAINT | 0.1 – 0.4 | both catalogs |
| BGS_ANY | 0.1 – 0.4 | file missing in both |
| ELG_LOP | 0.8 – 1.6 | file missing in both |

## Columns compared

`Z_not4clus`, `WEIGHT_ZFAIL`, `FRAC_TLOBS_TILES`, `COMP_TILE`,
`FRACZ_TILELOCID`, `mod_success_rate`, `TSNR2_LRG`, `TSNR2_ELG`, `TSNR2_QSO`

## Key findings

**Object counts:** matterhorn-v2 has ~0.2% fewer rows than daily across all
tracers, and 1.6–2.2% fewer after the z-cut.

**TSNR2 outliers:** The daily catalog has extreme TSNR2 outliers in LRG
(sigma = 522, max = 1.3×10⁶), ELG_LOPnotqso, and QSO that are absent in
matterhorn-v2 (sigma ≈ 13, max ≈ 250). BGS_BRIGHT and BGS_FAINT TSNR2 values
are clean in both catalogs. matterhorn-v2 has capped or removed outliers for
the dark-program tracers.

**Missing columns:** `WEIGHT_ZFAIL` and `mod_success_rate` are missing from
the daily catalog for ELG_LOPnotqso, ELGnotqso, BGS_BRIGHT, and BGS_FAINT.
Both columns are absent in both catalogs for QSO.

**Completeness columns** (`COMP_TILE`, `FRACZ_TILELOCID`): consistently lower
means in matterhorn-v2 across all tracers, consistent with stricter quality
cuts in that processing.

**Redshift distributions** (`Z_not4clus`): means and sigmas agree to four
decimal places across all tracers — the two catalogs select the same population.

**BGS tracers:** BGS_BRIGHT and BGS_FAINT are both present and consistent.
BGS_ANY is absent from both catalogs.
