# all-fibers-vs-time

*James W. Rohlf (rohlf@bu.edu), June 2026*

Per-fiber failure-rate vs time plots for DESI DA3 / matterhorn-v2.
Supports all four tracers: BGS, LRG, ELG, QSO.

---

## Quick start

```bash
python plot_fiber_time.py <fiber> <TRACER> [<TRACER> ...]
```

Examples:
```bash
python plot_fiber_time.py 725 BGS LRG ELG QSO
python plot_fiber_time.py 986 LRG ELG
python plot_fiber_time.py 2014 LRG
```

Output PNG is written to this directory, named automatically:
`fiber<N>_time_<tracers>.png`

---

## Files

### Data files (pre-built, load instantly)

| File | Tracer | Rows | Size |
|------|--------|------|------|
| `LRG_fiber_time.npz` | LRG | 8.0M | 19 MB |
| `ELG_fiber_time.npz` | ELG_LOPnotqso | 17.5M | 44 MB |
| `BGS_fiber_time.npz` | BGS_BRIGHT | 10.7M | 22 MB |
| `QSO_fiber_time.npz` | QSO | 4.6M | 10 MB |

Each NPZ contains:
- `fiber`  — int16 (N,): fiber number (0–4999)
- `mjd`    — float32 (N,): mean tile MJD
- `zsuc`   — bool (N,): True = good redshift
- `index`  — int32 (5000, 2): `[start, end)` row slice for each fiber; `[-1,-1]` if absent

Arrays are sorted by fiber. To load one fiber:
```python
d = np.load('LRG_fiber_time.npz')
s, e = d['index'][725]
mjd  = d['mjd'][s:e]
zsuc = d['zsuc'][s:e]
```

### Scripts

| Script | Purpose |
|--------|---------|
| `make_fiber_time_data.py` | Builds the four NPZ files from the raw LSS catalogs (~85s) |
| `plot_fiber_time.py` | Makes failure-rate vs time plot for one fiber (<1s) |
| `plot_fiber725_time_test.py` | Earlier development/test script (superseded) |

---

## Plot description

Each panel shows **failure rate (1 − z_suc / z_tot) vs mean tile MJD**, binned
into 32 equal-width bins over the [2nd, 98th] percentile of the fiber's MJD
range.  Error bars are binomial.  The dashed gray line is the overall mean
failure rate.

**Top axis:** tick marks only (no labels).  Long ticks = January 1 of each
year; short ticks = first of each month.

**Vertical dashed lines:** CCD swap dates for the petal containing the fiber
(petal = fiber // 500), color-coded by spectrograph channel:

| Color | Channel |
|-------|---------|
| Blue | b (blue CCD) |
| Red | r (red CCD) |
| Purple | z (NIR CCD) |

CCD swap dates from `desispec/etc/get_ccd_history.py` as of 2025-05-12.

---

## Selection cuts

| Tracer | z_tot | z_suc |
|--------|-------|-------|
| LRG | TSNR2_ELG > 80, FIBERSTATUS ∈ {0,8} | ZWARN=0, Δχ²>15, Z<1.5 |
| ELG | TSNR2_ELG > 80, FIBERSTATUS ∈ {0,8} | ZWARN=0, Δχ²>9, 0.6<Z<1.6 |
| BGS | TSNR2_BGS > 1000, FIBERSTATUS ∈ {0,8} | ZWARN=0, Δχ²>25, 0.001<Z<0.6 |
| QSO | TSNR2_ELG > 80, FIBERSTATUS ∈ {0,8} | ZWARN=0, Δχ²>25, 0.8<Z<3.5 |

---

## Source catalogs

All from `/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/matterhorn-v2/LSScats/test/`:

| Tracer | Catalog file |
|--------|-------------|
| LRG | `LRG_full_noveto.dat.fits` |
| ELG | `ELG_LOPnotqso_full_noveto.dat.fits` |
| BGS | `BGS_BRIGHT_full_noveto.dat.fits` |
| QSO | `QSO_full_noveto.dat.fits` |

Exposure timestamps from:
`/global/cfs/cdirs/desi/spectro/redux/matterhorn/exposures-matterhorn.fits`

LRG nsig values (MC significance of fiber badness) from:
`krolewski-pipeline/summaryscale/LRGfibersim_matterhorn-v2.txt`

---

## Interactive Focal Plane Viewer (LRG)

*J. Rohlf & Claude Sonnet 4.6 (2026)*

A self-contained HTML viewer showing all 5000 DESI fibers on the focal plane with clickable time plots.

- **Green dots** — active fibers; **orange dots** — Krolewski bad fibers (nsig ≤ −4)
- Hover over a dot to see fiber number, petal, and σ
- Click a dot to display its failure rate vs time plot (with survey average and model prediction)
- Goto box: jump to any fiber by number
- Save plot button: opens current time plot as PNG for saving

### Scripts

| Script | Purpose |
|--------|---------|
| `make_petal_focal_data.py` | Extracts per-petal time series + all focal plane positions into JSON |
| `make_focal_viewer_multi.py` | Builds self-contained HTML viewer for one or more petals |
| `plot_petal_mosaic.py` | 5×100 mosaic of all 500 fibers in a petal (PNG + PDF) |

### Build

```bash
# Extract focal data for all petals (~2 min each)
for p in 0 1 2 3 4 5 6 7 8 9; do
    python make_petal_focal_data.py --petal $p
done

# Interactive viewer (all 10 petals, ~3.9 MB HTML)
python make_focal_viewer_multi.py --petals 0 1 2 3 4 5 6 7 8 9

# Per-petal mosaics (~8 min each, PNG + PDF)
for p in 0 1 2 3 4 5 6 7 8 9; do
    python plot_petal_mosaic.py --petal $p
done
```

### Outputs (written to `/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/all-fibers-vs-time/`)

| File | Description |
|------|-------------|
| `lrg_petal{0..9}_focal_data.json` | Per-petal JSON (~470 KB each) |
| `lrg_petals0_1_2_3_4_5_6_7_8_9_focal_viewer.html` | All-10-petal interactive viewer |
| `lrg_petal{0..9}_mosaic.png` | 5×100 mosaic per petal |
| `lrg_petal{0..9}_mosaic.pdf` | 5×100 mosaic per petal (PDF, ~1.5 MB each) |
