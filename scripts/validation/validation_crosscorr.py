import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import fitsio
import healpy as hp
import treecorr
from astropy.table import join,Table,hstack

# import sys
# sys.path.append("/pscratch/sd/s/shreeb/shreeb/LSS/py")
import LSS.common_tools as common

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="Base directory for catalogs", default='/dvs_ro/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--outdir", help="Output directory for plots", default=None)
parser.add_argument("--version", help="Catalog version", default='test')
parser.add_argument("--extra_clusdir", help="extra directory in path for finding clustering catalogs", default='/nonKP')
parser.add_argument("--survey", help="Survey name (e.g. Y1, DA02)", default='Y1')
parser.add_argument("--data", help="Data type (LSS or mock)", default='LSS')
parser.add_argument("--verspec", help="Spectroscopic version", default='iron')
parser.add_argument("--tracers", help="Tracer type or 'all'", default='ELG_LOPnotqso')  # Currently works for only 1 tracer. 'all' will only return results for LRG. Needs to be updated later.
parser.add_argument("--nran",help="number of random files to use",default=1,type=int)
parser.add_argument("--sys_wts", help="Whether to use imaging systematic weights", default=False,type=bool)
parser.add_argument("--mapmd", help="set of maps to use",default='default')
parser.add_argument("--norm", help="whether to normalize the maps before cross correlation",default=False,type=bool)

args = parser.parse_args()

# --- Derived paths ---
indir = os.path.join(args.basedir, args.survey, args.data, args.verspec, 'LSScats',args.version)
outdir = args.outdir or os.path.join(indir, 'plots/imaging/').replace('dvs_ro', 'global')#.replace(args.extra_clusdir,'')
os.makedirs(outdir, exist_ok=True)

# --- Tracer list ---
if args.tracers == 'all':
    tracers = ['LRG', 'ELG_LOPnotqso', 'QSO', 'BGS_BRIGHT-21.5']
else:
    tracers = [args.tracers]

nside = 256
nest = True

# ---------- Map categories ----------
all_default_map_names = [
    'EBV_CHIANG_SFDcorr',
    'STARDENS',
    'HALPHA',
    'HALPHA_ERROR',
    'CALIB_G', 'CALIB_R', 'CALIB_Z',
    'EBV_MPF_Mean_FW15',
    'EBV_MPF_Mean_ZptCorr_FW15',
    'EBV_MPF_Var_FW15',
    'EBV_MPF_VarCorr_FW15',
    'EBV_MPF_Mean_FW6P1',
    'EBV_MPF_Mean_ZptCorr_FW6P1',
    'EBV_MPF_Var_FW6P1',
    'EBV_MPF_VarCorr_FW6P1',
    'EBV_SGF14',
    'BETA_ML',
    'BETA_MEAN',
    'BETA_RMS',
    'HI',
    'KAPPA_PLANCK',
    'EBV',
    'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z',
    'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z',
    'PSFDEPTH_W1', 'PSFDEPTH_W2',
    'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z',
    'EBV_DIFF_GR',
    'EBV_DIFF_RZ'
]

default_map_names = [
    'STARDENS', 'HI',
    'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z',
    'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z',
    'EBV_DIFF_GR', 'EBV_DIFF_RZ',
    'EBV', 'EBV_MPF_Var_FW15'
]

extradefault_map_names = [m for m in all_default_map_names if m not in default_map_names]

special_map_names = ['lrg_mask_frac', 'sagittarius', 'sky_g', 'sky_r', 'sky_z']
extraspecial_map_names = ['stellar_density_0_10', 'stellar_density_10_14', 'stellar_density_14_18', 'stellar_density_18_22'] + \
                         ['OII_3727', 'Hbeta_4861', 'OIII_4959', 'OIII_5007', 'NeIII_3869']

# ---------- Map containers ----------
default_maps = {}
extradefault_maps = {}
special_maps = {}
extraspecial_maps = {}

# ---------- Load default maps ----------
if args.mapmd in ['all', 'default']:
    filename_base = os.path.join(indir, 'hpmaps', f'{tracers[0]}_mapprops_healpix_nested_nside256')
    mf = {
        'N': hstack([Table.read(f'{filename_base}_N.fits'), common.get_debv()]),
        'S': hstack([Table.read(f'{filename_base}_S.fits'), common.get_debv()])
    }
    for mp in default_map_names:
        default_maps[mp] = {
            'N': mf['N'][mp],
            'S': mf['S'][mp]
        }

# ---------- Load EXTRA-default maps ----------
if args.mapmd in ['all', 'extradefault']:
    filename_base = os.path.join(indir, 'hpmaps', f'{tracers[0]}_mapprops_healpix_nested_nside256')
    mf = {
        'N': hstack([Table.read(f'{filename_base}_N.fits'), common.get_debv()]),
        'S': hstack([Table.read(f'{filename_base}_S.fits'), common.get_debv()])
    }
    for mp in extradefault_map_names:
        extradefault_maps[mp] = {
            'N': mf['N'][mp],
            'S': mf['S'][mp]
        }

# ---------- Load special maps ----------
if args.mapmd in ['all', 'special']:
    # Sky residuals maps
    def load_sky_residual_band(band):
        map_n = np.zeros(256*256*12)
        map_s = np.zeros(256*256*12)

        f = fitsio.read(f'/dvs_ro/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_256_north.fits')
        pix_nest = hp.ring2nest(256, f['HPXPIXEL'])
        map_n[pix_nest] = f[f'sky_median_{band}']

        f = fitsio.read(f'/dvs_ro/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_256_south.fits')
        pix_nest = hp.ring2nest(256, f['HPXPIXEL'])
        map_s[pix_nest] = f[f'sky_median_{band}']

        return {'N': map_n, 'S': map_s}

    for band in ['g', 'r', 'z']:
        special_maps[f'sky_{band}'] = load_sky_residual_band(band)

    # LRG mask fraction
    lrg_mask_frac = np.zeros(256*256*12)
    ranmap = np.zeros(256*256*12)
    ranmap_lmask = np.zeros(256*256*12)
    randir = '/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
    ran = fitsio.read(randir + 'randoms-1-0.fits', columns=['RA', 'DEC'])
    ran_lrgmask = fitsio.read('/dvs_ro/cfs/cdirs/desi/survey/catalogs/main/LSS/randoms-1-0lrgimask.fits')
    th, phi = np.radians(90 - ran['DEC']), np.radians(ran['RA'])
    ranpix = hp.ang2pix(256, th, phi, nest=True)
    for pix, mval in zip(ranpix, ran_lrgmask['lrg_mask']):
        ranmap[pix] += 1
        if mval > 1:
            ranmap_lmask[pix] += 1
    sel = ranmap > 0
    lrg_mask_frac[sel] = ranmap_lmask[sel] / ranmap[sel]
    special_maps['lrg_mask_frac'] = {'S': lrg_mask_frac,'N': lrg_mask_frac}

    # Sagittarius stream
    sag = np.load('/dvs_ro/cfs/cdirs/desi/survey/catalogs/extra_regressis_maps/sagittarius_stream_256.npy')
    special_maps['sagittarius'] = {'S':sag}

# ---------- Load EXTRA-special maps ----------
if args.mapmd in ['all', 'extraspecial']:
    # Stellar density maps
    star_bins = ['0_10', '10_14', '14_18', '18_22']
    for b in star_bins:
        data = np.load(f'/global/cfs/cdirs/desi/survey/catalogs/external_input_maps/stardens/stellar_density_maps_smoothed/stellar_density_map_data_{b}_smoothed.npy')
        extraspecial_maps[f'stellar_density_{b}'] = {'N': data,'S': data}

    # MWS emission line maps
    mws_maps = ['OII_3727', 'Hbeta_4861', 'OIII_4959', 'OIII_5007', 'NeIII_3869']
    mwsf = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/mws_emline_maps/MWS_emission_line_fluxes_combined.fits')
    for mp in mws_maps:
        extraspecial_maps[mp] = {'N': mwsf[mp],'S': mwsf[mp]}

# ---------- Combine all maps ----------
all_maps = {}
if args.mapmd in ['all', 'default']:
    all_maps.update(default_maps)
if args.mapmd in ['all', 'extradefault']:
    all_maps.update(extradefault_maps)
if args.mapmd in ['all', 'special']:
    all_maps.update(special_maps)
if args.mapmd in ['all', 'extraspecial']:
    all_maps.update(extraspecial_maps)

for name, val in all_maps.items():
    if isinstance(val, dict):
        for reg, v in val.items():
            print(f"[CHECK] {name} region {reg} → type: {type(v)}, shape: {getattr(v, 'shape', 'no shape')}")
    else:
        print(f"[ERROR] {name} is NOT a dict: type={type(val)}, value={val}")


#Define functions
def get_region_pixels(data, nside=256, nest=True,reg_split='NS'):
    """
    Get pixel indices for N and S regions using PHOTSYS and DES mask.

    Returns:
        pix_north, pix_south : arrays of pixel indices for N and S
    """
    
    data_pix = hp.ang2pix(nside, data['RA'], data['DEC'], nest=nest, lonlat=True)

    # Load DES mask pixels
    mode = 'nest' if nest else 'ring'
    des_mask_path = f'/global/cfs/cdirs/desi/users/rongpu/useful/in_des/hp_in_des_{nside}_{mode}.fits.gz'
    des_mask = Table.read(des_mask_path)['in_des']
    des_pixels = np.where(des_mask)[0]

    # Region masks
    mask_bm = data['PHOTSYS'] == 'N'
    mask_des = np.in1d(data_pix, des_pixels)
    mask_decals = (~mask_bm) & (~mask_des)
    
    # Pixels for each region
    pix_north = data_pix[mask_bm]
    pix_south = data_pix[~mask_bm]
    pix_des = data_pix[mask_des]  # strictly DES here
    pix_decals = data_pix[mask_decals]
    if reg_split == 'NS':
        return pix_north, pix_south, mask_bm, ~mask_bm
    if reg_split == 'NSdes':
        return pix_north, pix_decals, pix_des, mask_bm, mask_decals,mask_des


def compute_overdensity_north_south(
    catalog_dir,
    data_filename=args.tracers+'_clustering.dat.fits',
    random_prefix=args.tracers+'_',
    random_suffix='_clustering.ran.fits',
    n_randoms=18,
    nside=256,
    min_frac_area=0.2,
    sys_wts = False
):
    """
    Compute overdensity maps in North and South using PHOTSYS-defined regions.
    """

    npix = hp.nside2npix(nside)
    data_counts_n = np.zeros(npix)
    data_counts_s = np.zeros(npix)
    rand_counts_n = np.zeros(npix)
    rand_counts_s = np.zeros(npix)
    allsky_counts_n = np.zeros(npix)
    allsky_counts_s = np.zeros(npix)

    # ---- Load data ----
    data_path = os.path.join(catalog_dir, data_filename)
    print(f"Reading data: {data_path}")
    data = fitsio.read(data_path)

    # Region-based pixels from PHOTSYS
    pix_north, pix_south, mask_n, mask_s = get_region_pixels(data, nside=nside, nest=True)

    weights_data = data['WEIGHT'] if 'WEIGHT' in data.dtype.names else np.ones_like(data['RA'])
    if sys_wts == False: weights_data = weights_data/data['WEIGHT_SYS']

    # Accumulate data counts (sum of weighted galaxies) - used for overdensity computation
    np.add.at(data_counts_n, pix_north, weights_data[mask_n])
    np.add.at(data_counts_s, pix_south, weights_data[mask_s])

    # ---- Load randoms and accumulate ----
    for i in range(n_randoms):
        rand_file = os.path.join(catalog_dir, f"{random_prefix}{i}{random_suffix}")
        # print(f"Reading random: {rand_file}")
        rand = fitsio.read(rand_file)
    
        rand_pix = hp.ang2pix(nside, rand['RA'], rand['DEC'], nest=True, lonlat=True)
        weights_rand = rand['WEIGHT'] if 'WEIGHT' in rand.dtype.names else np.ones_like(rand['RA'])
        if sys_wts == False: weights_rand = weights_rand/rand['WEIGHT_SYS']

        # Apply region pixel mask (data mask to randoms)
        in_n = np.in1d(rand_pix, pix_north)
        in_s = np.in1d(rand_pix, pix_south)

        # sum of weighted randoms - used for overdensity computation
        np.add.at(rand_counts_n, rand_pix[in_n], weights_rand[in_n])
        np.add.at(rand_counts_s, rand_pix[in_s], weights_rand[in_s])
        # sum of unweighted randoms - used for fractional area masking
        np.add.at(allsky_counts_n, rand_pix[in_n], 1)
        np.add.at(allsky_counts_s, rand_pix[in_s], 1)

    print("Available columns in data file:")
    print(data.dtype.names)
    print("Available columns in randoms file:")
    print(rand.dtype.names)

    # ---- Compute mask ----
    frac_area_n = rand_counts_n / (allsky_counts_n + 1e-10)
    frac_area_s = rand_counts_s / (allsky_counts_s + 1e-10)

    mask_n = frac_area_n > min_frac_area
    mask_s = frac_area_s > min_frac_area

    # print(f"Sum rand_counts_n (masked): {np.sum(rand_counts_n[mask_n])}")
    # print(f"Sum rand_counts_s (masked): {np.sum(rand_counts_s[mask_s])}")

    # ---- Compute overdensity ----
    mean_dataran_n = np.sum(data_counts_n[mask_n]) / np.sum(rand_counts_n[mask_n])
    mean_dataran_s = np.sum(data_counts_s[mask_s]) / np.sum(rand_counts_s[mask_s])

    overdensity_n = ((data_counts_n / (rand_counts_n + 1e-10)) / mean_dataran_n) - 1
    overdensity_s = ((data_counts_s / (rand_counts_s + 1e-10)) / mean_dataran_s) - 1

    # ---- Mask result ----
    overdensity_n_masked = np.full_like(overdensity_n, hp.UNSEEN)
    overdensity_s_masked = np.full_like(overdensity_s, hp.UNSEEN)
    overdensity_n_masked[mask_n] = overdensity_n[mask_n]
    overdensity_s_masked[mask_s] = overdensity_s[mask_s]

    return overdensity_n_masked, overdensity_s_masked, mask_n, mask_s


def compute_overdensity(
    catalog_dir,
    data_filename=args.tracers+'_clustering.dat.fits',
    random_prefix=args.tracers+'_',
    random_suffix='_clustering.ran.fits',
    n_randoms=18,
    nside=256,
    min_frac_area=0.2,
    sys_wts = False
):
    """Compute overdensity map with WEIGHT_COMP and optional imaging weight."""

    npix = hp.nside2npix(nside)
    data_counts = np.zeros(npix)
    rand_counts = np.zeros(npix)
    allsky_rand_counts = np.zeros(npix)

    # --- 1. Load data ---
    data = fitsio.read(os.path.join(catalog_dir, data_filename))
    ra_data, dec_data = data['RA'], data['DEC']

    weights_data = data['WEIGHT'] if 'WEIGHT' in data.dtype.names else np.ones_like(data['RA'])
    if sys_wts == False: weights_data = weights_data/data['WEIGHT_SYS']

    pix_data = hp.ang2pix(nside, ra_data, dec_data, lonlat=True)
    np.add.at(data_counts, pix_data, weights_data)

    # --- 2. Load and accumulate randoms ---
    for i in range(n_randoms):
        rand_file = os.path.join(catalog_dir, f"{random_prefix}{i}{random_suffix}")
        rand = fitsio.read(rand_file)
    
        ra_rand, dec_rand = rand['RA'], rand['DEC']
        rand_pix = hp.ang2pix(nside, ra_rand, dec_rand, lonlat=True)
    
        weights_rand = rand['WEIGHT'] if 'WEIGHT' in rand.dtype.names else np.ones_like(rand['RA'])
        if sys_wts == False: weights_rand = weights_rand/rand['WEIGHT_SYS']
    
        # Accumulate weighted and unweighted randoms
        np.add.at(rand_counts, rand_pix, weights_rand)
        np.add.at(allsky_rand_counts, rand_pix, 1)

    print("Available columns in data file:")
    print(data.dtype.names)
    print("Available columns in randoms file:")
    print(rand.dtype.names)

    # --- 3. Apply fractional area mask ---
    frac_area = rand_counts / (allsky_rand_counts + 1e-10)
    mask = frac_area > min_frac_area

    # --- 4. Mean density (masked) ---
    mean_dataran = np.sum(data_counts[mask]) / np.sum(rand_counts[mask])

    # --- 5. Overdensity map ---
    dat_over_ran = data_counts / (rand_counts + 1e-10)
    overdensity = (dat_over_ran / mean_dataran) - 1

    # Apply mask
    overdensity_masked = np.full_like(overdensity, hp.UNSEEN)
    overdensity_masked[mask] = overdensity[mask]

    return overdensity_masked, mask

def compute_cross_correlation(dens_map, sys_map, label, ax=None, plot=True, color = None, normalize=False, valid_mask = None):
    """
    Compute and plot cross-correlation between tracer overdensity and a systematics map.
    This version uses strict masking and avoids TreeCorr seeing any UNSEEN values.
    """
        
    assert dens_map.ndim == 1, f"density map is not 1D: shape = {dens_map.shape}"
    assert sys_map.ndim == 1, f"systematics map is not 1D: shape = {sys_map.shape}"

    nside = hp.get_nside(dens_map)
    npix = hp.nside2npix(nside)

    # Get RA/DEC for all pixels
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    ra = np.degrees(phi)
    dec = 90 - np.degrees(theta)

    # Valid pixels only (exclude UNSEEN and NaNs)
    if valid_mask is None:
        valid = (
            (dens_map != hp.UNSEEN) & np.isfinite(dens_map) &
            (sys_map != hp.UNSEEN) & np.isfinite(sys_map)
        )
    else:
        valid = valid_mask

    # Extract only valid pixel positions and values
    ra_valid = ra[valid]
    dec_valid = dec[valid]
    dens_valid = dens_map[valid]
    sys_valid = sys_map[valid]

    if normalize:
        dens_valid = (dens_valid - np.mean(dens_valid)) / (np.std(dens_valid) + 1e-10)
        sys_valid = (sys_valid - np.mean(sys_valid)) / (np.std(sys_valid) + 1e-10)

    # Build TreeCorr catalogs from just valid pixels
    cat1 = treecorr.Catalog(ra=ra_valid, dec=dec_valid, k=dens_valid,
                            ra_units='deg', dec_units='deg')
    cat2 = treecorr.Catalog(ra=ra_valid, dec=dec_valid, k=sys_valid,
                            ra_units='deg', dec_units='deg')

    # Correlate
    kk = treecorr.KKCorrelation(nbins=20, min_sep=0.2, max_sep=10, sep_units='degrees')
    kk.process(cat1, cat2)

    if plot and ax is not None:
        ax.plot(kk.meanr, kk.xi, label=label, color = color, marker='o')
        # ax.set_xscale("log")
        ax.set_xlabel(r"$\theta$ (degrees)", fontsize = 14)
        ax.set_ylabel(r"Cross-Correlation w$_{fg}$($\theta$)", fontsize = 14)
        ax.legend(fontsize = 15)
            
    return kk.meanr, kk.xi


def estimate_A_from_R(theta, R, tmin=0.5, tmax=10, method='median'):
    """
    Estimate amplitude A from R(θ) = w_fg / w_ff over a stable θ range.

    Args:
        theta : array of angular separations (degrees)
        R     : array of R(θ) values
        tmin, tmax : angular range (in degrees) for averaging
        method : 'median' or 'mean'

    Returns:
        A : estimated contamination amplitude
    """
    mask = (theta >= tmin) & (theta <= tmax)
    if not np.any(mask):
        raise ValueError(r"No $\theta$ values within specified range")

    R_selected = R[mask]

    if method == 'mean':
        return np.mean(R_selected)
    elif method == 'median':
        return np.median(R_selected)
    else:
        raise ValueError("method must be 'mean' or 'median'")


def compute_sys_overdensity_map(sys_map, mask, unseen_fill=True):
    """
    Convert a systematics map to overdensity: δ = map / ⟨map⟩ - 1
    using the masked pixels.

    Parameters
    ----------
    sys_map : np.ndarray
        Input HEALPix systematics map (1D array of shape [npix])
    mask : np.ndarray of bool
        Boolean mask indicating valid pixels (same length as sys_map)
    unseen_fill : bool
        If True, return map with hp.UNSEEN outside mask; else leave unchanged

    Returns
    -------
    overdensity_map : np.ndarray
        Overdensity version of the input map
    """

    if not isinstance(sys_map, np.ndarray):
        raise TypeError("Systematics map must be a NumPy array.")
    if sys_map.shape != mask.shape:
        raise ValueError("Map and mask must have the same shape.")

    valid = mask & np.isfinite(sys_map) & (sys_map != hp.UNSEEN)

    if not np.any(valid):
        raise ValueError("No valid pixels found for overdensity calculation.")

    mean_value = np.mean(sys_map[valid])
    overdensity = (sys_map / mean_value) - 1.0

    if unseen_fill:
        overdensity_masked = np.full_like(sys_map, hp.UNSEEN)
        overdensity_masked[valid] = overdensity[valid]
        return overdensity_masked
    else:
        return overdensity

# ----------- Loop over Tracers -------------------
colors = {'N': 'blue', 'S': 'red', 'NS': 'green'}

for tracer in tracers:
    print(f"\n🔍   Running validation for tracer: {tracer}")
    catalog_file = f"{tracer}_clustering.dat.fits"

    # Compute overdensity tracer maps
    dens_n, dens_s, mask_n, mask_s = compute_overdensity_north_south(indir+args.extra_clusdir, data_filename=catalog_file, random_prefix=f'{tracer}_', n_randoms= args.nran, sys_wts = args.sys_wts)
    #dens_ns, mask_ns = compute_overdensity(indir+args.extra_clusdir, data_filename=catalog_file, random_prefix=f'{tracer}_', n_randoms= args.nran, sys_wts = args.sys_wts)
    dens_map = {'N': dens_n, 'S': dens_s}#, 'NS': dens_ns}
    masks = {'N': mask_n, 'S': mask_s}#, 'NS': mask_ns}

    for region in ['N', 'S']:#, 'NS']:
        if region not in dens_map: continue

        tracer_map = dens_map[region]
        region_maps_to_plot = {k: v[region] for k, v in all_maps.items() if region in v}
        if not region_maps_to_plot:
            print(f"  ⚠ No systematics maps found for region {region} — skipping.")
            continue

        nplots = len(region_maps_to_plot)
        fig, axes = plt.subplots(nrows=nplots, ncols=3, figsize=(18, 5 * nplots), squeeze=False)
        row_idx = 0

        for sys_name, sys_map in region_maps_to_plot.items():

            # Convert to overdensity using region masks from tracer density
            sys_map = compute_sys_overdensity_map(sys_map, masks[region])
            
            print(f"🔎   Running {sys_name} | sys_map shape: {sys_map.shape}")

            try:
                valid_mask = ((tracer_map != hp.UNSEEN) & np.isfinite(tracer_map) & (sys_map != hp.UNSEEN) & np.isfinite(sys_map))
                
                # Cross-correlation w_fg
                theta_fg, w_fg = compute_cross_correlation(tracer_map, sys_map, label = sys_name, ax=axes[row_idx][0], plot=True, color = colors[region], normalize=args.norm, valid_mask=valid_mask)

                # Auto-correlation w_ff
                theta_ff, w_ff = compute_cross_correlation(sys_map, sys_map, f"{sys_name} auto", plot=False, normalize=args.norm, valid_mask=valid_mask)

                # R(theta) and Systematic Impact
                if np.allclose(theta_fg, theta_ff):
                    epsilon = 1e-50
                    R_theta = w_fg / (w_ff + epsilon)
                    A = estimate_A_from_R(theta_fg, R_theta)

                    # R plot
                    ax1 = axes[row_idx][1]
                    ax1.plot(theta_fg, R_theta, marker='o', color = colors[region])
                    # ax1.set_xscale("log")
                    ax1.set_xlabel(r"$\theta$ (degrees)", fontsize = 14)
                    ax1.set_ylabel(r"R($\theta$) = w$_{fg}$/w$_{ff}$", fontsize = 14)

                    # A² × w_ff
                    ax2 = axes[row_idx][2]
                    ax2.plot(theta_fg, A**2 * w_ff, marker='o', color = colors[region], label = f'A = {round(A,4)}')
                    ax2.plot(theta_fg, R_theta**2 * w_ff,':k', label = r'using R($\theta$) instead of A')
                    # ax2.set_xscale("log")
                    ax2.set_xlabel(r"$\theta$ (degrees)", fontsize = 14)
                    ax2.set_ylabel(r"Systematic Impact = A$^2$ × w$_{ff}$", fontsize = 14)
                    ax2.legend(fontsize = 15)

                else:
                    print(f"[WARN] $\\theta$ bins do not match — skipping R/impact plot for {sys_name}")
                    axes[row_idx][1].set_visible(False)
                    axes[row_idx][2].set_visible(False)

            except Exception as e:
                print(f"    ❌ Failed for {sys_name}: {e}")
                for i in range(3):
                    axes[row_idx][i].set_visible(False)

            row_idx += 1

        # Hide unused rows
        for r in range(row_idx, nplots):
            for c in range(3):
                axes[r][c].set_visible(False)

        if args.norm == True:
            normtype = 'norm_'
        else:
            normtype = ''
        str_wt = '_nowsys'
        if args.sys_wts:
            str_wt = '_wsys'
        pdf_file = os.path.join(
            outdir, f'crosscorr_{tracer}_{region}_{args.mapmd}_{normtype}full_validation{str_wt}.pdf'
        )
        axes[0,1].set_title(f"{tracer} - {region}", fontsize = 20)
        fig.tight_layout()
        fig.savefig(pdf_file)
        plt.close(fig)
        print(f"✅    Saved combined validation plot: {pdf_file}")

print(f"🎉    Finished all regions for tracer {tracer}")

