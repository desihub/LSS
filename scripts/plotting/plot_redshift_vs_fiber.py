# Make 2-D density plots of redshift vs fiber
# Example:
# python plot_redshift_vs_fiber.py --tracer LRG --plot_dir /pscratch/sd/r/rongpu/tmp/

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc, argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

params = {'figure.facecolor': 'w'}
plt.rcParams.update(params)

parser = argparse.ArgumentParser()
parser.add_argument("--tracer", help="tracer type (BGS, LRG, ELG, or QSO)", required=True)
parser.add_argument("--fn", help="path of the catalog file", default=None)
parser.add_argument("--plot_dir", help="directory to save the plots", default=None)

args = parser.parse_args()
tracer = args.tracer.upper()
fn = args.fn
plot_dir = args.plot_dir

z_bins_dict = {'BGS': np.arange(-0.025, 0.7, 0.02), 'LRG': np.arange(-0.025, 1.5, 0.025), 'ELG': np.arange(0.6, 1.8, 0.025), 'QSO': np.arange(-0.025, 4.5, 0.1)}
vmax_dict = {'BGS': 1/3e4, 'LRG': 1/4e4, 'ELG': 1/6e4, 'QSO': 1/4e4}

if fn is None:
    fn_dict = {'BGS': 'BGS_ANY_full_noveto.dat.fits', 'LRG': 'LRG_full_noveto.dat.fits', 'ELG': 'ELG_full_noveto.dat.fits', 'QSO': 'QSO_full_noveto.dat.fits'}
    fn = os.path.join('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/daily/LSScats/test/', fn_dict[tracer])

min_nobs = 100

cat = Table(fitsio.read(fn))
print(len(cat))

if 'Z_not4clus' in cat.colnames:
    cat.rename_column('Z_not4clus', 'Z')

cat['EFFTIME_BGS'] = 0.1400 * cat['TSNR2_BGS']
cat['EFFTIME_LRG'] = 12.15 * cat['TSNR2_LRG']

# Remove FIBERSTATUS!=0 fibers
mask = cat['COADD_FIBERSTATUS']==0
print('FIBERSTATUS   ', np.sum(~mask), np.sum(mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# GOODHARDLOC
mask = cat['GOODHARDLOC']==True
print('GOODHARDLOC   ', np.sum(~mask), np.sum(mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# Remove "no data" fibers
mask = cat['ZWARN'] & 2**9==0
print('No data   ', np.sum(~mask), np.sum(mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# Require a minimum depth for the cat coadd
if tracer=='BGS':
    min_depth = 160
    mask = cat['EFFTIME_BGS']>min_depth
else:
    min_depth = 800.
    mask = cat['EFFTIME_LRG']>min_depth
print('Min depth   ', np.sum(~mask), np.sum(mask), np.sum(~mask)/len(mask))
cat = cat[mask]

if tracer=='LRG':
    # Apply LRG mask
    mask = cat['lrg_mask']==0
    print('LRG mask', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]
elif tracer=='ELG':
    # Apply maskbits
    maskbits = [1, 11, 12, 13]
    mask = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask &= (cat['MASKBITS'] & 2**bit)==0
    print('MASKBITS  ', np.sum(~mask), np.sum(mask), np.sum(~mask)/len(mask))
    cat = cat[mask]

if tracer=='ELG':
    # Redshift quality cut
    cat['q'] = (cat['OII_FLUX']>0) & (cat['OII_FLUX_IVAR']>0)
    cat['q'] &= np.log10(cat['OII_FLUX'] * np.sqrt(cat['OII_FLUX_IVAR'])) > 0.9 - 0.2 * np.log10(cat['DELTACHI2'])
    print('ELG redshift quality', np.sum(~cat['q'])/len(cat))

print(len(cat))

fiberstats = Table()
fiberstats['FIBER'], fiberstats['n_tot'] = np.unique(cat['FIBER'], return_counts=True)
fiberstats['weight'] = 1/fiberstats['n_tot']*np.median(fiberstats['n_tot'])
cat = join(cat, fiberstats[['FIBER', 'weight']], join_type='outer')
too_few_fibers = fiberstats['FIBER'][fiberstats['n_tot']<=min_nobs]

fig, ax = plt.subplots(10, 1, figsize=(16, 20))
for index in range(10):
    fiber_min, fiber_max = index*500-0.5, (index+1)*500-0.5
    mask = (cat['FIBER']>fiber_min) & (cat['FIBER']<fiber_max)
    xbins, ybins = np.linspace(fiber_min, fiber_max, 501), z_bins_dict[tracer]
    ybins = ybins - np.diff(ybins)[0]/2  # center one of the bins at z=0
    ax[index].hist2d(cat['FIBER'][mask], cat['Z'][mask], bins=[xbins, ybins], vmin=0, vmax=len(cat)*vmax_dict[tracer])
    ax[index].set_xlim(fiber_min, fiber_max)
plt.tight_layout()
if plot_dir is not None:
    plt.savefig(os.path.join(plot_dir, 'redshift_vs_fiber_{}.png'.format(tracer.lower())))
plt.show()

# Normalize by the number of objects in each fiber
fig, ax = plt.subplots(10, 1, figsize=(16, 20))
for index in range(10):
    fiber_min, fiber_max = index*500-0.5, (index+1)*500-0.5
    mask = (cat['FIBER']>fiber_min) & (cat['FIBER']<fiber_max)
    mask &= ~np.in1d(cat['FIBER'], too_few_fibers)  # Do not plot fibers with too few objects
    xbins, ybins = np.linspace(fiber_min, fiber_max, 501), z_bins_dict[tracer]
    ybins = ybins - np.diff(ybins)[0]/2  # center one of the bins at z=0
    ax[index].hist2d(cat['FIBER'][mask], cat['Z'][mask], weights=cat['weight'][mask], bins=[xbins, ybins], vmin=0, vmax=len(cat)*vmax_dict[tracer])
    ax[index].set_xlim(fiber_min, fiber_max)
plt.tight_layout()
if plot_dir is not None:
    plt.savefig(os.path.join(plot_dir, 'redshift_vs_fiber_{}_norm.png'.format(tracer.lower())))
plt.show()

if tracer=='ELG':

    fig, ax = plt.subplots(10, 1, figsize=(16, 20))
    for index in range(10):
        fiber_min, fiber_max = index*500-0.5, (index+1)*500-0.5
        mask = (cat['FIBER']>fiber_min) & (cat['FIBER']<fiber_max)
        mask &= cat['q']
        xbins, ybins = np.linspace(fiber_min, fiber_max, 501), z_bins_dict[tracer]
        ybins = ybins - np.diff(ybins)[0]/2  # center one of the bins at z=0
        ax[index].hist2d(cat['FIBER'][mask], cat['Z'][mask], bins=[xbins, ybins], vmin=0, vmax=len(cat)*vmax_dict[tracer])
        ax[index].set_xlim(fiber_min, fiber_max)
    plt.tight_layout()
    if plot_dir is not None:
        plt.savefig(os.path.join(plot_dir, 'redshift_vs_fiber_{}_goodz.png'.format(tracer.lower())))
    plt.show()

    # Normalize by the number of objects in each fiber
    fig, ax = plt.subplots(10, 1, figsize=(16, 20))
    for index in range(10):
        fiber_min, fiber_max = index*500-0.5, (index+1)*500-0.5
        mask = (cat['FIBER']>fiber_min) & (cat['FIBER']<fiber_max)
        mask &= ~np.in1d(cat['FIBER'], too_few_fibers)  # Do not plot fibers with too few objects
        mask &= cat['q']
        xbins, ybins = np.linspace(fiber_min, fiber_max, 501), z_bins_dict[tracer]
        ybins = ybins - np.diff(ybins)[0]/2  # center one of the bins at z=0
        ax[index].hist2d(cat['FIBER'][mask], cat['Z'][mask], weights=cat['weight'][mask], bins=[xbins, ybins], vmin=0, vmax=len(cat)*vmax_dict[tracer])
        ax[index].set_xlim(fiber_min, fiber_max)
    plt.tight_layout()
    if plot_dir is not None:
        plt.savefig(os.path.join(plot_dir, 'redshift_vs_fiber_{}_goodz_norm.png'.format(tracer.lower())))
    plt.show()

