import numpy as np
import healpy as hp
import fitsio
from mockfactory import Catalog
import desimodel.footprint

def radec2thphi(ra,dec):
    return (-dec+90.)*np.pi/180.,ra*np.pi/180.

def cut_randoms(cat_fns, tls_fn):

    tls = fitsio.read(tls_fn)
    cat = Catalog.read(fns)
    inds = desimodel.footprint.is_point_in_desi(tls,cat['RA'], cat['DEC'])
    cat_cut = cat[inds]

    return cat_cut

def map_count(cat, nside = 1024):

    theta,phi = radec2thphi(cat['RA'],cat['DEC'])
    pix = hp.ang2pix(nside,theta,phi,nest=True)
    pix_non_zero,pix_count = np.unique(pix,return_counts=True)
    map_c = np.zeros(12*nside**2)
    map_c[pix_non_zero] = pix_count

    return map_c

nb_ran = 18
tracer = 'QSO'
save_name = ''

target_fns = [f'/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-{idx}.fits' for idx in range(nb_ran)]
fibered_fns = [f'//global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/2/{tracer}zdone_{idx}_full.ran.fits' for idx in range(nb_ran)]

tile_fn = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/everest/DA02tiles-DARK.fits'

target = cut_randoms(target_fns, tile_fn)
fibered = cut_randoms(fibered_fns, tile_fn)

map_fib = map_count(fibered)
map_target = map_count(target)

map_weight = np.divide(map_fib, map_tar, out=np.zeros_like(map_fib), where=map_fib!=0)

hp.write_map(save_name, map_weight, overwrite=True)
