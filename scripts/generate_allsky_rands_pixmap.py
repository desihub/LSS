# script that reads allsky randoms and generates a 
# healpix map with the counts of randoms per pixel
import os
import fitsio 
import numpy as np
from astropy.table import Table
import healpy as hp
import LSS.common_tools as common

randir = "/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/"
outdir = "/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/"
nside = 256
nran  = 18
nest  = False
if nest:
    hp_order = 'nest'
else:
    hp_order = 'ring'
    
def radec2hpix(nside, ra, dec, nest=False):
    hpix = hp.ang2pix(nside, np.radians(90 - dec), np.radians(ra), nest=nest)
    return hpix

def hpixsum(nside, ra, dec, weights=None, nest=False):
    hpix = radec2hpix(nside, ra, dec, nest=nest)
    npix = hp.nside2npix(nside)
    weight_hp = np.bincount(hpix, weights=weights, minlength=npix)
    return weight_hp

ranl = []
for i in range(0,nran):
    print(f"reading allsky randoms {i+1}/{nran}")
    ran = fitsio.read(randir+f'randoms-allsky-1-{i}.fits',columns=['RA','DEC','PHOTSYS'])
    ranl.append(ran)
allsky_rands = np.concatenate(ranl)

for reg in ['N','S']: 
    selr = allsky_rands['PHOTSYS'] == reg
    rallpix = hpixsum(nside,allsky_rands[selr]['RA'],allsky_rands[selr]['DEC'],nest=nest) # count numbers of randoms per pixel

    fnout = outdir+f"allsky_rpix_{reg}_nran{nran}_nside{nside}_{hp_order}.fits"
    tab = Table()
    tab['HPXPIXEL'] = np.arange(12*nside*nside)
    tab['RANDS_HPIX'] = rallpix # random count per pixel
    #tab.write(out_fn)
    common.write_LSS(tab,fnout)
