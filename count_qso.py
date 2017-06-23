import os
from astropy.io import fits
from astropy.table import Table, vstack, unique
from astropy.table import Column, join
import numpy as np

tile_info_file = '/project/projectdirs/desi/software/edison/desimodel/master/data/footprint/desi-tiles.fits'
pass_data = fits.open(tile_info_file)[1].data
mask_desi = pass_data['IN_DESI']==1
pass_desi =pass_data[mask_desi]

allqso = fits.open('/scratch1/scratchdirs/angela/data_targets/qso_all_pass.fits')[1].data
tqso = Table(allqso)
qsotableall = join(tqso, pass_desi, keys='TILEID')
qsotableall['index'] = np.arange(len(qsotableall))
first_obs=unique(qsotableall, keys='TARGETID')
array_remove = first_obs['index']
qsotableall.remove_rows([array_remove])
qsotableall['index'][:] = np.arange(len(qsotableall))
second_obs = unique(qsotableall,keys='TARGETID')
array_remove2 = second_obs['index'] 
qsotableall.remove_rows([array_remove2])
qsotableall['index'][:] = np.arange(len(qsotableall))
third_obs = unique(qsotableall, keys='TARGETID')
array_remove3 = third_obs['index']
qsotableall.remove_rows([array_remove3])
qsotableall['index'][:] = np.arange(len(qsotableall))
forth_obs = unique(qsotableall, keys='TARGETID')
array_remove4 = forth_obs['index']
qsotableall.remove_rows([array_remove4])
qsotableall['index'][:] = np.arange(len(qsotableall))
fifth_obs = unique(qsotableall, keys='TARGETID')



outfile = "/scratch1/scratchdirs/angela/data_targets/first_obs_qso_all_pass.fits"
hdu = fits.table_to_hdu(unique_qso)
hdu.writeto(outfile)
