import os
from astropy.io import fits
from astropy.table import Table, vstack, unique
from astropy.table import Column, join
import numpy as np
from random import shuffle

tile_info_file = '/project/projectdirs/desi/software/edison/desimodel/master/data/footprint/desi-tiles.fits'
pass_data = fits.open(tile_info_file)[1].data
mask_desi = pass_data['IN_DESI']==1
pass_desi =pass_data[mask_desi]
#mask_dark = pass_desi['PROGRAM']=='DARK'
#pass_desi_dark = pass_desi[mask_dark]
#len_pass_file = pass_desi_dark.shape[0]
telg = Table( names=('TARGETID','TILEID'))
tlrg = Table( names=('TARGETID','TILEID'))
tqso = Table( names=('TARGETID','TILEID'))
for file_index in range(0, 101):
    infile = '/scratch1/scratchdirs/angela/data_targets/alltargs_data_%d.fits' % file_index
    data = fits.open(infile)[1].data
    maskelg = data['DESI_TARGET']==2
    add_rowselg = Table([data[maskelg]['TARGETID'],data[maskelg]['TILEID']])
    telg=vstack([telg, add_rowselg])

    masklrg = data['DESI_TARGET']==1
    add_rowslrg = Table([data[masklrg]['TARGETID'],data[masklrg]['TILEID']])
    tlrg=vstack([tlrg, add_rowslrg])

    maskqso = data['DESI_TARGET']==4
    add_rowsqso = Table([data[maskqso]['TARGETID'],data[maskqso]['TILEID']])
    tqso=vstack([tqso, add_rowsqso])

#bug in code that leaves TILEID and TARGETID columns empty but adds two new col0:1 instead -> quick work around below
del telg['TARGETID']
del telg['TILEID']
del tlrg['TARGETID']
del tlrg['TILEID']
del tqso['TARGETID']
del tqso['TILEID']

telg['col0'].name = 'TARGETID'
telg['col1'].name = 'TILEID'

tlrg['col0'].name = 'TARGETID'
tlrg['col1'].name = 'TILEID'

tqso['col0'].name = 'TARGETID'
tqso['col1'].name = 'TILEID'
elgtableall = join(telg, pass_desi, keys='TILEID')
elgtableall['index'] = np.arange(len(elgtableall))
outfile = "/scratch1/scratchdirs/angela/data_targets/elg_all_pass.fits"
hdu = fits.table_to_hdu(elgtableall)
hdu.writeto(outfile)
unique_elg=unique(elgtableall, keys='TARGETID')
array_remove = unique_elg['index']
other_elg =elgtableall.remove_rows([array_remove])
outfile = "/scratch1/scratchdirs/angela/data_targets/first_obs_elg_all_pass.fits"
hdu = fits.table_to_hdu(unique_elg)
hdu.writeto(outfile)

lrgtableall = join(tlrg, pass_desi, keys='TILEID')
lrgtableall['index'] = np.arange(len(lrgtableall))
outfile = "/scratch1/scratchdirs/angela/data_targets/lrg_all_pass.fits"
hdu = fits.table_to_hdu(lrgtableall)
hdu.writeto(outfile)
unique_lrg=unique(lrgtableall, keys='TARGETID')
array_remove = unique_lrg['index']
other_lrg =lrgtableall.remove_rows([array_remove])
outfile = "/scratch1/scratchdirs/angela/data_targets/first_obs_lrg_all_pass.fits"
hdu = fits.table_to_hdu(unique_lrg)
hdu.writeto(outfile)

qsotableall = join(tqso, pass_desi, keys='TILEID')
qsotableall['index'] = np.arange(len(qsotableall))
#outfile = "/scratch1/scratchdirs/angela/data_targets/qso_all_pass.fits"
#hdu = fits.table_to_hdu(qsotableall)
#hdu.writeto(outfile)
first_obs=unique(qsotableall, keys='TARGETID')
array_remove = first_obs['index']
second_qso =qsotableall.remove_rows([array_remove])
second_obs = unique(second_qso,keys='TARGETID')
array_remove2 = np.concatenate((second_obs['index'],array_remove) 
thirdqso = qsotableall.remove_rows([array_remove2])
third_obs = unique(thirdqso, keys='TARGETID')
array_remove3 = np.concatenate((third_obs['index'],array_remove2)
forthqso = qsotableall.remove_rows([array_remove3])
forth_obs = unique(forthqso, keys='TARGETID')
array_remove4 = np.concatenate((forth_obs['index'],array_remove3)
fifthqso = qsotableall.remove_rows([array_remove4])
fifth_obs = unique(fifthqso, keys='TARGETID')

outfile = "/scratch1/scratchdirs/angela/data_targets/first_obs_qso_all_pass.fits"
hdu = fits.table_to_hdu(unique_qso)
hdu.writeto(outfile)
