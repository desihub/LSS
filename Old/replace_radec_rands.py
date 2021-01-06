from astropy.io import fits
from astropy.table import Table
from astropy.table import Column, vstack
import numpy as np

mtl_data = fits.open('/global/project/projectdirs/desi/datachallenge/quicksurvey2017/output/dark/0/mtl.fits')[1].data
mask = (mtl_data['DESI_TARGET']==2)
masked_mtl = mtl_data[mask]
nrows_tot = masked_mtl.shape[0]
nrows_tot_all = mtl_data.shape[0]
o_mask = mtl_data['DESI_TARGET']!=2
o_masked_mtl = mtl_data[o_mask]

for file_index in range(0,9):
    print(file_index)
    first_file = '/global/project/projectdirs/desi/datachallenge/LSScat/quicksurvey2016/elg_ran_%d.fits' % file_index
    if file_index < 9:
        extra_data = '/global/project/projectdirs/desi/datachallenge/LSScat/quicksurvey2016/elg_ran_%d.fits' % (file_index+1)
    else:
        extra_data = '/global/project/projectdirs/desi/datachallenge/LSScat/quicksurvey2016/elg_ran_0.fits'
    rd_data = fits.open(first_file)
    rd_extra = fits.open(extra_data)
    nrows1 = rd_data[1].data.shape[0]
    if nrows_tot > nrows1:
        nrows_diff = nrows_tot-nrows1
    hdu = fits.BinTableHDU.from_columns(rd_data[1].columns, nrows=nrows_tot)
    if nrows_tot > nrows1:
        for colname in rd_data[1].columns.names:
            hdu.data[colname][nrows1:]=rd_extra[1].data[colname][0:nrows_diff]

    masked_mtl.OBSCONDITIONS = masked_mtl.OBSCONDITIONS*0
    masked_mtl.OBSCONDITIONS = masked_mtl.OBSCONDITIONS +3
    t = Table([masked_mtl.TARGETID,hdu.data.ra,hdu.data.dec, masked_mtl.DESI_TARGET,masked_mtl.BGS_TARGET,masked_mtl.MWS_TARGET,masked_mtl.SUBPRIORITY,masked_mtl.OBSCONDITIONS,masked_mtl.BRICKNAME,masked_mtl.DECAM_FLUX,masked_mtl.SHAPEDEV_R,masked_mtl.SHAPEEXP_R,masked_mtl.DEPTH_R,masked_mtl.GALDEPTH_R,masked_mtl.NUMOBS_MORE,masked_mtl.PRIORITY],names=( 'TARGETID','RA', 'DEC', 'DESI_TARGET', 'BGS_TARGET','MWS_TARGET', 'SUBPRIORITY', 'OBSCONDITIONS','BRICKNAME','DECAM_FLUX', 'SHAPEDEV_R', 'SHAPEEXP_R', 'DEPTH_R', 'GALDEPTH_R', 'NUMOBS_MORE','PRIORITY'), dtype = ('>i8', '>f8', '>f8', '>i8', '>i8', '>i8', '>f8', '>i2', 'S8', '>f4','>f4', '>f4', '>f4', '>f4', '>i4', '>i8'))

    full_table =vstack([t2, t])
    hdu =fits.table_to_hdu(full_table)

    file_out ='/scratch1/scratchdirs/angela/LSS/LSS_rand_pre_fibassign_target_%d.fits' % file_index
    hdu.writeto(file_out)
