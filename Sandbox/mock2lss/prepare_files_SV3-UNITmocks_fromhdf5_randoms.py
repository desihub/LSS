from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import numpy as np
import glob
import os
import h5py

mtlfile = '/global/project/projectdirs/desi/users/arroyoc/mocks/darksky-v1.0.1-v3-mtlz-files/targets-0-mtlz.fits'
#mtlfile = '/global/cscratch1/sd/mmagana/for_EZmocks/EZ_UNIT/targets-UNIT-mtlz_noELG.fits'

hdu2 = fits.open(mtlfile) # Open the file.
print(hdu2.info())
print(hdu2[1].data.names)
hdu_mtl = hdu2[1].data 
print(np.max(hdu_mtl['TARGETID']))

sky_file='/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/skies.fits'
sky = fits.open(sky_file)
print(sky.info())
hdu_data_sky = sky[1].data 
skyid=hdu_data_sky['TARGETID']

#catalogs=['UNIT_lightcone_multibox_ELG_SV3_footprint_nz_cutsky_snap102.fits',\
#          'UNIT_lightcone_multibox_LRG_SV3_footprint_nz_cutsky_snap102_elg17.fits',\
#          'UNIT_lightcone_multibox_QSO_SV3_footprint_nz_cutsky_snap102.fits']
#catalogs=['UNIT_lightcone_multibox_ELG_SV3_ZL0.6_footprint_nz_cutsky_snap102.fits',\
#          'UNIT_lightcone_multibox_LRG_SV3_footprint_nz_cutsky_snap102_elg17.fits',\
#          'UNIT_lightcone_multibox_QSO_SV3_footprint_nz_cutsky_snap102.fits']
#dirf='/global/project/projectdirs/desi/users/avariu/UNIT/fixedAmp_001_lowres/Gcat-wpmax-v3/sv3_cutsky/'

path = '/global/cscratch1/sd/avariu/desi/ABACUS_2GPC'

mockpath = 'ran_1_cutsky_shells_sv3_nz_radecz_xyz_{TARGET}_h5py' # 'cutsky_shells_sv3_nz_radecz_xyz_{TARGET}_h5py'

nametosave = 'TEST_targets-UNIT-mtlz_{TARGET}_sv3bits_random1.fits'

types = ['ELG', 'LRG', 'QSO']

desi_target=[131074,65537,262148]
sv3_desi_target=[2,1,4]
priority=[3000,3200,3400]
tar=['ELG_SV', 'LRG_SV', 'QSO_SV']
maxall=0

def mask(nz=0, Y5=0, sv3=0):
    return nz * (2**0) + Y5 * (2**1) + sv3 * (2**2)



fits_tables = []


for ii in range(3):
    print('Starting with {TARGET}'.format(TARGET = types[ii]))

    thepath = os.path.join(path, mockpath.format(TARGET = types[ii]))
    print(thepath)

    files_in_path = glob.glob(os.path.join(thepath,'*.h5py'))
    print('There are {size} files in {path}'.format(size=len(files_in_path), path = thepath))

    f0 = h5py.File(files_in_path[0], 'r')
    keys = f0['galaxy'].keys()
    
    hdu_data_ez = {}
    for k in keys:
        hdu_data_ez[k] = np.array([])

    for thefile in files_in_path:
        print(thefile)
        f = h5py.File(thefile, 'r')
        data = f['galaxy']
        print('There are %d objects'%len(data['STATUS']))
        
        status = data['STATUS'][()]
        idx = np.arange(len(status))
        mask_ = mask(nz=1, Y5=0, sv3=1)
        idx_tmp_SV3 = idx[(status & (mask_))==mask_]

        print('In SV3 there are %d objects'%len(status[(idx_tmp_SV3)]))

        for key in keys:
            hdu_data_ez[key] = np.concatenate((hdu_data_ez[key], data[key][(idx_tmp_SV3)]), axis=None)
        f.close()
    '''
    mock_file=dirf+catalogs[ii]
    ez = fits.open(mock_file) # Open the file.
    print(ez.info())
    hdu_data_ez = ez[1].data # Point to the ZCATALOG data.
    print(ez[1].data.names)
    '''

    n=len(hdu_data_ez['RA'])
    targets = Table()
    targets['RA']=hdu_data_ez['RA']
    targets['DEC']=hdu_data_ez['DEC']
    #targets['DESI_TARGET'] = np.zeros(n, dtype='i8')+int(desi_target[ii])
    targets['SV3_DESI_TARGET'] = np.zeros(n, dtype='i8')+int(sv3_desi_target[ii])
    targets['BGS_TARGET'] = np.zeros(n, dtype='i8')
    targets['MWS_TARGET'] = np.zeros(n, dtype='i8')
    targets['PRIORITY'] =np.zeros(n, dtype='i8')+int(priority[ii]) #hdu_mtl['PRIORITY'][mask_elg ][0]
    targets['SUBPRIORITY'] = np.random.uniform(0, 1, n)
    targets['BRICKNAME'] = np.full(n, '000p0000')    #- required !?!
    targets['OBSCONDITIONS'] = np.zeros(n, dtype='i8')+int(3) #hdu_mtl['OBSCONDITIONS'][mask_elg ][0]
    targets['NUMOBS_MORE'] = np.zeros(n, dtype='i8')+int(1) #hdu_mtl['NUMOBS_MORE'][mask_elg ][0]
    targets['NUMOBS_INIT'] = np.zeros(n, dtype='i8')+int(0)
    targets['TRUEZ'] = hdu_data_ez['Z_COSMO']
    #targets['TARGETID'] = np.max(hdu_mtl['TARGETID'])+np.linspace(1,n+1,n,dtype=type(np.max(hdu_mtl['TARGETID'])))
    targets['TARGETID'] = np.max(hdu_mtl['TARGETID'])+np.max(skyid)+np.linspace(1,n+1,n,dtype=type(np.max(skyid)))+maxall
    #targets['TARGETID'] = np.max(hdu_mtl['TARGETID'])+np.linspace(0,n,n,dtype=type(288230427687173549))
    hdu_to_write = targets
    hdu1 = fits.BinTableHDU(data=hdu_to_write, header=hdu2[1].header) # Convert data into BinTableHDU object.
    # Write FITS file to our DESI user directory.
    #hdu.writeto('/global/cscratch1/sd/mmagana/for_EZmocks/EZ_UNIT/mtlz_EZ00'+str(ii)+'.fits', overwrite=True)
    hdu1.writeto(nametosave.format(TARGET = types[ii]),  overwrite=True)
    maxall=n+maxall+1
    print(maxall, n)
    fits_tables.append(nametosave.format(TARGET = types[ii]))
    #print('Mock'+str(ii)+' DONE')


t1 = Table.read(fits_tables[0])
t2 = Table.read(fits_tables[1])
t3 = Table.read(fits_tables[2])

new = vstack([t1, t2, t3])
new.write('targets-UNIT-mtlz_SV3_alltracers_sv3bits_v2_random1.fits')

'''
fits_table_filename1='targets-UNIT-mtlz_ELG_SV_sv3bits.fits'
fits_table_filename2='targets-UNIT-mtlz_LRG_SV_sv3bits.fits'
fits_table_filename3='targets-UNIT-mtlz_QSO_SV_sv3bits.fits'
with fits.open(fits_table_filename1) as hdul1:
    with fits.open(fits_table_filename2) as hdul2:
        with fits.open(fits_table_filename3) as hdul3:
            nrows1 = hdul1[1].data.shape[0]
            nrows2 = hdul2[1].data.shape[0]
            nrows3 = hdul3[1].data.shape[0]
            print(nrows1, nrows2,nrows3)
            nrows = nrows1 + nrows2+nrows3
            hdu = fits.BinTableHDU.from_columns(hdul1[1].columns, nrows=nrows)
            for colname in hdu1.columns.names:
                hdu.data[colname][0:nrows1] = hdul1[1].data[colname]
                hdu.data[colname][nrows1:nrows1+nrows2] = hdul2[1].data[colname]
                hdu.data[colname][nrows1+nrows2:] = hdul3[1].data[colname]
hdu.writeto('targets-UNIT-mtlz_SV3_alltracers_sv3bits.fits', overwrite=True)

'''
