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


path = '/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky'

#mockpath = 'cutsky_shells_sv3_nz_radecz_xyz_{TARGET}_h5py'

#nametosave = 'mockTargets_{PH}_FirstGen_CutSky_{TARGET}_sv3bits.fits'

types = {'ELG':['z1.100', 131074, 34, 3000], 'LRG':['z0.800', 65537, 1, 3200], 'QSO':['z1.400', 262148, 4, 3400]}
#zname, desi_target, sv3_desi_target, priority

#desi_target=[131074,65537,262148]
#sv3_desi_target=[2,1,4]
#priority=[3000,3200,3400]
#tar=['ELG_SV', 'LRG_SV', 'QSO_SV']
#maxall=0



def mask(main=0, nz=0, Y5=0, sv3=0):
    return main * (2**3) + sv3 * (2**2) + Y5 * (2**1) + nz * (2**0)

#def mask(nz=0, Y5=0, sv3=0):
#    return nz * (2**0) + Y5 * (2**1) + sv3 * (2**2)

file_name = 'cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits'

def permock(i):
    print('Realization {PH}'.format(PH = i))
    fits_tables = []
    sizes = []
    maxall = 0
    output_name = '/global/cscratch1/sd/acarnero/SV3/mockTargets_{PH}_FirstGen_CutSky_alltracers_sv3bits.fits'.format(PH = "%03d" % i)
    for type_ in types: #ii in range(3):
#        print('Starting FirstGen CutSky with {TARGET}'.format(TARGET = type_))
    

        thepath = os.path.join(path, type_, types[type_][0], file_name.format(TYPE = type_, Z = types[type_][0], PH = "%03d" % i)) #mockpath.format(TARGET = types[ii]))
#        print(thepath)

#        file_type = os.path.join(thepath,'cutsky_*_AbacusSummit_base_c000_ph*.fits')
#    print('There are {size} files in {path}'.format(size=len(files_in_path), path = thepath))
#    assert len(files_in_path) == 25, "Different number of files for {TARGET}".format(TARGET=type_)

        f = fits.open(thepath)
        data = f[1].data
#        print('There are %d objects'%len(data['STATUS']))
        
        status = data['STATUS'][()]
        idx = np.arange(len(status))

#        if type_ != 'LRG':
        mask_SV = mask(main=0, nz=1, Y5=0, sv3=1)
#        else:
#            mask_SV = mask(main=1, nz=0, Y5=0, sv3=1)
        
        idx_SV = idx[(status & (mask_SV))==mask_SV]
        data = data[idx_SV]

#        print('In SV3 there are %d objects for realization %3d'%(len(status[(idx_SV)]), i))

#        for key in keys:
#            hdu_data_ez[key] = np.concatenate((hdu_data_ez[key], data[key][list(idx_SV)]), axis=None)
        f.close()
        
        n=len(data)
        targets = Table()
        targets['RA']=data['RA']
        targets['DEC']=data['DEC']
        targets['TRUEZ'] = data['Z_COSMO']
        targets['RSDZ'] = data['Z']
        targets['NZ'] = data['NZ']

        #targets['DESI_TARGET'] = np.zeros(n, dtype='i8')+int(desi_target[ii])
        targets['SV3_DESI_TARGET'] = np.zeros(n, dtype='i8')+int(types[type_][2])
        targets['BGS_TARGET'] = np.zeros(n, dtype='i8')
        targets['MWS_TARGET'] = np.zeros(n, dtype='i8')
        targets['PRIORITY_INIT'] =np.zeros(n, dtype='i8')+int(types[type_][3]) #hdu_mtl['PRIORITY'][mask_elg ][0]
        targets['PRIORITY'] =np.zeros(n, dtype='i8')+int(types[type_][3]) #hdu_mtl['PRIORITY'][mask_elg ][0]
        targets['SUBPRIORITY'] = np.random.uniform(0, 1, n)
        targets['BRICKNAME'] = np.full(n, '000p0000')    #- required !?!
        targets['OBSCONDITIONS'] = np.zeros(n, dtype='i8')+int(3) #hdu_mtl['OBSCONDITIONS'][mask_elg ][0]
        targets['NUMOBS_MORE'] = np.zeros(n, dtype='i8')+int(1) #hdu_mtl['NUMOBS_MORE'][mask_elg ][0]
        targets['NUMOBS_INIT'] = np.zeros(n, dtype='i8')+int(9)
        targets['SV3_SCND_TARGET'] = np.zeros(n, dtype='i8')+int(0)
        targets['ZWARN'] = np.zeros(n, dtype='i8')+int(0)
    #    targets['TRUEZ'] = hdu_data_ez['Z_COSMO']
        #targets['TARGETID'] = np.max(hdu_mtl['TARGETID'])+np.linspace(1,n+1,n,dtype=type(np.max(hdu_mtl['TARGETID'])))
#        targets['TARGETID'] = np.max(hdu_mtl['TARGETID'])+np.max(skyid)+np.linspace(1,n+1,n,dtype=type(np.max(skyid)))+maxall
        targets['TARGETID'] = np.linspace(1,n+1,n,dtype=type(np.max(hdu_mtl['TARGETID'])))+maxall
        #targets['TARGETID'] = np.max(hdu_mtl['TARGETID'])+np.linspace(0,n,n,dtype=type(288230427687173549))
        fits_tables.append(targets)
    #print('Mock'+str(ii)+' DONE')
        maxall += n+1
        sizes.append([type_, n])
    print('sizes: %s = %d, %s = %d, %s = %d'%(sizes[0][0], sizes[0][1], sizes[1][0], sizes[1][1], sizes[2][0], sizes[2][1]))
    new = vstack(fits_tables)
    new.write(output_name, overwrite = True)

    fits.setval(output_name, 'EXTNAME', value='TARGETS', ext=1)
    fits.setval(output_name, 'OBSCON', value='DARK', ext=1)


if __name__ == '__main__':
    par = True
    N = 25
    if par:
        from multiprocessing import Pool
        import sys
        #N = int(sys.argv[2])
        inds = []
        for i in range(0,N):
            inds.append(i)
        with Pool(25) as pool:
            pool.map(permock,inds)
    else:
        for i in range(0,N):
            permock(i)
#permock

