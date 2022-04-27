import astropy.io.fits as pf
import os
import numpy as np

types_ = ['LRG','ELG','QSO']
d_bite = [1,34,4]
r_bit = [1,2,4]


path_mock = '/global/cscratch1/sd/acarnero/SV3/LSS_MTL_rea000_univ1/fuji/LSScats/test'
path_data = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/fuji/LSScats/3'
path_input = '/global/cscratch1/sd/acarnero/SV3'

area_sv3 = 207.5

nran = 18
SNUMS = np.linspace(100, 5000, num=50, dtype=np.int)



for j,ty in enumerate(types_):
    print('Numbers for ',ty)
    #Check density data as Ngal/Nran*2500 in the full catalog
    datacat = os.path.join(path_data, ty+'_full.dat.fits')
    size_datacat = float(len(pf.open(datacat)[1].data))
    dens_datacat = []
    for i in range(nran):
        rancat = os.path.join(path_data, ty+'_%d_full.ran.fits' %i)
        dens_datacat.append(size_datacat*2500/float(len(pf.open(rancat)[1].data)))
    DENS_datacat = np.mean(dens_datacat)
    DENS_datacat_std = np.std(dens_datacat)
    print('Density in SV3 fuji version 3 data: ',DENS_datacat,' +- ',DENS_datacat_std,' gals/deg2')


    #Check density input mock as N#target / area_sv3
    inputcat = os.path.join(path_input,'mockTargets_000_FirstGen_CutSky_alltracers_sv3bits.fits')
    temp = pf.open(inputcat)[1].data
    mask = (temp['SV3_DESI_TARGET']==d_bite[j])

    size_inputcat = float(len(temp[mask]))
    dens_inputcat = size_inputcat/area_sv3
    print('Density of input FirstGen MOCK in SV3 footprint as Ngal/SV3: ',dens_inputcat,' gals/deg2')

    #Check density full input mock, weighted by Nran_full/Nran_input

    fullmock = os.path.join(path_mock, ty+'_full.dat.fits')
    size_fullmock = float(len(pf.open(fullmock)[1].data))

    dens_fullmock = size_fullmock/area_sv3

    mydens = []
    for i in range(nran):
        ranfull = os.path.join(path_mock, ty+'_%d_full.ran.fits' %i)
        nranfull = float(len(pf.open(ranfull)[1].data))
        ranfull = os.path.join(path_input, 'mockRandom_%d_FirstGen_CutSky_alltracers_sv3bits.fits'%SNUMS[i])
        temp = pf.open(ranfull)[1].data
        mask = (temp['SV3_DESI_TARGET']==r_bit[j])
        nranin = float(len(temp[mask]))
        mydens.append(dens_fullmock*nranin/nranfull)
    MYDENS = np.mean(mydens)
    MYSTD = np.std(mydens)
    print('Density of full mock, weighted by Nran_input/Nran_full: ',MYDENS,' +- ',MYSTD,' gals/deg2')





