import astropy.io.fits as pf
import os
import numpy as np

types_ = ['LRG','ELG','QSO']
d_bite = [1,34,4]
r_bit = [1,2,4]


path_mock = '/global/cscratch1/sd/acarnero/SV3/LSS_rea000_univ1/fuji/LSScats/test'
path_data = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/fuji/LSScats/3'
path_input = '/global/cscratch1/sd/acarnero/SV3'

area_sv3 = 207.5

n_ran = 18
SNUMS = np.linspace(100, 5000, num=50, dtype=np.int)



for j,ty in enumerate(types_):
    print('Numbers for ',ty)
    #Check density data as Ngal/Nran*2500 in the full NOVETO catalog
    datacat = os.path.join(path_data, ty+'_full_noveto.dat.fits')
    size_datacat = float(len(pf.open(datacat)[1].data))
    dens_datacat = []
    nran = []
    for i in range(n_ran):
        rancat = os.path.join(path_data, ty+'_%d_full_noveto.ran.fits' %i)
        dens_datacat.append(size_datacat*2500/float(len(pf.open(rancat)[1].data)))
        nran.append(float(len(pf.open(rancat)[1].data)))
    DENS_datacat = np.mean(dens_datacat)
    DENS_datacat_std = np.std(dens_datacat)
    print('Density in DATA FULL NOVETO: ',"{0:.1f}".format(DENS_datacat),' +- ',"{0:.1f}".format(DENS_datacat_std),' gals/deg2')
    print('Ngal full noveto = ',"{0:.1f}".format(size_datacat))
    print('Nran full noveto = ',"{0:.1f}".format(np.mean(nran)),' +- ',"{0:.1f}".format(np.std(nran)))
    print('Ngal/Nran =',"{0:.1f}".format(size_datacat/np.mean(nran)))
    print('--------------------------')
    #Check density data as Ngal/Nran*2500 in the full catalog
    datacat = os.path.join(path_data, ty+'_full.dat.fits')
    size_datacat = float(len(pf.open(datacat)[1].data))
    dens_datacat = []
    nran = []
    for i in range(n_ran):
        rancat = os.path.join(path_data, ty+'_%d_full.ran.fits' %i)
        dens_datacat.append(size_datacat*2500/float(len(pf.open(rancat)[1].data)))
        nran.append(float(len(pf.open(rancat)[1].data)))
    DENS_datacat = np.mean(dens_datacat)
    DENS_datacat_std = np.std(dens_datacat)
    print('Density in DATA FULL: ',"{0:.1f}".format(DENS_datacat),' +- ',"{0:.1f}".format(DENS_datacat_std),' gals/deg2')
    print('Ngal full = ',"{0:.1f}".format(size_datacat))
    print('Nran full = ',"{0:.1f}".format(np.mean(nran)),' +- ',"{0:.1f}".format(np.std(nran)))
    print('Ngal/Nran =',"{0:.1f}".format(size_datacat/np.mean(nran)))
    print('--------------------------')
    


    #Check density input mock as N#target / area_sv3
    inputcat = os.path.join(path_input,'mockTargets_000_FirstGen_CutSky_alltracers_sv3bits.fits')
    temp = pf.open(inputcat)[1].data
    mask = (temp['SV3_DESI_TARGET']==d_bite[j])

    size_inputcat = float(len(temp[mask]))
    dens_inputcat = size_inputcat/area_sv3
    
    nran = []
    for i in range(n_ran):
        rancat = os.path.join(path_input, 'mockRandom_{ID}_FirstGen_CutSky_alltracers_sv3bits.fits'.format(ID=SNUMS[i]))
        temp = pf.open(rancat)[1].data
        mask = (temp['SV3_DESI_TARGET']==r_bit[j])
        nran.append(float(len(temp[mask])))

    print('Density input FirstGen MOCK: ',"{0:.1f}".format(dens_inputcat),' gals/deg2')
    print('Ngal Firsgen MOCK =',"{0:.1f}".format(float(len(temp[mask]))))
    print('Nran Firsgen MOCK = ',"{0:.1f}".format(np.mean(nran)),' +- ',"{0:.1f}".format(np.std(nran)))
    print('Ngal/Nran =',"{0:.1f}".format(float(len(temp[mask]))/np.mean(nran)))
    print('--------------------------')
    
    #Check density full noveto input mock, weighted by Nran_full/Nran_input

    fullmock = os.path.join(path_mock, ty+'_full_noveto.dat.fits')
    size_fullmock = float(len(pf.open(fullmock)[1].data))

    dens_fullmock = size_fullmock/area_sv3

    mydens = []
    nran = []
    for i in range(n_ran):
        ranfull = os.path.join(path_mock, ty+'_%d_full_noveto.ran.fits' %i)
        nranfull = float(len(pf.open(ranfull)[1].data))
        ranfull = os.path.join(path_input, 'mockRandom_%d_FirstGen_CutSky_alltracers_sv3bits.fits'%SNUMS[i])
        temp = pf.open(ranfull)[1].data
        mask = (temp['SV3_DESI_TARGET']==r_bit[j])
        nranin = float(len(temp[mask]))
        nran.append(nranfull)
        mydens.append(dens_fullmock*nranin/nranfull)
    MYDENS = np.mean(mydens)
    MYSTD = np.std(mydens)
    print('Density in MOCK FULL NOVETO, weighted by Nran_input/Nran_full: ',"{0:.1f}".format(MYDENS),' +- ',"{0:.1f}".format(MYSTD),' gals/deg2')
    print('Ngal FULL NOVETO MOCK =',"{0:.1f}".format(size_fullmock))
    print('Nran FULL NOVETO MOCK =',"{0:.1f}".format(np.mean(nran)),' +- ', "{0:.1f}".format(np.std(nran)))
    print('Ngal/Nran =',"{0:.1f}".format(size_fullmock/np.mean(nran)))
    print('--------------------------')


    #Check density full input mock, weighted by Nran_full/Nran_input

    fullmock = os.path.join(path_mock, ty+'_full.dat.fits')
    size_fullmock = float(len(pf.open(fullmock)[1].data))

    dens_fullmock = size_fullmock/area_sv3

    mydens = []
    nran = []
    for i in range(n_ran):
        ranfull = os.path.join(path_mock, ty+'_%d_full.ran.fits' %i)
        nranfull = float(len(pf.open(ranfull)[1].data))
        ranfull = os.path.join(path_input, 'mockRandom_%d_FirstGen_CutSky_alltracers_sv3bits.fits'%SNUMS[i])
        temp = pf.open(ranfull)[1].data
        mask = (temp['SV3_DESI_TARGET']==r_bit[j])
        nranin = float(len(temp[mask]))
        nran.append(nranfull)
        mydens.append(dens_fullmock*nranin/nranfull)
    MYDENS = np.mean(mydens)
    MYSTD = np.std(mydens)
    print('Density in MOCK FULL, weighted by Nran_input/Nran_full: ',"{0:.1f}".format(MYDENS),' +- ',"{0:.1f}".format(MYSTD),' gals/deg2')
    print('Ngal FULL MOCK =',"{0:.1f}".format(size_fullmock))
    print('Nran FULL MOCK =',"{0:.1f}".format(np.mean(nran)),' +- ',"{0:.1f}".format(np.std(nran)))
    print('Ngal/Nran =',"{0:.1f}".format(size_fullmock/np.mean(nran)))

    print('*************************')
    print('*************************')
    print('*************************')



