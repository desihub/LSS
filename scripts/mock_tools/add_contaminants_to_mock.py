import astropy.io.fits as fits
from astropy.table import Table,vstack
import numpy as np
from numpy.random import Generator, PCG64

rng = Generator(PCG64())

for i in range(478,500):
    if i == '477':
        continue
    
    num = rng.integers(0, 21)

    type_ = 'QSO'

    filein = f'/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/seed0{i}/QSO/forFA0.fits'    
    #filein = 'GLAM_rea0_QSO.fits'

    contaminants = f'/global/cfs/projectdirs/desi/mocks/cai/contaminants/DA2/loa-v1/v2/QSO/noveto/contaminants_rea{num}.fits'


    dat_ = Table.read(filein)

    highz = dat_['RSDZ'] > 2.1

    normal = dat_[~highz]
    qso_highz = dat_[highz]


    datC = Table.read(contaminants)

    size_cont = len(datC)

    newtargid = (np.arange(1,size_cont+1)+1e8*2**22).astype(int)
    print(newtargid)
    print('size contaminats/size sample', size_cont/float(len(dat_)))
    print('size contaminats/size sample normal', size_cont/float(len(normal)))

    replace_idx = rng.choice(len(normal), size=size_cont, replace=False)


    normal['RA'][replace_idx] = datC['RA']
    normal['DEC'][replace_idx] = datC['DEC']
    normal['TARGETID'][replace_idx] = newtargid


    thecat = vstack([normal, qso_highz])

    thecat.write(filein.replace('.fits', '_withcontaminants.fits'), overwrite=True)
    '''
    type_ = 'ELG'
    num = rng.integers(0, 8)

    #filein = 'GLAM_rea0_ELG.fits'
    filein = f'/global/cfs/cdirs/desi/mocks/cai/holi/v5.0/seed0{i}/ELG/forFA0.fits' 
    contaminants = f'/global/cfs/projectdirs/desi/mocks/cai/contaminants/DA2/loa-v1/v2/ELGnotqso/noveto/contaminants_rea{num}.fits'

    dat_ = Table.read(filein)

    lop = dat_['DESI_TARGET'] & 2**5 == 2**5
    elg_vlo = dat_[~lop]
    elg_lop = dat_[lop]



    datC = Table.read(contaminants)
    size_cont = len(datC)

    newtargid = (np.arange(1,size_cont+1)+1e8*2**23).astype(int)
    print(newtargid)
    print('size contaminats/size sample', size_cont/float(len(dat_)))
    print('size contaminats/size sample lop', size_cont/float(len(elg_lop)))

    
    replace_idx = rng.choice(len(elg_lop), size=size_cont, replace=False)

    elg_lop['RA'][replace_idx] = datC['RA']
    elg_lop['DEC'][replace_idx] = datC['DEC']
    elg_lop['TARGETID'][replace_idx] = newtargid

    thecat = vstack([elg_lop, elg_vlo])

    thecat.write(filein.replace('.fits', '_withcontaminants.fits'), overwrite=True)
    '''
