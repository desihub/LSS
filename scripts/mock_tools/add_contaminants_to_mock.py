import astropy.io.fits as fits
from astropy.table import Table,vstack
import numpy as np
from numpy.random import Generator, PCG64

rng = Generator(PCG64())


type_ = 'QSO'

if type_ == 'QSO':
    filein = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/QSO_v5/z1.400/forFA0_Y3_noimagingmask_applied.fits'

    contaminants = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/QSO/z1.400/contaminants_forFA0_Y3_noimagingmask_applied.fits'


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

    thecat.write('/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/QSO_v5/z1.400/forFA0_Y3_noimagingmask_applied_withcontaminants.fits')

elif type_ == 'ELG':

    filein = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/ELG_v5/z0.950/forFA0_Y3_noimagingmask_applied.fits'
    contaminants = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/ELG_v5/z0.950/contaminants_forFA0_Y3_noimagingmask_applied.fits'

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

    thecat.write('/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/ELG_v5/z0.950/forFA0_Y3_noimagingmask_applied_withcontaminants.fits')

