import astropy.io.fits as fits
from astropy.table import Table,vstack
import numpy as np
from numpy.random import Generator, PCG64
import os

rng = Generator(PCG64())

type_ = 'QSO'


list_to_be_done = [201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 477, 617, 647]

if type_ == 'QSO':
    for i in list_to_be_done:
        realization = str(i).zfill(4)
        num = rng.integers(0, 99)


        filein = f'/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/seed{realization}/QSO/forFA0.fits'    
        
        if not os.path.isfile(filein):
            print('dont',i)
            continue

        if os.path.isfile(filein.replace('.fits', '_withcontaminants.fits')):
            print('already done', i)
            continue

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

type_ = 'ELG'

if type_ == 'ELG':
    for i in list_to_be_done:
        realization = str(i).zfill(4)
        num = rng.integers(0, 15)

        filein = f'/global/cfs/cdirs/desi/mocks/cai/holi/v5.0_Y5/seed{realization}/ELG/forFA0.fits' 

        contaminants = f'/global/cfs/projectdirs/desi/mocks/cai/contaminants/DA2/loa-v1/v2/ELGnotqso/noveto/contaminants_rea{num}.fits'

        if not os.path.isfile(filein):
            print('dont',i)
            continue

        if os.path.isfile(filein.replace('.fits', '_withcontaminants.fits')):
            print('already done', i)
            continue


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
