from astropy.table import Table,vstack
import numpy as np
from numpy.random import Generator, PCG64
import os
import argparse
import sys

rng = Generator(PCG64())


parser = argparse.ArgumentParser()
parser.add_argument("--tracer", default=None)
parser.add_argument("--realization", default=None)
parser.add_argument("--mock_version", default=None)
parser.add_argument("--mock", default='holi')
parser.add_argument("--file_name", default='forFA0.fits')

args = parser.parse_args()

if args.mock == 'holi':
    realization = str(args.realization).zfill(4)
    filein = os.path.join("/global/cfs/cdirs/desi/mocks/cai", args.mock, args.mock_version, "seed%s" % realization, args.tracer, args.file_name)
elif args.mock == 'abacus':

    realization = str(args.realization).zfill(3)
    filein = os.path.join("/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph%s/CutSky" % realization, args.tracer, args.file_name)
    #/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph024/CutSky/QSO/imforFA0_Y3_noimagingmask_applied.fits

elif args.mock == 'uchuuref':
    realization = str(args.realization).zfill(4)
    filein = os.path.join("/global/cfs/cdirs/desi/mocks/cai/Uchuu-SHAM/Y3-v2.0/0000/prep_altmtl", args.tracer, args.file_name)


#realization = str(args.realization).zfill(4)
#filein = os.path.join("/global/cfs/cdirs/desi/mocks/cai", args.mock, args.mock_version, "seed%s" % realization, args.tracer, args.file_name)

fileout = filein.replace('.fits', '_withcontaminants.fits')

if not os.path.isfile(filein):
    print('input file', filein, 'dont exist')
    sys.exit()

if os.path.isfile(fileout):
    print('output file', fileout, 'already exist')
    sys.exit()

dat_ = Table.read(filein)

if args.tracer == 'QSO':
    
    num = rng.integers(0, 99)

    contaminants = f'/global/cfs/projectdirs/desi/mocks/cai/contaminants/DA2/loa-v1/v2/QSO/noveto/contaminants_rea{num}.fits'


    highz = dat_['RSDZ'] > 2.1

    normal = dat_[~highz]
    qso_highz = dat_[highz]

    datC = Table.read(contaminants)

    size_cont = len(datC)

    newtargid = (np.arange(1,size_cont+1)+1e8*2**22).astype(int)
    
    print('size contaminants/size sample', size_cont/float(len(dat_)))

    replace_idx = rng.choice(len(normal), size=size_cont, replace=False)

    normal['RA'][replace_idx] = datC['RA']
    normal['DEC'][replace_idx] = datC['DEC']
    normal['TARGETID'][replace_idx] = newtargid

    thecat = vstack([normal, qso_highz])


if args.tracer == 'ELG':
    
    num = rng.integers(0, 33)

    contaminants = f'/global/cfs/projectdirs/desi/mocks/cai/contaminants/DA2/loa-v1/v2/ELGnotqso/noveto/contaminants_rea{num}.fits'

    lop = dat_['DESI_TARGET'] & 2**5 == 2**5
    elg_vlo = dat_[~lop]
    elg_lop = dat_[lop]

    datC = Table.read(contaminants)
    size_cont = len(datC)

    newtargid = (np.arange(1,size_cont+1)+1e8*2**23).astype(int)

    print('size contaminants/size sample', size_cont/float(len(dat_)))

    
    replace_idx = rng.choice(len(elg_lop), size=size_cont, replace=False)

    elg_lop['RA'][replace_idx] = datC['RA']
    elg_lop['DEC'][replace_idx] = datC['DEC']
    elg_lop['TARGETID'][replace_idx] = newtargid

    thecat = vstack([elg_lop, elg_vlo])

print('saving file', fileout)
thecat.write(fileout, overwrite=True)
