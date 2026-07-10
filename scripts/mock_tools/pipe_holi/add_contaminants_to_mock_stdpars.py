#!/usr/bin/env python3

import astropy.io.fits as fits
from astropy.table import Table,vstack
import numpy as np
from numpy.random import Generator, PCG64
import os
import argparse

rng = Generator(PCG64())


def process_one_file(filein, out_f):
    if filein.find('QSO') > 0:
        num = rng.integers(0, 99)
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
        thecat.write(out_f, overwrite=True)

    elif filein.find('ELG') > 0:
        num = rng.integers(0, 15)

        contaminants = f'/global/cfs/projectdirs/desi/mocks/cai/contaminants/DA2/loa-v1/v2/ELGnotqso/noveto/contaminants_rea{num}.fits'
        dat_ = Table.read(filein)

        lop = dat_['DESI_TARGET'] & 2**5 == 2**5
        elg_vlo = dat_[~lop]
        elg_lop = dat_[lop]
        num = rng.integers(0, 15)

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

        thecat.write(out_f, overwrite=True)
        filein.replace('.fits', '_withcontaminants.fits')
    elif filein.find('LRG') > 0:
        os.system(f"touch {out_f}")


# Main
#
parser = argparse.ArgumentParser()

parser.add_argument("--inputs", nargs="+")
parser.add_argument("--outputs", nargs="+")
args = parser.parse_args()

print("adding contaminants to mock catalogs")
print(f"inputs: {args.inputs}")
print(f"outputs: {args.outputs}")

for i, o in zip(args.inputs, args.outputs):
    process_one_file(i, o)
