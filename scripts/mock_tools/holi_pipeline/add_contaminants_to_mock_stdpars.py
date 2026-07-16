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
        if len(normal) > 0:
            print('size contaminats/size sample normal', size_cont/float(len(normal)))
        else:
            print('size contaminats/size sample normal', 'inf (normal sample is empty)')
        n_replace = min(size_cont, len(normal))
        if n_replace < size_cont:
            print(f'Warning: truncating contaminants from {size_cont} to {n_replace} for QSO.')
        if n_replace > 0:
            replace_idx = rng.choice(len(normal), size=n_replace, replace=False)
            cont_idx = rng.choice(size_cont, size=n_replace, replace=False)
            normal['RA'][replace_idx] = datC['RA'][cont_idx]
            normal['DEC'][replace_idx] = datC['DEC'][cont_idx]
            normal['TARGETID'][replace_idx] = newtargid[cont_idx]
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
        if len(elg_lop) > 0:
            print('size contaminats/size sample lop', size_cont/float(len(elg_lop)))
        else:
            print('size contaminats/size sample lop', 'inf (elg_lop sample is empty)')

        n_replace = min(size_cont, len(elg_lop))
        if n_replace < size_cont:
            print(f'Warning: truncating contaminants from {size_cont} to {n_replace} for ELG.')
        if n_replace > 0:
            replace_idx = rng.choice(len(elg_lop), size=n_replace, replace=False)
            cont_idx = rng.choice(size_cont, size=n_replace, replace=False)

            elg_lop['RA'][replace_idx] = datC['RA'][cont_idx]
            elg_lop['DEC'][replace_idx] = datC['DEC'][cont_idx]
            elg_lop['TARGETID'][replace_idx] = newtargid[cont_idx]

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
