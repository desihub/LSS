from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import fitsio
import numpy as np
import os, sys
import argparse
import random
import json
from desitarget.targetmask import desi_mask, bgs_mask, obsconditions
from desimodel.footprint import is_point_in_desi
from multiprocessing import Pool
import time
import LSS.common_tools as common
from LSS.imaging import get_pixel_bitmasknobs as bitmask #get_nobsandmask
from LSS.main.cattools import count_tiles_better
from LSS.globals import main
tiletab = Table.read('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/tiles-DARK.fits')
#just using 1 random file for now
ranf = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/randoms/rands_intiles_DARK_with_imagingmask_{ID}.fits'
ranNOMASK= '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/randoms/rands_intiles_DARK_NO_imagingmask_{ID}.fits'
mainp = main(tp = 'LRG', specver = 'loa-v1')

for i in range(18):
    input_ran = fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-allsky-1-{ID}.fits'.format(ID=i),
                            columns=['RA','DEC','BRICKNAME','BRICKID','NOBS_R','NOBS_G','NOBS_Z','MASKBITS'])
    sel_tiles = is_point_in_desi(tiletab,input_ran['RA'],input_ran['DEC'])
    input_ran = input_ran[sel_tiles]
    common.write_LSS_scratchcp(input_ran, ranNOMASK.format(ID=i))
    continue
    print(len(input_ran))
    targets = common.cutphotmask(input_ran, bits=mainp.imbits)
    print(len(targets))
    common.write_LSS_scratchcp(targets, ranf.format(ID=i))

