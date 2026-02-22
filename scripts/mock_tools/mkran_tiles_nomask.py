#make sure add the LSS repo to your python path
from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import fitsio
import numpy as np
import os, sys
import argparse
import random
import json
from desimodel.footprint import is_point_in_desi
from multiprocessing import Pool
import time
import LSS.common_tools as common


parser = argparse.ArgumentParser()
parser.add_argument("--survey", help="e.g., Y1, DA2",default='DA2')
parser.add_argument("--program", choices=['DARK','BRIGHT'],default='DARK')
parser.add_argument("--ranmin", help="minimum random index to process",default=1,type=int)
parser.add_argument("--ranmax", help="maximum-1 random index to process",default=18,type=int)

args = parser.parse_args()

tiletab = Table.read(f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/tiles-{args.program}.fits')

for i in range(args.ranmin,args.ranmax):
    ranf = f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/rands_intiles_{args.program}_nomask_'+str(i)+'.fits'

    print('making '+ranf)
    input_ran = fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-allsky-1-0.fits',columns=['RA','DEC'])
    sel_tiles = is_point_in_desi(tiletab,input_ran['RA'],input_ran['DEC'])
    input_ran = input_ran[sel_tiles]
    common.write_LSS_scratchcp(input_ran,ranf)
    print('made '+ranf)
print('all done')
