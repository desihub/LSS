import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import join,Table
import healpy as hp

from LSS.imaging import densvar

parser = argparse.ArgumentParser()
parser.add_argument("--version", help="catalog version",default='EDAbeta')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='SV3')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='fuji')
args = parser.parse_args()


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'+args.verspec+'/LSScats/'+args.version+'/'
zcol = 'Z'
nran = 18

tps = [args.tracers]
if args.tracers == 'all':
    tps = ['QSO','LRG','BGS_BRIGHT','ELG_LOPnotqso']

zdw = ''#'zdone'

regl = ['_N','_S']

if args.survey == 'SV3' and args.tracers == 'all':
    tps = ['QSO','LRG','BGS_ANY','ELGnotqso']

tot = 0 
for tp in tps:
    
    for nr in range(0,nran):
        rffh = fitsio.read_header(indir+tp+zdw+'_'+str(nr)+'_full.ran.fits',ext=1)   
        print(tp+' area is '+str(rffh['NAXIS2']/2500)+' deg2, using random '+str(nr))

    tot_tp = 0
    for reg in regl:
        dtf = fitsio.read_header(indir+tp+zdw+reg+'_clustering.dat.fits',ext=1)
        ncat = dtf['NAXIS2']
        print('number for '+tp+' in '+reg +' is '+str(ncat))
        tot_tp += ncat
    print('number for '+tp+' is '+str(tot_tp))
    tot += tot_tp
print('total number for '+args.survey +' is '+str(tot))