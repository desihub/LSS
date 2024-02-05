import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import join,Table

import LSS.common_tools as common
from LSS.globals import main


parser = argparse.ArgumentParser()
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
args = parser.parse_args()


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'

lssmapdirout = indir+'/hpmaps/'



tps = [args.tracers]
if args.tracers == 'all':
    tps = ['QSO','LRG','ELG_LOPnotqso']#'BGS_BRIGHT'

for tp in tps:
    mainp = main(tp,args.verspec,survey=args.survey)
    df_cutdisk = fitsio.read(indir+tp+'_full_HPmapcut.dat.fits')
    df = fitsio.read(indir+tp+'_full.dat.fits')
    mapn = fitsio.read(lssmapdirout+tp+'_mapprops_healpix_nested_nside'+str(nside)+'_N.fits')
    maps = fitsio.read(lssmapdirout+tp+'_mapprops_healpix_nested_nside'+str(nside)+'_S.fits')
    mapcuts = mainp.mapcuts
    df_cut = common.apply_map_veto_arrays(df,mapn,maps,mapcuts)
    sel_idmatch = np.isin(df_cut['TARGETID'],df_cutdisk['TARGETID'])
    df_cutnomatch = df_cut[~sel_idmatch]
    df_comb = np.concatenate((df_cutdisk,df_cutnomatch))
    print(tp,len(df_comb),len(np.unique(df_comb['TARGETID']))) 
