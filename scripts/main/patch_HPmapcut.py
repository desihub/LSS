import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import join,Table,vstack

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

nside = 256

tps = [args.tracers]
if args.tracers == 'all':
    tps = ['QSO','LRG','ELG_LOPnotqso']#'BGS_BRIGHT'

for tp in tps:
    mainp = main(tp,args.verspec,survey=args.survey)
    df_cutdisk = Table(fitsio.read(indir+tp+'_full_HPmapcut.dat.fits'))
    df_cutdisk_cols = list(df_cutdisk.dtype.names)
    #if 'BITWEIGHTS' in df_cutdisk_cols:
    #    df_cutdisk.remove_columns(['BITWEIGHTS','PROB_OBS'])
    #if 'BITWEIGHTS_1' in df_cutdisk_cols:
    #    df_cutdisk.remove_columns(['BITWEIGHTS_1','PROB_OBS_1','BITWEIGHTS_2','PROB_OBS_2'])
    df_cutdisk_cols = list(df_cutdisk.dtype.names)
    df = Table(fitsio.read(indir+tp+'_full.dat.fits'))
    df_cols = list(df.dtype.names)
    rem_cols = []
    for name in df_cols:
        if name not in df_cutdisk_cols:
            print(name+' not in HPmapcut file')
            rem_cols.append(name)
    if len(rem_cols) > 0:
        df.remove_columns(rem_cols)
    for name in df_cutdisk_cols:
        if name not in df_cols:
            df[name] = np.ones(len(df))
            print(name+' added to not file as 1')
        
    mapn = fitsio.read(lssmapdirout+tp.replace('-21.5','')+'_mapprops_healpix_nested_nside'+str(nside)+'_N.fits')
    maps = fitsio.read(lssmapdirout+tp.replace('-21.5','')+'_mapprops_healpix_nested_nside'+str(nside)+'_S.fits')
    mapcuts = mainp.mapcuts
    df_cut = common.apply_map_veto_arrays(df,mapn,maps,mapcuts)
    sel_idmatch = np.isin(df_cut['TARGETID'],df_cutdisk['TARGETID'])
    df_cutnomatch = df_cut[~sel_idmatch]
    #df_comb = np.concatenate((df_cutdisk,df_cutnomatch))
    print(df_cutdisk.dtype.names,df_cutnomatch.dtype.names)
    df_comb = vstack((df_cutdisk,df_cutnomatch))
    print(tp,len(df_comb),len(np.unique(df_comb['TARGETID']))) 
    #if tp[:3] != 'BGS':
    #    bitf = fitsio.read(mainp.darkbitweightfile)
    #else:
    #    bitf = fitsio.read(mainp.brightbitweightfile)
    #df_comb = join(df_comb,bitf,keys=['TARGETID'],join_type='left')
    print(len(df_comb['TARGETID']))
    print(df_comb.dtype.names)
    common.write_LSS(df_comb,indir+tp+'_full_HPmapcut.dat.fits')
