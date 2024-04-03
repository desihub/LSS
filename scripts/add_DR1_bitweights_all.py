import numpy as np
import fitsio
import glob
from astropy.table import Table,join

import LSS.common_tools as common
from LSS.globals import main
from LSS.bitweights import pack_bitweights

import logging

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--prog", choices=['DARK','BRIGHT'],default='DARK')
parser.add_argument("--cat_version",default='test')
parser.add_argument("--tracers",default='all')
args = parser.parse_args()



lssdir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/'+args.cat_version+'/'


if args.prog == 'DARK':
    if args.tracers == 'all':
        tpl = ['LRG','QSO','ELG_LOPnotqso']
    else:
        tpl = [args.tracers]

if args.prog == 'BRIGHT':
    if args.tracers == 'all':
        tpl = ['BGS_BRIGHT-21.5','BGS_BRIGHT','BGS_ANY']
    else:
        tpl = [args.tracers]


bitf = fitsio.read(lssdir+args.prog+'_bitweights.fits')
fl = ['full_noveto','full','full_HPmapcut','clustering','NGC_clustering','SGC_clustering']
for tp in tpl:
    fll = fl
    if tp == 'BGS_BRIGHT-21.5':
        fll = ['clustering','NGC_clustering','SGC_clustering'] #full catalogs don't really make sense for abs mag cut because they require a redshift
    for ft in fll:
        inflnm = lssdir+tp+'_'+ft+'.dat.fits'
        infl = fitsio.read(lssdir+tp+'_'+ft+'.dat.fits')
        li = len(infl)
        infl = join(infl,bitf,keys=['TARGETID'],join_type='left')
        lij = len(infl)
        if li == lij:
            common.write_LSS(infl,inflnm)
        else:
            print('mismatch after join!')
            print(tp,li,lij)    

