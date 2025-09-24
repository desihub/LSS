import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import join,Table

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected",default='all')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='main')
parser.add_argument("--verspec",help="version for redshifts",default='daily')

args = parser.parse_args()



outdir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'+args.verspec+'/LSScats/'+args.version+'/plots/'
qt = 'COMP_TILE'

if args.type == 'all':
    tps = ['LRG','BGS_ANY','BGS_BRIGHT','QSO','ELG','ELG_LOP']
else:
    tps = args.type.split(',')    
for tp in tps:
    tars = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+tp+'targetsDR9v1.1.1.fits',columns=['TARGETID','RA','DEC'])
    dat = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'+args.verspec+'/LSScats/test/'+tp+'_full_noveto.dat.fits',columns=['TARGETID',qt])
    dat = join(tars,dat,keys=['TARGETID'],join_type='left')
    mrows = dat['COMP_TILE'].mask
    dat['OBS'] = np.zeros(len(dat))
    dat['OBS'][~mrows] = 1
    dat['COMP_TILE'][mrows] = 0



    sel = dat['OBS'] == 1
    dato = dat[sel]
    datno = dat[~sel]
    if tp == 'ELG':
         rns = np.random.random_sample(len(datno))
         sel = rns < 0.1
         #sel |= dat['COMP_TILE'] > 0
         datno = datno[sel]
    if tp == 'BGS_ANY':
         rns = np.random.random_sample(len(datno))
         sel = rns < 0.2
         #sel |= dat['COMP_TILE'] > 0
         datno = datno[sel]
    if tp == 'LRG':
         rns = np.random.random_sample(len(datno))
         sel = rns < 0.5
         #sel |= dat['COMP_TILE'] > 0
         datno = datno[sel]



    print(len(dato),len(datno))

    
    vm = np.min(dat[qt])
    vx = np.max(dat[qt])
    ras = datno['RA']
    sel = ras > 300
    ras[sel] -= 360
    plt.scatter(ras,np.sin(datno['DEC']*np.pi/180.),c=datno[qt],edgecolor='none',s=.1,vmax=vx,vmin=vm,zorder=1)
    ras = dato['RA']
    sel = ras > 300
    ras[sel] -= 360    
    plt.scatter(ras,np.sin(dato['DEC']*np.pi/180.),c=dato[qt],edgecolor='none',s=.1,vmax=vx,vmin=vm,zorder=1000)
    plt.ylim(-0.5,1)

    plt.colorbar()
    plt.title(tp+' '+qt)
    plt.xlabel('RA')
    plt.ylabel('sin(DEC)')
    plt.savefig(outdir+tp+'_flat'+qt+'.png',dpi=1000)
    plt.clf()
    del dat
    del datno
    del dato
    print(tp+' done')