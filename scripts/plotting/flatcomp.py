import matplotlib.pyplot as plt
import numpy as np
import os
import sys

import fitsio
from astropy.table import join,Table


outdir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/plots/'
qt = 'COMP_TILE'

tps = ['ELG','ELG_LOP','LRG','QSO','BGS_ANY','BGS_BRIGHT']
for tp in tps:
    tars = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+tp+'targetsDR9v1.1.1.fits',columns=['TARGETID','RA','DEC'])
    dat = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'zdone_full_noveto.dat.fits',columns=['TARGETID',qt])
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
    
    plt.scatter(datno['RA'],np.sin(datno['DEC']*np.pi/180.),c=datno[qt],edgecolor='none',s=.1,vmax=vx,vmin=vm,zorder=1)
    plt.scatter(dato['RA'],np.sin(dato['DEC']*np.pi/180.),c=dato[qt],edgecolor='none',s=.1,vmax=vx,vmin=vm,zorder=1000)
    plt.ylim(-0.5,1)

    plt.colorbar()
    plt.title(tp+' '+qt)
    plt.xlabel('RA')
    plt.ylabel('sin(DEC)')
    plt.savefig(outdir+tp+'_flat'+qt+'.png')
    plt.clf()
    del dat
    del datno
    del dato