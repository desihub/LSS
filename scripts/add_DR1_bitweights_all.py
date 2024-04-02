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
args = parser.parse_args()



lssdir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/'+args.cat_version+'/'

if args.prog == 'BRIGHT':
    sys.exit('needs to be written')
    #alltids = fitsio.read(lssdir+'BGS_ANY_full_noveto.dat.fits',columns=['TARGETID'])
    #alltids = np.unique(alltids['TARGETID'])


if args.prog == 'DARK':
    tpl = ['LRG','QSO','ELG_LOPnotqso']

bitf = fitsio.read(lssdir+args.prog+'_bitweights.fits')
fl = ['full_noveto','full','full_HPmapcut','clustering','NGC_clustering','SGC_clustering']
for tp in tpl:
    for ft in fl:
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

