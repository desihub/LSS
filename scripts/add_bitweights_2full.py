import numpy as np
import fitsio
import glob
from astropy.table import Table,join

import LSS.common_tools as common
from LSS.globals import main
from LSS.bitweights import pack_bitweights

import logging
logname = 'add_bitweights'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--prog", choices=['DARK','BRIGHT'],default='DARK')
parser.add_argument("--cat_version",default='test')
parser.add_argument("--tracers",default='all')
parser.add_argument("--survey",default='DA2')
parser.add_argument("--specrel",default='kibo-v1')
parser.add_argument("--replace",default='y')
args = parser.parse_args()



lssdir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'+args.specrel+'/LSScats/'+args.cat_version+'/'


if args.prog == 'DARK':
    if args.tracers == 'all':
        tpl = ['LRG','QSO','ELG_LOPnotqso']
    else:
        tpl = [args.tracers]

if args.prog == 'BRIGHT':
    if args.tracers == 'all':
        tpl = ['BGS_BRIGHT','BGS_ANY']
    else:
        tpl = [args.tracers]


bitf = fitsio.read(lssdir+args.prog+'_bitweights.fits')
fl = ['full_noveto','full','full_HPmapcut']
for tp in tpl:
    fll = fl
    for ft in fll:
        inflnm = lssdir+tp+'_'+ft+'.dat.fits'
        infl = Table(fitsio.read(inflnm.replace('global','dvs_ro')))
        cols = list(infl.dtype.names)
        dojoin = 'y'
        if 'PROB_OBS' in cols:
            if args.replace == 'y':
                common.printlog('removing columns before adding info back',logger)
                infl.remove_columns(['PROB_OBS','BITWEIGHTS'])
            else:
                dojoin = 'n'
                common.printlog('PROB_OBS is in original and replace is set to n, so just moving to next file',logger)
        if dojoin == 'y':
            li = len(infl)
            infl = join(infl,bitf,keys=['TARGETID'],join_type='left')
            lij = len(infl)
            if li == lij:
                common.write_LSS_scratchcp(infl,inflnm,logger=logger)
            else:
                common.printlog('mismatch after join! '+tp,logger)
                print(tp,li,lij)    

