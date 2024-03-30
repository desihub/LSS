import numpy as np
import fitsio
import glob
from astropy.table import Table

import LSS.common_tools as common
from LSS.globals import main

import logging

# create logger
logname = 'bitweights'
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


#just start with the mock 1 v3_1 altmtl as an example

alltids = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/forFA1.fits',columns=['TARGETID'])
alltids = np.unique(alltids['TARGETID'])

def removeLeadingZeros(num): 
  
    # traverse the entire string 
    for i in range(len(num)): 
  
        # check for the first non-zero character 
        if num[i] != '0': 
            # return the remaining string 
            res = num[i::]; 
            return res; 
          
    # If the entire string is traversed 
    # that means it didn't have a single 
    # non-zero character, hence return "0" 
    return "0";

def get_all_asgn(indir):
    fls = glob.glob(indir+'/*/fba*.fits')
    assign_list = []
    for fl in fls:
        asgn = Table(fitsio.read(fl,columns=['FIBER', 'TARGETID', 'LOCATION']))
        sp = fl.split('-')
        tid = int(removeLeadingZeros(sp[-1].strip('.fits')))
        #print(tid)
        asgn['TILEID'] = tid
        sel = asgn['TARGETID'] > 0
        assign_list.append(asgn[sel])
    all_asgn = np.concatenate(assign_list)
    return all_asgn

#get the list of good tilelocid

mainp = main('LRG','iron')
pdir = 'dark'
mt = mainp.mtld
tiles = mainp.tiles


tsnrcut = mainp.tsnrcut
dchi2 = mainp.dchi2
tnsrcol = mainp.tsnrcol        
badfib = mainp.badfib


wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == pdir
wd &=mt['ZDATE'] < 20220900 #Y1 cutoff

mtld = mt[wd]
ldirspec = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/'
specfo = ldirspec+'datcomb_'+pdir+'_spec_zdone.fits'
logger.info('loading specf file '+specfo)
specf = Table(fitsio.read(specfo))
sel = np.isin(specf['TILEID'],mtld['TILEID'])
specf = specf[sel]
specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    
logger.info('loaded specf file '+specfo)
specfc = common.cut_specdat(specf,badfib=mainp.badfib,tsnr_min=tsnrcut,tsnr_col=tnsrcol)
gtl = np.unique(specfc['TILELOCID'])

assign_real_dic = {}

def get_good_real(real_num):
    indir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/altmtl1_R64/Univ'+str(real_num).zfill(3)+'/fa/MAIN'
    all_asgn = get_all_asgn(indir)
    asgn_tloc = 10000*all_asgn['TILEID'] +all_asgn['LOCATION']
    good_asgn = np.isin(asgn_tloc,gtl)
    good_tids = all_asgn['TARGETID'][gtl]
    asgn_real = np.isin(alltids,good_tids)
    assign_real_dic[real_num] = asgn_real
    del asgn_real

from multiprocessing import Pool
Nreal = 64
inds = np.arange(0,Nreal)
pool = sharedmem.MapReduce()
logger.info('about to get '+str(Nreal)+' realizations in parallel')
with pool:
    pool.map(get_good_real,inds)

logger.info('got all realizations')

probl = np.zeros(len(alltids))
for i in range(0,len(alltids)):
    nt = 0
    for real in inds:
        nt += assign_real_dic[real][i]
    prob = nt/Nreal
    probl[i] = prob
