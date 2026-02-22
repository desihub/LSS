import numpy as np
import fitsio
import glob
from astropy.table import Table

import LSS.common_tools as common
from LSS.globals import main
from LSS.bitweights import pack_bitweights

import logging

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--prog", choices=['DARK','BRIGHT'])
parser.add_argument("--amtl_version",default='JL_Y1Run2')
parser.add_argument("--cat_version",default='test')
parser.add_argument("--nreal",default=128,type=int)
args = parser.parse_args()


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

logger.info('script is starting')
#just start with the mock 1 v3_1 altmtl as an example

lssdir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/'+args.cat_version+'/'

if args.prog == 'BRIGHT':
    alltids = fitsio.read(lssdir+'BGS_ANY_full_noveto.dat.fits',columns=['TARGETID'])
    alltids = np.unique(alltids['TARGETID'])


if args.prog == 'DARK':
    tpl = ['LRG','QSO','ELG_LOPnotqso']
    tidl = []
    for tp in tpl:
        f = fitsio.read(lssdir+tp+'_full_noveto.dat.fits',columns=['TARGETID'])
        tidl.append(f)
    alltids = np.concatenate(tidl)
    alltids = np.unique(alltids['TARGETID'])

logger.info(str(len(alltids))+ ' TARGETID will get their number of assignments tracked')

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

if args.prog == 'DARK':
    mainp = main('LRG','iron')
    pdir = 'dark'
else:
    mainp = main('BGS','iron')
    pdir = 'bright'

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



def get_good_real(dic,real_num):
    indir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/altmtl/'+args.amtl_version+args.prog+'/Univ'+str(real_num).zfill(3)+'/fa/MAIN'
    all_asgn = get_all_asgn(indir)
    asgn_tloc = 10000*all_asgn['TILEID'] +all_asgn['LOCATION']
    good_asgn = np.isin(asgn_tloc,gtl)
    good_tids = all_asgn['TARGETID'][good_asgn]
    asgn_real = np.isin(alltids,good_tids)
    assign_real_dic[real_num] = asgn_real
    del asgn_real

from multiprocessing import Pool
Nreal = args.nreal
inds = np.arange(0,Nreal)
#pool = sharedmem.MapReduce()
logger.info('about to get '+str(Nreal)+' realizations in parallel')
#with Pool() as pool:
#    #pool.map(get_good_real,inds)
#    pool.map(test,inds)

from multiprocessing import Process, Manager
manager = Manager()
assign_real_dic = manager.dict()
job = [Process(target=get_good_real, args=(assign_real_dic, i)) for i in inds]
_ = [p.start() for p in job]
_ = [p.join() for p in job]

logger.info('got all realizations')
#logger.info('dictionary keys are '+str(d.keys()))
import sys
#sys.exit()
logger.info('dictionary keys are '+str(assign_real_dic.keys()))
bool_list = []
for real in inds:
    bool_list.append(assign_real_dic[real])
bool_2d = np.vstack(bool_list).transpose()
logger.info('about to pack bitweights from array of shape '+str(np.shape(bool_2d)))
bitweights = pack_bitweights(bool_2d)

probl = np.zeros(len(alltids))
for real in inds:
    probl += assign_real_dic[real]*1.
probl = probl/Nreal   

outf = lssdir+args.prog+'_bitweights.fits'
out_tab = Table()
out_tab['TARGETID'] = alltids
out_tab['BITWEIGHTS'] = bitweights
out_tab['PROB_OBS'] = probl

common.write_LSS(out_tab,outf)

#h = np.histogram(probl)
#print(h) 
#for i in range(0,len(alltids)):
#    nt = 0
#    for real in inds:
#        nt += assign_real_dic[real][i]
#    prob = nt/Nreal
#    probl[i] = prob
#    if i%1e5 == 0:
#        logger.info(str(i)+' '+str(prob))