#standard python
import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import glob
import argparse
import gc
#gc.enable()
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
import healpy as hp

#import tracemalloc

#tracemalloc.start()

from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desitarget import targetmask
from desitarget.internal import sharedmem
from desimodel.footprint import is_point_in_desi
import desimodel.footprint as foot

#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.common_tools as common
import LSS.mkCat_singletile.fa4lsscat as fa
from LSS.globals import main

if os.environ['NERSC_HOST'] == 'cori':
    scratch = '$CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = '$PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 

import logging

# create logger
logname = 'LSSran'
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


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=scratch)
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='daily')
parser.add_argument("--ranmtl", help="make a random mtl file for the tile",default='n')
parser.add_argument("--rfa", help="run randoms through fiberassign",default='n')
parser.add_argument("--combhp", help="combine the random tiles together but in separate  healpix",default='n')
parser.add_argument("--combr", help="combine the random healpix files together",default='n')
parser.add_argument("--fullr", help="make the random files with full info, divided into healpix",default='n')
parser.add_argument("--refullr", help="make the full files from scratch rather than only updating pixels with new tiles",default='n')
parser.add_argument("--combfull", help="combine the full files in healpix into one file",default='n')
parser.add_argument("--clus", help="make the data/random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--maskz", help="apply sky line mask to redshifts?",default='n')
parser.add_argument("--faver", help="version of fiberassign code to use for random; versions for main should be 5.0.0 or greater",default='5.0.0')
parser.add_argument("--par", help="run different random number in parallel?",default='y')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=18,type=int) 
parser.add_argument("--ran_ind",help='index for the input randoms, just 1 by default',default=1,type=int)
parser.add_argument("--redos",help="whether or not to redo match to spec data (e.g., to add in a new column)",default='n')
parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')

args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
version = args.version
faver = args.faver
specrel = args.verspec
rm = int(args.minr)
rx = int(args.maxr)
par = True
if args.par == 'n':
    par = False


mkranmtl = False
if args.ranmtl == 'y':
    mkranmtl = True
runrfa = True#run randoms through fiberassign
if args.rfa == 'n':
    runrfa = False
combr = True
if args.combr == 'n':
    combr = False   
combhp = True
if args.combhp == 'n':
    combhp = False    

redos = False
if args.redos == 'y':
    redos = True

mkfullr = True #make the random files associated with the full data files
if args.fullr == 'n':
    mkfullr = False
mkclus = True #make the data/random clustering files; these are cut to a small subset of columns
mkclusran = True
if args.clus == 'n':
    mkclus = False
    mkclusran = False

if type == 'bright' or type == 'dark':
    #don't do any of the actual catalogs in this case
    mkclus = False
    mkclusran = False
    mkfullr = False

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'


if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    pr = 'BRIGHT'
    pdir = 'bright'
else:
    pr = 'DARK'
    pdir = 'dark'

pd = pdir

globtype = args.type
if 'dark' in args.type:
    globtype = 'LRG'
if 'bright' in args.type:
    print('changing globtype to BGS')
    globtype == 'BGS'
print(globtype)
mainp = main(globtype,args.verspec)

mt = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied


wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == pdir
if specrel != 'daily':
    if specrel == 'everest':
        specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+pdir+'-cumulative.fits')
        wd &= np.isin(mt['TILEID'],np.unique(specf['TILEID']))
mtld = mt[wd]
#print('found '+str(len(mtd))+' '+prog+' time main survey tiles that are greater than 85% of goaltime')
print('found '+str(len(mtld))+' '+pdir+' time main survey tiles with zdone true for '+specrel+' version of reduced spectra')

selt = np.isin(tiles['TILEID'],mtld['TILEID'])
ta = Table()
ta['TILEID'] = tiles[selt]['TILEID']
ta['RA'] = tiles[selt]['RA']
ta['DEC'] =tiles[selt]['DEC']

#tiles4comb = Table()
#tiles4comb['TILEID'] = mtld['TILEID']
#tiles4comb['ZDATE'] = mtld['LASTNIGHT']

#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/main/LSS/'




if not os.path.exists(maindir+'/logs'):
    os.mkdir(maindir+'/logs')
    print('made '+maindir+'/logs')



randir = maindir+'random'
for ii in range(args.minr,args.maxr):
    #logf.write('using random files '+str(rm)+ ' through '+str(rx)+' (this is python, so max is not inclusive)\n')
    if not os.path.exists(randir+str(ii)):
        os.mkdir(randir+str(ii))
        print('made '+randir+str(ii)+' random directory')
    if not os.path.exists(randir+str(ii)+'/healpix'):
        os.mkdir(randir+str(ii)+'/healpix')
        print('made '+randir+str(ii)+'/healpix'+' random directory')
        

ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    print('made '+ldirspec)


dirout = ldirspec+'LSScats/'+version+'/'
if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)


#construct a table with the needed tile information
# if len(mtld) > 0:
#     tilel = []
#     ral = []
#     decl = []
#     mtlt = []
#     fal = []
#     obsl = []
#     pl = []
#     fver = []
#     fahal = []
#     
#     #for tile,pro in zip(mtld['TILEID'],mtld['PROGRAM']):
#     for tile in mtld['TILEID']:
#         ts = str(tile).zfill(6)
#         try:
#             fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
#             tilel.append(tile)
#             ral.append(fht['TILERA'])
#             decl.append(fht['TILEDEC'])
#             mtlt.append(fht['MTLTIME'])
#             fal.append(fht['FA_RUN'])
#             obsl.append(fht['OBSCON'])
#             fav = fht['FA_VER']
#             try:
#                 if int(fav[:1]) >= 5:
#                     fav = '5.0.0'
#                 #else:
#                 #    print(fav)    
#             except:
#                 print(fav)        
#             pl.append(pr)
#         except:
#             print('failed to find and/or get info for tile '+ts)    
#     ta = Table()
#     ta['TILEID'] = tilel
#     ta['RA'] = ral
#     ta['DEC'] = decl
#     ta['MTLTIME'] = mtlt
#     ta['FA_RUN'] = fal
#     ta['OBSCON'] = obsl
#     ta['PROGRAM'] = pl
#     #ta['FA_HA'] = fahal
#     #ta['FA_VER'] = fver
#     #print(np.unique(fver,return_counts=True))
#     #wfv = (np.array(fver) == faver)
#     #mtld =  mtld[wfv]
#     #ta = ta[wfv]
# else:
#     print('no done tiles in the MTL')

print(len(ta))

print(specrel)


#print(tracemalloc.get_traced_memory())

hpxs = foot.tiles2pix(8, tiles=ta)

if combhp or mkfullr:
    
    if specrel == 'daily':
        specfo = ldirspec+'datcomb_'+pdir+'_spec_zdone.fits'
        specf = Table.read(specfo)
        specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    

    if specrel == 'everest':    

        #specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+type+'-cumulative.fits')
        #wt = np.isin(mtld['TILEID'],specf['TILEID'])
        #above two lines already done above
        specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+pdir+'-cumulative.fits')
        wt = np.isin(specf['TILEID'],mtld['TILEID']) #cut spec file to dark or bright time tiles
        specf = specf[wt]
        print('number of TILEID in spec data being used:')
        print(len(np.unique(specf['TILEID'])))
        specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    print('loaded specf file '+specfo)
    specfc = ct.cut_specdat(specf)
    if combhp:
        ntls = len(np.unique(specf['TILEID']))
        if ntls != len(ta):
            print(ntls,len(ta))
            sys.exit('mismatch in number of tileids NOT PROCEEDING')
    gtl = np.unique(specfc['TILELOCID'])
    del specfc

if type != 'dark' and type != 'bright' and mkfullr:
    if type == 'BGS_BRIGHT':
        bit = targetmask.bgs_mask[type]
        desitarg='BGS_TARGET'
    else:
        bit = targetmask.desi_mask[type]    
        desitarg='DESI_TARGET'
    del specf
    print('loading '+ldirspec+'datcomb_'+type+notqso+'_tarspecwdup_zdone.fits')
    specf = fitsio.read(ldirspec+'datcomb_'+type+notqso+'_tarspecwdup_zdone.fits')#,columns=['TARGETID','ZWARN','TILELOCID'])
    
    wg = np.isin(specf['TILELOCID'],gtl)
    specf = Table(specf[wg])
    print('length after selecting type and good hardware '+str(len(specf)))
    lznp = common.find_znotposs(specf)
    del specf
    print('finished finding znotposs')




#print(tracemalloc.get_traced_memory())
# tsnrcol = 'TSNR2_ELG'
# tsnrcut = 0
# if type[:3] == 'ELG':
#     #dchi2 = 0.9 #This is actually the OII cut criteria for ELGs
#     tsnrcut = 80
# if type == 'LRG':
#     #dchi2 = 16  
#     tsnrcut = 80          
# if type[:3] == 'BGS':
#     tsnrcol = 'TSNR2_BGS'
#     dchi2 = 40
#     tsnrcut = 1000

tsnrcut = mainp.tsnrcut
dchi2 = mainp.dchi2
tsnrcol = mainp.tsnrcol        


def doran(ii):
    ran_out = (args.ran_ind-1)*20+ii
    common.printlog('ran_out is '+str(ran_out),logger)

    #dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
    #dirrt = '/global/cscratch1/sd/adamyers/forashley/dr9/2.3.0.dev5334/randoms/resolve/'  
    dirrt = '/global/cfs/cdirs/desi/target/catalogs/dr9/2.4.0/randoms/resolve/'

    if mkranmtl:
        #print('making random mtl files for each tile')
        #ct.randomtiles_allmain_pix(ta,imin=ii,imax=ii+1,dirrt=dirrt+'randoms-1-'+str(ii))
        
        ct.randomtiles_allmain_pix_2step(ta,ii=ran_out,dirrt=dirrt+'randoms-'+str(args.ran_ind)+'-'+str(ii),logger=logger)


    
    if runrfa:
        print('DID YOU DELETE THE OLD FILES!!!')
        nd = 0
        for it in range(0,len(ta)):
            #print(it,len(mtld))    
            tile = ta['TILEID'][it]
            ts = str(tile).zfill(6)
            testfbaf = randir+str(ii)+'/fba-'+str(tile).zfill(6)+'.fits'
            if os.path.isfile(testfbaf):
                #print('fba file already made')
                pass
            else:
                fbah = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
                dt = fbah['RUNDATE'][:19]
                fav = fbah['FA_VER']
                if np.isin(fav,['2.2.0.dev2811','2.3.0','2.3.0.dev2838']):#2.3.0 confirmed to work for these
                    fav = '2.3.0'
                try:
                    if int(fav[:1]) >= 5:
                        fav = '5.0.0'
                except:
                    print(fav)        

                if fav == faver:
                    ttemp = Table(ta[it])
                    ttemp['OBSCONDITIONS'] = 516
                    ttemp['IN_DESI'] = 1
                    ttemp['MTLTIME'] = fbah['MTLTIME']
                    ttemp['FA_RUN'] = fbah['FA_RUN']
                    ttemp['PROGRAM'] = pr
                    try:
                        ttemp['FA_PLAN'] = fbah['FA_PLAN']
                        ttemp['FA_HA'] = fbah['FA_HA']
                        ttemp['FIELDROT'] = fbah['FIELDROT']
                    except:
                        print('did not add FA_PLAN and FIELDROT')
                #for i in range(rm,rx):
                    ttemp.write('tiletemp'+str(ii)+'.fits',format='fits', overwrite=True)
                    fa.getfatiles(randir+str(ii)+'/tilenofa-'+str(tile)+'.fits','tiletemp'+str(ii)+'.fits',dirout=randir+str(ii)+'/',dt = dt,faver=faver)
                    nd += 1
                    print('completed '+str(nd))
                    #del ttemp
                    #del fbah
                    #gc.collect()
                else:                   
                    #print(ttemp)
                    print('mismatch in fiberassign version, not doing fiberassign for '+str(tile)+' version is ' +fav)

    if combhp:
        if type == 'dark' or type == 'bright':
            
            npx = 0
            kc = ['ZWARN','LOCATION','FIBER','COADD_FIBERSTATUS','TILEID','TILELOCID','FIBERASSIGN_X','FIBERASSIGN_Y','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT'\
            ,'MEAN_DELTA_X','MEAN_DELTA_Y','RMS_DELTA_X','RMS_DELTA_Y','MEAN_PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B'\
            ,'TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
            'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
            'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','PRIORITY']
            

            for px in hpxs:
                print('combining target data for pixel '+str(px)+' '+str(npx)+' out of '+str(len(hpxs)))
                new = ct.combran_wdup_hp(px,ta,ii,randir,type,ldirspec,specf,keepcols=kc,redos=redos)
                if new:
                    tc = ct.count_tiles_better_px('ran',type,gtl,ii,specrel=specrel,px=px)
                    tc.write(ldirspec+'/healpix/rancomb_'+str(ii)+type+'_'+str(px)+'__Alltilelocinfo.fits',format='fits', overwrite=True)
                npx += 1
           
  

    if combr:
        print(len(mtld['TILEID']))
        #ct.combran(mtld,ii,randir,dirout,type,sv3_targetmask.desi_mask)
        if type == 'dark' or type == 'bright':
            
            s = 0
            npx =0 
            for px in hpxs:                
                tarfo = ldirspec+'healpix/rancomb_'+str(ii)+type+'_'+str(px)+'_wdupspec_zdone.fits'
                if os.path.isfile(tarfo):
                    tarf = Table.read(tarfo)
                    if s == 0:
                        tarfn = tarf
                        s = 1
                    else:
                        tarfn = vstack([tarfn,tarf],metadata_conflicts='silent')
                    print(len(tarfn),npx,len(hpxs))
                else:
                    print('file '+tarfo+' not found')
                npx += 1    
            tarfn.write(ldirspec+'rancomb_'+str(ii)+type+'wdupspec_zdone.fits',format='fits', overwrite=True)            
            tc = ct.count_tiles_better('ran',type,ii,specrel=specrel)
            tc.write(ldirspec+'/rancomb_'+str(ii)+type+'_Alltilelocinfo.fits',format='fits', overwrite=True)

        
    if mkfullr:
        maxp = 3400
        if type[:3] == 'LRG' or notqso == 'notqso':
            maxp = 3200
        if type[:3] == 'BGS':
            maxp = 2100

        npx = 0
        uhpxs = []
        if args.refullr == 'y':
            uhpxs = hpxs
        else:
            cf = dirout+type+notqso+'_'+str(ii)+'_full_noveto.ran.fits'
            dosel = False
            try:
                tls = fitsio.read(cf,columns=['TILEID'])
                dosel = True
            except:
                print('problem reading '+cf+' redoing all')
                uhpxs = hpxs
                
            if dosel:
                otls = np.unique(tls['TILEID'])
                print('got tileids currently in '+dirout+type+notqso+'_'+str(ii)+'_full_noveto.ran.fits')
                selt = ~np.isin(ta['TILEID'].astype(int),otls.astype(int))
                if len(ta[selt]) > 0:
                    uhpxs = foot.tiles2pix(8, tiles=ta[selt])
        for px in uhpxs:
            outf = ldirspec+'/healpix/'+type+notqso+'zdone_px'+str(px)+'_'+str(ii)+'_full.ran.fits'
            print(outf,npx,len(uhpxs))
            ct.mkfullran_px(ldirspec+'/healpix/',ii,imbits,outf,type,pdir,gtl,lznp,px,dirrt+'randoms-1-'+str(ii),maxp=maxp,min_tsnr2=tsnrcut)
            npx += 1  
        npx = 0
        s = 0
        #del lznp
        #del gtl
        

    if args.combfull == 'y':
        s = 0
        npx =0 
        outf = dirout+type+notqso+'_'+str(ii)+'_full_noveto.ran.fits'
        print('now combining to make '+outf)
        cols = ['GOODHARDLOC','ZPOSSLOC','PRIORITY','LOCATION', 'FIBER', 'TARGETID', 'RA', 'DEC', 'TILEID', 'ZWARN', 'FIBERASSIGN_X', 'FIBERASSIGN_Y', 'TSNR2_ELG_B', 'TSNR2_LYA_B', 'TSNR2_BGS_B', 'TSNR2_QSO_B', 'TSNR2_LRG_B', 'TSNR2_ELG_R', 'TSNR2_LYA_R', 'TSNR2_BGS_R', 'TSNR2_QSO_R', 'TSNR2_LRG_R', 'TSNR2_ELG_Z', 'TSNR2_LYA_Z', 'TSNR2_BGS_Z', 'TSNR2_QSO_Z', 'TSNR2_LRG_Z', 'TSNR2_ELG', 'TSNR2_LYA', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG', 'COADD_FIBERSTATUS', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'MEAN_DELTA_X', 'RMS_DELTA_X', 'MEAN_DELTA_Y', 'RMS_DELTA_Y', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'TILELOCID', 'NTILE', 'TILES','NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS']
        pl = []
        for px in hpxs:
            po = ldirspec+'/healpix/'+type+notqso+'zdone_px'+str(px)+'_'+str(ii)+'_full.ran.fits'
            if os.path.isfile(po):
                #pf = Table.read(po)
                try:
                    pf = fitsio.read(po,columns=cols)
                except:
                    print(po+' was corrupted')
                    ct.mkfullran_px(ldirspec+'/healpix/',ii,imbits,po,type,pdir,gtl,lznp,px,dirrt+'randoms-1-'+str(ii),maxp=maxp,min_tsnr2=tsnrcut)
                    pf = fitsio.read(po,columns=cols)
                pl.append(pf)
                print(npx,len(hpxs))
                #ptls = Table.read(po)
                #ptls.keep_columns(['TARGETID','TILES'])
                #if s == 0:
                #    pn = pf
                #    #ptlsn = ptls
                #    s = 1
                #else:
                    #pn = vstack([pn,pf],metadata_conflicts='silent')
                #    pn = np.hstack((pn,pf))
                    #ptlsn = vstack([ptlsn,ptls],metadata_conflicts='silent')
                #    print(len(pn),npx,len(hpxs))
            else:
                print('file '+po+' not found')
            npx += 1
        #pn = join(pn,ptlsn,keys=['TARGETID'],join_type='left')
        #pn.write(outf,overwrite=True,format='fits')
        print('stacking pixel arrays')
        #pn = vstack(pl,metadata_conflicts='silent')
        pn = np.hstack(pl)
        print('writing out')
        fitsio.write(outf,pn,clobber=True)
        del pn
    #logf.write('ran mkfullran\n')
    #print('ran mkfullran\n')


    if mkclusran:

        ct.mkclusran(dirout+type+notqso+'_',ii,zmask=zma,tsnrcut=tsnrcut,tsnrcol=tsnrcol)
    print('done with random '+str(ii))
    return True
        #ct.mkclusran(dirout+type+'Alltiles_',ii,zmask=zma)
    #logf.write('ran mkclusran\n')
    #print('ran mkclusran\n')
    
if __name__ == '__main__':
    if par:
        from multiprocessing import Pool
        import sys
        #N = int(sys.argv[2])
        #N = 32
        N = rx-rm+1
        #p = Pool(N)
        inds = []
        for i in range(rm,rx):
            inds.append(i)
        #with sharedmem.MapReduce() as pool:
        pool = sharedmem.MapReduce(np=N)
        with pool:
        
            def reduce( r):
                print('chunk done')
                return r
            pool.map(doran,inds,reduce=reduce)

        #p.map(doran,inds)
    else:
        for i in range(rm,rx):
            doran(i)