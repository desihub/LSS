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
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desitarget import targetmask
from desitarget.internal import sharedmem
from desimodel.footprint import is_point_in_desi

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


#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.common_tools as common
#import LSS.mkCat_singletile.fa4lsscat as fa
from LSS.globals import main

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    logger.info('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=scratch)
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='main')
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='daily')
parser.add_argument("--ranmtl", help="make a random mtl file for the tile",default='n')
parser.add_argument("--rfa", help="run randoms through fiberassign",default='n')
parser.add_argument("--combr", help="combine the random tiles together",default='n')
parser.add_argument("--combwspec", help="combine the random potential assignment info with spec info",default='n')
parser.add_argument("--counttiles", help="get NTILE, etc. counts",default='n')

parser.add_argument("--fullr", help="make the random files associated with the full data files",default='n')
parser.add_argument("--fullr_mode", help="if prog, noveto files are only split dark/bright",default='prog')
parser.add_argument("--mkdupranmasked",help="make duplicate randoms but with masks applied, to be used for randoms",default='n')
parser.add_argument("--hpmapcut", help="string indicating whether healpix map cut gets applied",default='_HPmapcut')
parser.add_argument("--add_veto", help="add veto column to the full files",default='n')
parser.add_argument("--fillran", help="add columns",default='n')
parser.add_argument("--apply_veto", help="apply vetos to the full files",default='n')
parser.add_argument("--add_tlcomp", help="add completeness FRAC_TLOBS_TILES to randoms",default='n')

parser.add_argument("--clus", help="make the data/random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--maskz", help="apply sky line mask to redshifts?",default='n')
parser.add_argument("--faver", help="version of fiberassign code to use for random; versions for main should be 5.0.0 or greater",default='5.0.0')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=4) 
parser.add_argument("--par", help="run different random number in parallel?",default='y')
parser.add_argument("--nproc", help="number to run in parallel",default=9)

parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
parser.add_argument("--newspec",help="if y, merge in redshift info even if no new tiles",default='n')

args = parser.parse_args()
logger.info(args)



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

#if specrel == 'daily':
#    sys.exit('support for daily needs to be written back in')
    
mkranmtl = False
if args.ranmtl == 'y':
    mkranmtl = True
runrfa = True#run randoms through fiberassign
if args.rfa == 'n':
    runrfa = False
combr = True
if args.combr == 'n':
    combr = False   


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
    #mkfullr = False

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




# tiles4comb = Table()
# tiles4comb['TILEID'] = mtld['TILEID']
# tiles4comb['ZDATE'] = mtld['LASTNIGHT']

#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/'+args.survey+'/LSS/'

randir = basedir +'/main/LSS/random'



if not os.path.exists(maindir+'/logs'):
    os.mkdir(maindir+'/logs')
    logger.info('made '+maindir+'/logs')

if not os.path.exists(maindir+'/LSScats'):
    os.mkdir(maindir+'/LSScats')
    logger.info('made '+maindir+'/LSScats')




#randir = maindir+'random'
#logf.write('using random files '+str(rm)+ ' through '+str(rx)+' (this is python, so max is not inclusive)\n')
#for i in range(rm,rx):
#    if not os.path.exists(maindir+'random'+str(i)):
#        os.mkdir(maindir+'random'+str(i))
#        logger.info('made '+str(i)+' random directory')

ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    logger.info('made '+ldirspec)


dirout = ldirspec+'LSScats/'+version+'/'
if not os.path.exists(dirout):
    os.makedirs(dirout)
    logger.info('made '+dirout)

mainp = main(type,args.verspec)

mt = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied


tsnrcut = mainp.tsnrcut
dchi2 = mainp.dchi2
tnsrcol = mainp.tsnrcol        
zmin = mainp.zmin
zmax = mainp.zmax
badfib = mainp.badfib


wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == pdir
if args.survey == 'Y1':
    wd &=mt['ZDATE'] < 20220900

if args.survey == 'DA2':
    wd &=mt['ZDATE'] < 20240410


mtld = mt[wd]
#logger.info('found '+str(len(mtd))+' '+prog+' time main survey tiles that are greater than 85% of goaltime')
logger.info('found '+str(len(mtld))+' '+pdir+' time main survey tiles with zdone true for '+specrel+' version of reduced spectra')

selt = np.isin(tiles['TILEID'],mtld['TILEID'])
ta = Table()
ta['TILEID'] = tiles[selt]['TILEID']
ta['RA'] = tiles[selt]['RA']
ta['DEC'] =tiles[selt]['DEC']



#if mkfullr or combr:
specfo = ldirspec+'datcomb_'+pdir+'_spec_zdone.fits'
logger.info('loading specf file '+specfo)
specf = Table(fitsio.read(specfo.replace('global','dvs_ro')))
sel = np.isin(specf['TILEID'],mtld['TILEID'])
specf = specf[sel]
specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    
logger.info('loaded specf file '+specfo)
#specfc = common.cut_specdat(specf,badfib=mainp.badfib,tsnr_min=tsnrcut,tsnr_col=tnsrcol,fibstatusbits=mainp.badfib_status)
specfc = common.cut_specdat(specf,badfib=mainp.badfib_td,tsnr_min=tsnrcut,tsnr_col=tnsrcol,fibstatusbits=mainp.badfib_status,remove_badfiber_spike_nz=True,mask_petal_nights=True,logger=logger)
gtl = np.unique(specfc['TILELOCID'])
del specfc

if mkfullr and args.fullr_mode != 'prog':
    logger.info('loading '+ldirspec+'datcomb_'+type+notqso+'_tarspecwdup_zdone.fits')
    specft = fitsio.read(ldirspec+'datcomb_'+type+notqso+'_tarspecwdup_zdone.fits')#,columns=['TARGETID','ZWARN','TILELOCID'])

    wg = np.isin(specft['TILELOCID'],gtl)
    specft = Table(specft[wg])
    logger.info('length after selecting type and good hardware '+str(len(specf)))

    lznp = common.find_znotposs(specft,logname=logname)    
    logger.info('finished finding znotposs')
    del specft

    if type == 'BGS_BRIGHT':
        bit = targetmask.bgs_mask[type]
        desitarg='BGS_TARGET'
    else:
        bit = targetmask.desi_mask[type]    
        desitarg='DESI_TARGET'


kc = ['LOCATION','FIBER','TILEID','TILELOCID','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','PRIORITY']
specf.keep_columns(kc)


logger.info(len(ta))
logger.info('done with preliminaries')

def doran(ii):
    logger.info('doing random '+str(ii))
    #dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
    #dirrt = '/global/cscratch1/sd/adamyers/forashley/dr9/2.3.0.dev5334/randoms/resolve/'  
    dirrt = '/global/cfs/cdirs/desi/target/catalogs/dr9/2.4.0/randoms/resolve/'

    if mkranmtl:
        logger.info('making random mtl files for each tile')
        #ct.randomtiles_allmain_pix(ta,imin=ii,imax=ii+1,dirrt=dirrt+'randoms-1-'+str(ii))
        ct.randomtiles_allmain_pix_2step(ta,ii=ii,dirrt=dirrt+'randoms-1-'+str(ii))


    
    if runrfa:
        logger.info('DID YOU DELETE THE OLD FILES!!!')
        for it in range(0,len(ta)):
            #logger.info(it,len(mtld))    
            tile = ta['TILEID'][it]
            ts = str(tile).zfill(6)
            fbah = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
            dt = fbah['RUNDATE'][:19]
            fav = fbah['FA_VER']
            if np.isin(fav,['2.2.0.dev2811','2.3.0','2.3.0.dev2838']):#2.3.0 confirmed to work for these
                fav = '2.3.0'
            try:
                if int(fav[:1]) >= 5:
                    fav = '5.0.0'
            except:
                logger.info(fav)        

            if fav == faver:
                ttemp = Table(ta[it])
                ttemp['OBSCONDITIONS'] = 516
                ttemp['IN_DESI'] = 1
                try:
                    ttemp['FA_PLAN'] = fbah['FA_PLAN']
                    ttemp['FA_HA'] = fbah['FA_HA']
                    ttemp['FIELDROT'] = fbah['FIELDROT']
                except:
                    logger.info('did not add FA_PLAN and FIELDROT')
                #for i in range(rm,rx):
                testfbaf = randir+str(ii)+'/fba-'+str(tile).zfill(6)+'.fits'
                if os.path.isfile(testfbaf):
                    logger.info('fba file already made')
                else:                   
                    logger.info(ttemp)
                    logger.info(fav,dt)
                    ttemp.write('tiletemp'+str(ii)+'.fits',format='fits', overwrite=True)
                    fa.getfatiles(randir+str(ii)+'/tilenofa-'+str(tile)+'.fits','tiletemp'+str(ii)+'.fits',dirout=randir+str(ii)+'/',dt = dt,faver=faver)

 

    if combr:
        logger.info(len(mtld['TILEID']))
        #ct.combran(mtld,ii,randir,dirout,type,sv3_targetmask.desi_mask)
        if type == 'dark' or type == 'bright':
            
            #kc = ['ZWARN','LOCATION','FIBER','COADD_FIBERSTATUS','TILEID','TILELOCID','FIBERASSIGN_X','FIBERASSIGN_Y','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT'\
            #,'MEAN_DELTA_X','MEAN_DELTA_Y','RMS_DELTA_X','RMS_DELTA_Y','MEAN_PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B'\
            #,'TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
            #'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
            #'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','PRIORITY']
            outf = maindir+'random'+str(ii)+'/rancomb_'+type+'wdup_Alltiles.fits'
            new = ct.combran_wdup(mtld,ii,randir,outf,keepcols=kc)
            #tiles,rann,randir,outf,keepcols=[]
            #ct.combran_wdup(mtld,ii,randir,type,ldirspec,specf,outf,keepcols=kc,collf=maindir+'random'+str(ii)+'collisions-'+pr+'.fits')
            if new or args.newspec == 'y':
                ct.combran_wdupspec(ii,type,ldirspec,specf,outf,keepcols=kc,mask_coll=True,collf=maindir+'random'+str(ii)+'collisions-'+pr+'.fits')
                tc = ct.count_tiles_better('ran',type,ii,specrel=specrel,survey=args.survey)
                tc.write(ldirspec+'/rancomb_'+str(ii)+type+'_Alltilelocinfo.fits',format='fits', overwrite=True)

    if args.combwspec == 'y':
        logger.info('combining with spec info')
        infile = maindir+'random'+str(ii)+'/pota-'+type.upper()+'.fits'
        mask_coll = False
        if args.survey != 'main':
            mask_coll = True
        ct.combran_wdupspec(ii,type,ldirspec,specf,infile,keepcols=kc,mask_coll=mask_coll,logger=logger)
    
    if args.counttiles == 'y':    
        logger.info('counting tiles')
        tc = ct.count_tiles_better('ran',type,ii,specrel=specrel,survey=args.survey,gtl=gtl)
        common.write_LSS_scratchcp(tc,ldirspec+'/rancomb_'+str(ii)+type+'_Alltilelocinfo.fits')
        #tc.write(ldirspec+'/rancomb_'+str(ii)+type+'_Alltilelocinfo.fits',format='fits', overwrite=True)

    if args.mkdupranmasked == 'y':
        outf = dirout+type+notqso+'_'+str(ii)+'_dupran_masked'+args.hpmapcut+'.fits'
        mapn = None
        maps = None
        mapcuts = None
        if args.hpmapcut == '_HPmapcut':
            lssmapdirout = dirout+'/hpmaps/'
            tracer_clus = type+notqso
            tracer_clushp = tracer_clus
            if tracer_clus == 'BGS_ANY':
                tracer_clushp = 'BGS_BRIGHT'
            if 'ELG' in tracer_clus:
                tracer_clushp = 'ELG_LOPnotqso'
            mapn = fitsio.read(lssmapdirout+tracer_clushp+'_mapprops_healpix_nested_nside256_N.fits')
            maps = fitsio.read(lssmapdirout+tracer_clushp+'_mapprops_healpix_nested_nside256_S.fits')
            mapcuts = mainp.mapcuts

        ct.mk_maskedran_wdup(gtl,ldirspec,ii,imbits,outf,pdir,ebits,notqso='',hpmapcut=args.hpmapcut,ftiles=None,mapn=mapn,maps=maps,mapcuts=mapcuts,reccircmasks=mainp.reccircmasks)        

    if mkfullr:
        maxp = 3400
        if type[:3] == 'LRG' or notqso == 'notqso':
            maxp = 3200
        if type[:3] == 'BGS':
            maxp = 2100
        if args.fullr_mode == 'prog':
            outf = dirout+pdir+'_'+str(ii)+'_full_noveto.ran.fits'
            logger.info('about to make full ran '+outf)
            ct.mkfullran_prog(gtl,ldirspec,ii,imbits,outf,pdir)
        
        else:
            outf = dirout+type+notqso+'_'+str(ii)+'_full_noveto.ran.fits'
            logger.info('about to make full ran '+outf)
            ct.mkfullran(gtl,lznp,ldirspec,ii,imbits,outf,type,pdir,notqso=notqso,maxp=maxp,min_tsnr2=tsnrcut)
        
    #logf.write('ran mkfullran\n')
    #logger.info('ran mkfullran\n')
    if args.add_veto == 'y':
        if args.fullr_mode == 'prog':
            fin = dirout+pdir+notqso+'_'+str(ii)+'_full_noveto.ran.fits'
        else:
            fin = dirout+type+notqso+'_'+str(ii)+'_full_noveto.ran.fits'
        logger.info('adding veto column to '+fin)
        common.add_veto_col(fin,ran=True,tracer_mask=type[:3].lower(),rann=ii,logger=logger)

    if args.fillran == 'y':
        
        if args.fullr_mode == 'prog':
            fn = dirout+pdir+'_'+str(ii)+'_full_noveto.ran.fits'
        else:
            fn = dirout+type+notqso+'_'+str(ii)+'_full_noveto.ran.fits'
        #ct.addcol_ran(fn,ii)
        logger.info('filling randoms with imaging properties to '+fn)
        new_cols=mainp.new_cols
        fid_cols=mainp.fid_cols
        common.add_map_cols(fn,ii,new_cols=new_cols,fid_cols=fid_cols,logger=logger)
        logger.info('done with '+str(ii))


    if args.apply_veto == 'y':
        logger.info('applying vetos')
        maxp = 3400
        if type[:3] == 'LRG' or notqso == 'notqso':
            maxp = 3200
        if type[:3] == 'BGS':
            maxp = 2100
        if args.fullr_mode == 'prog':
            fin = dirout+pdir+'_'+str(ii)+'_full_noveto.ran.fits'
        else:
            fin = dirout+type+notqso+'_'+str(ii)+'_full_noveto.ran.fits'
        fout = dirout+type+notqso+'_'+str(ii)+'_full.ran.fits'
        common.apply_veto(fin,fout,ebits=ebits,zmask=False,maxp=maxp,logger=logger,reccircmasks=mainp.reccircmasks)
        #logger.info('random veto '+str(ii)+' done')

    if args.add_tlcomp == 'y':
        fl = dirout+type+notqso+'_'
        ct.add_tlobs_ran(fl,ii,logger=logger)

    if mkclusran:
#         tsnrcol = 'TSNR2_ELG'
#         tsnrcut = 0
#         if type[:3] == 'ELG':
#             #dchi2 = 0.9 #This is actually the OII cut criteria for ELGs
#             tsnrcut = 80
#         if type == 'LRG':
#             #dchi2 = 16  
#             tsnrcut = 80          
#         if type[:3] == 'BGS':
#             tsnrcol = 'TSNR2_BGS'
#             dchi2 = 40
#             tsnrcut = 1000

        ct.mkclusran(dirout+type+notqso+'_',ii,zmask=zma,tsnrcut=tsnrcut,tsnrcol=tsnrcol)
    logger.info('done with random '+str(ii))
    return True
        #ct.mkclusran(dirout+type+'Alltiles_',ii,zmask=zma)
    #logf.write('ran mkclusran\n')
    #logger.info('ran mkclusran\n')
    
if __name__ == '__main__':
    if par:
        from multiprocessing import Pool
        import sys
        #N = int(sys.argv[2])
        #N = 32
        N = rx-rm#+1
        #p = Pool(N)
        inds = []
        for i in range(rm,rx):
            inds.append(i)
        #with sharedmem.MapReduce() as pool:
        #pool = sharedmem.MapReduce(np=6)
        #with pool:
        with Pool(processes=int(args.nproc)) as pool:
            #def reduce(ii, r):
            #    logger.info('chunk done '+str(ii))
            #    return r
            pool.map(doran,inds)#,reduce=reduce)

        #p.map(doran,inds)
    else:
        for i in range(rm,rx):
            doran(i)