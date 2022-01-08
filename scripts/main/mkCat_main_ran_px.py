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
import healpy as hp

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
import LSS.mkCat_singletile.fa4lsscat as fa
from LSS.globals import main

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='daily')
parser.add_argument("--ranmtl", help="make a random mtl file for the tile",default='n')
parser.add_argument("--rfa", help="run randoms through fiberassign",default='y')
parser.add_argument("--combhp", help="combine the random tiles together in healpix",default='y')
parser.add_argument("--combr", help="combine the random healpix files together",default='n')
parser.add_argument("--fullr", help="make the random files associated with the full data files",default='n')
parser.add_argument("--clus", help="make the data/random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--maskz", help="apply sky line mask to redshifts?",default='n')
parser.add_argument("--faver", help="version of fiberassign code to use for random; versions for main should be 5.0.0 or greater",default='5.0.0')
parser.add_argument("--par", help="run different random number in parallel?",default='y')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=18) 

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

mainp = main(type)

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
if len(mtld) > 0:
    tilel = []
    ral = []
    decl = []
    mtlt = []
    fal = []
    obsl = []
    pl = []
    fver = []
    fahal = []
    
    #for tile,pro in zip(mtld['TILEID'],mtld['PROGRAM']):
    for tile in mtld['TILEID']:
        ts = str(tile).zfill(6)
        try:
            fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
            tilel.append(tile)
            ral.append(fht['TILERA'])
            decl.append(fht['TILEDEC'])
            mtlt.append(fht['MTLTIME'])
            fal.append(fht['FA_RUN'])
            obsl.append(fht['OBSCON'])
            fav = fht['FA_VER']
            try:
                if int(fav[:1]) >= 5:
                    fav = '5.0.0'
                #else:
                #    print(fav)    
            except:
                print(fav)        
            pl.append(pr)
        except:
            print('failed to find and/or get info for tile '+ts)    
    ta = Table()
    ta['TILEID'] = tilel
    ta['RA'] = ral
    ta['DEC'] = decl
    ta['MTLTIME'] = mtlt
    ta['FA_RUN'] = fal
    ta['OBSCON'] = obsl
    ta['PROGRAM'] = pl
    #ta['FA_HA'] = fahal
    #ta['FA_VER'] = fver
    print(np.unique(fver,return_counts=True))
    #wfv = (np.array(fver) == faver)
    #mtld =  mtld[wfv]
    #ta = ta[wfv]
else:
    print('no done tiles in the MTL')

print(len(ta))

if type == 'dark' or type == 'bright':
    if specrel == 'daily':
        specfo = ldirspec+'datcomb_'+type+'_spec_zdone.fits'
        specf = Table.read(specfo)
        specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']

    if specrel == 'everest':    

        #specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+type+'-cumulative.fits')
        #wt = np.isin(mtld['TILEID'],specf['TILEID'])
        #above two lines already done above
        specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+type+'-cumulative.fits')
        wt = np.isin(specf['TILEID'],mtld['TILEID']) #cut spec file to dark or bright time tiles
        specf = specf[wt]
        print('number of TILEID in spec data being used:')
        print(len(np.unique(specf['TILEID'])))
        specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']


def doran(ii):
    #dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
    dirrt = '/global/cscratch1/sd/adamyers/forashley/dr9/2.3.0.dev5334/randoms/resolve/'  

    if mkranmtl:
        print('making random mtl files for each tile')
        ct.randomtiles_allmain_pix(ta,imin=ii,imax=ii+1,dirrt=dirrt+'randoms-1-'+str(ii))
    
    if runrfa:
        print('DID YOU DELETE THE OLD FILES!!!')
        for it in range(0,len(ta)):
            #print(it,len(mtld))    
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
                print(fav)        

            if fav == faver:
                ttemp = Table(ta[it])
                ttemp['OBSCONDITIONS'] = 516
                ttemp['IN_DESI'] = 1
                try:
                    ttemp['FA_PLAN'] = fbah['FA_PLAN']
                    ttemp['FA_HA'] = fbah['FA_HA']
                    ttemp['FIELDROT'] = fbah['FIELDROT']
                except:
                    print('did not add FA_PLAN and FIELDROT')
                #for i in range(rm,rx):
                testfbaf = randir+str(ii)+'/fba-'+str(tile).zfill(6)+'.fits'
                if os.path.isfile(testfbaf):
                    print('fba file already made')
                else:                   
                    print(ttemp)
                    print(fav,dt)
                    ttemp.write('tiletemp'+str(ii)+'.fits',format='fits', overwrite=True)
                    fa.getfatiles(randir+str(ii)+'/tilenofa-'+str(tile)+'.fits','tiletemp'+str(ii)+'.fits',dirout=randir+str(ii)+'/',dt = dt,faver=faver)

    if combhp:
        if type == 'dark' or type == 'bright':
            hpxs = foot.tiles2pix(8, tiles=ta)
            npx = 0
            for px in hpxs:
                print('combining target data for pixel '+str(px)+' '+str(npx)+' out of '+str(len(hpxs)))
                ct.combran_wdup_hp(mtld,ii,randir,type,ldirspec,specf,keepcols=kc)
                npx += 1
           
  

    if combr:
        print(len(mtld['TILEID']))
        #ct.combran(mtld,ii,randir,dirout,type,sv3_targetmask.desi_mask)
        if type == 'dark' or type == 'bright':

            
            kc = ['ZWARN','LOCATION','FIBER','COADD_FIBERSTATUS','TILEID','TILELOCID','FIBERASSIGN_X','FIBERASSIGN_Y','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT'\
            ,'MEAN_DELTA_X','MEAN_DELTA_Y','RMS_DELTA_X','RMS_DELTA_Y','MEAN_PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B'\
            ,'TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
            'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
            'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG']
            
            tc = ct.count_tiles_better('ran',type,ii,specrel=specrel)
            tc.write(ldirspec+'/rancomb_'+str(ii)+type+'_Alltilelocinfo.fits',format='fits', overwrite=True)




        
    if mkfullr:

        if specrel == 'everest':
            #specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+pdir+'-cumulative.fits')
            #wt = np.isin(specf['TILEID'],ta['TILEID']) #cut spec file to dark or bright time tiles
            #specf = specf[wt]
            fbcol = 'COADD_FIBERSTATUS'
        if specrel == 'daily':
            #specf = Table.read(ldirspec+'datcomb_'+pdir+'_specwdup_Alltiles.fits')
            fbcol = 'FIBERSTATUS'

        outf = dirout+type+notqso+'zdone_'+str(ii)+'_full.ran.fits'
        if type == 'BGS_BRIGHT':
            bit = targetmask.bgs_mask[type]
            desitarg='BGS_TARGET'
        else:
            bit = targetmask.desi_mask[type]    
            desitarg='DESI_TARGET'
        
        ct.mkfullran(ldirspec,ii,imbits,outf,type,pdir,bit,desitarg=desitarg,fbcol=fbcol,notqso=notqso)
        
    #logf.write('ran mkfullran\n')
    #print('ran mkfullran\n')


    if mkclusran:
        tsnrcol = 'TSNR2_ELG'
        tsnrcut = 0
        if type[:3] == 'ELG':
            #dchi2 = 0.9 #This is actually the OII cut criteria for ELGs
            tsnrcut = 80
        if type == 'LRG':
            #dchi2 = 16  
            tsnrcut = 80          
        if type[:3] == 'BGS':
            tsnrcol = 'TSNR2_BGS'
            dchi2 = 40
            tsnrcut = 1000

        ct.mkclusran(dirout+type+notqso+'zdone_',ii,zmask=zma,tsnrcut=tsnrcut,tsnrcol=tsnrcol)
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