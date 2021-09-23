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
from desitarget import targetmask

#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.imaging.select_samples as ss
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='everest')
parser.add_argument("--redotar", help="remake the target file for the particular type (needed if, e.g., the requested columns are changed)",default='n')
parser.add_argument("--fulld", help="make the 'full' catalog containing info on everything physically reachable by a fiber",default='y')
parser.add_argument("--fillran", help="add imaging properties to randoms",default='n')
parser.add_argument("--clusd", help="make the 'clustering' catalog intended for paircounts",default='y')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='y')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=1) 
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')

parser.add_argument("--ntile",help="add any constraint on the number of overlapping tiles",default=0,type=int)
parser.add_argument("--ccut",help="add some extra cut based on target info; should be string that tells cattools what to ",default=None)


args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
version = args.version
specrel = args.verspec
ntile = args.ntile
ccut = args.ccut
rm = int(args.minr)
rx = int(args.maxr)


print('running catalogs for tracer type '+type)

redotar = False
if args.redotar == 'y':
    redotar = True

mkfulld = True #make the 'full' catalog containing info on everything physically reachable by a fiber
if args.fulld == 'n':
    mkfulld = False
        
    
if mkfulld:
    print('making "full" catalog file for data')    
    
    
#mkfullr = True #make the random files associated with the full data files
#if args.fullr == 'n':
#    mkfullr = False
    
#if mkfullr:
#    print('making full catalog for randoms, files '+str(rm)+ ' through '+str(rx))
#    print('(if running all, consider doing in parallel)')    
    
mkclusdat = False
mkclusran = False
if args.clusd == 'y':
    mkclusdat = True
    
if mkclusdat:
    print('making clustering catalog for data')
    
if args.clusran == 'y':
    mkclusran = True
    
if mkclusran:
    print('making clustering catalog for randoms, files '+str(rm)+ ' through '+str(rx))
    print('(if running all, consider doing in parallel)')  
    

if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    prog = 'BRIGHT'

else:
    prog = 'DARK'

progl = prog.lower()

mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
wt = mtld['FAPRGRM'] == progl
wt &= mtld['SURVEY'] == 'main'
wt &= mtld['ZDONE'] == 'true'
mtld = mtld[wt]
print('there are '+str(len(mtld))+' tiles')

imbits = [1,8,9,11,12,13]


#location of targets
tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/'+progl+'/' 

#columns to select from target sample
keys = ['RA', 'DEC', 'BRICKID', 'BRICKNAME','MORPHTYPE','DCHISQ','FLUX_G', 'FLUX_R', 'FLUX_Z','FLUX_W1','FLUX_W2','MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2','FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2','NOBS_G', 'NOBS_R', 'NOBS_Z','PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R',\
	   'GALDEPTH_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z',\
	   'MASKBITS','WISEMASK_W1','WISEMASK_W2', 'EBV', 'PHOTSYS','TARGETID','DESI_TARGET','BGS_TARGET','SHAPE_R']



#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/main/LSS/'

if not os.path.exists(maindir+'/logs'):
    os.mkdir(maindir+'/logs')
    print('made '+maindir+'/logs')

if not os.path.exists(maindir+'/LSScats'):
    os.mkdir(maindir+'/LSScats')
    print('made '+maindir+'/LSScats')

ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    print('made '+ldirspec)
    
if not os.path.exists(ldirspec+'LSScats'):
    os.mkdir(ldirspec+'LSScats')
    print('made '+ldirspec+'LSScats')

dirout = ldirspec+'LSScats/'+version+'/'
if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)    

tarver = '1.1.1'
tardir = '/global/cfs/cdirs/desi/target/catalogs/dr9/'+tarver+'/targets/main/resolve/'
tarf = maindir+type +'targetsDR9v'+tarver.strip('.')+'.fits'

mktar = True
if os.path.isfile(tarf) and redotar == False:
    mktar = False
#if type == 'BGS_BRIGHT':
#    mktar = False    

if mktar: #concatenate target files for given type, with column selection hardcoded
    ss.gather_targets(type,tardir,maindir,tarver,'main',progl,keys=keys)
        
        
if mkfulld:
    if specrel == 'everest':
        specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+progl+'-cumulative.fits')
        wt = np.isin(specf['TILEID'],mtld['TILEID']) #cut spec file to dark or bright time tiles
        specf = specf[wt]
        zmtlf = fitsio.read()
    if specrel == 'daily':
        specf = Table.read(ldirspec+'datcomb_'+type+'_specwdup_Alltiles.fits')
 
    ftar = fitsio.read(tarf)
    #azf = '/global/homes/r/raichoor/main/main-elg-daily-thru20210629.fits' #file generated by Anand with OII flux info
    azf = '/global/cfs/cdirs/desi/users/raichoor/everest/main-elg-everest-tiles.fits'

    azf=''
    if type[:3] == 'ELG':
        azf = '/global/cfs/cdirs/desi/users/raichoor/everest/main-elg-everest-tiles.fits'
    if type[:3] == 'QSO':
        azf = '/global/cscratch1/sd/edmondc/SHARE/QSO_CATALOG/QSO_catalog_MAIN.fits'

    dz = ldirspec+'datcomb_'+progl+'_tarspecwdup_zdone.fits' #new
    if type == 'BGS_BRIGHT':
        bit = targetmask.bgs_mask[type]
        desitarg='BGS_TARGET'
    else:
        bit = targetmask.desi_mask[type]
        desitarg='DESI_TARGET'
    
    ct.mkfulldat(specf,dz,imbits,ftar,type,bit,dirout+type+'zdone_full.dat.fits',ldirspec+'Alltiles_'+progl+'_tilelocs.dat.fits',azf=azf,desitarg=desitarg,specver=specrel)

#needs to happen before randoms so randoms can get z and weights
if mkclusdat:
    dchi2 = 9
    tsnrcut = 0
    if type[:3] == 'ELG':
        dchi2 = 0.9 #This is actually the OII cut criteria for ELGs
        tsnrcut = 80
    if type == 'LRG':
        dchi2 = 16  
        tsnrcut = 80  
    if type[:3] == 'BGS':
        dchi2 = 40
        tsnrcut = 1000
    ct.mkclusdat(dirout+type+'zdone_',tp=type,dchi2=dchi2,tsnrcut=tsnrcut)#,ntilecut=ntile,ccut=ccut)

if args.fillran == 'y':
    print('filling randoms with imaging properties')
    for ii in range(rm,rx):
        fn = dirout+type+'zdone_'+str(ii)+'_full.ran.fits'
        ct.addcol_ran(fn,ii)
        print('done with '+str(ii))

if mkclusran:
    print('doing clustering randoms')
    tsnrcol = 'TSNR2_ELG'
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

    for ii in range(rm,rx):
        ct.mkclusran(dirout+type+'zdone_',ii,tsnrcut=tsnrcut,tsnrcol=tsnrcol)#,ntilecut=ntile,ccut=ccut)

if args.nz == 'y':
    wzm = ''
#     if zmask:
#         wzm = 'zmask_'
#     if rcut is not None:
#         wzm += '_rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
#     if ntile > 0:
#         wzm += '_ntileg'+str(ntilecut)+'_'    
#     if ccut is not None:
#         wzm += '_'+ccut #you could change this to however you want the file names to turn out

    regl = ['','_N','_S']
    
    for reg in regl:
        fb = dirout+type+'zdone'+wzm+reg
        fcr = fb+'_0_clustering.ran.fits'
        fcd = fb+'_clustering.dat.fits'
        fout = fb+'_nz.dat'
        if type == 'QSO':
            zmin = 0.6
            zmax = 4.5
            dz = 0.05
            
        else:    
            dz = 0.02
            zmin = 0.01
            zmax = 1.61
        ct.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
        ct.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax)

        