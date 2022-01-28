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
import LSS.common_tools as common
import LSS.imaging.select_samples as ss
from LSS.globals import main
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
parser.add_argument("--apply_veto", help="apply vetos for imaging, priorities, and hardware failures",default='n')
parser.add_argument("--fillran", help="add imaging properties to randoms",default='n')
parser.add_argument("--clusd", help="make the 'clustering' catalog intended for paircounts",default='n')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=18) 
parser.add_argument("--imsys",help="add weights for imaging systematics?",default='n')
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')

parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
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

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'

    

if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    prog = 'BRIGHT'

else:
    prog = 'DARK'

progl = prog.lower()

mainp = main(args.type)
mdir = mainp.mdir+progl+'/' #location of ledgers
tdir = mainp.tdir+progl+'/' #location of targets
#mtld = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied

#wt = mtld['FAPRGRM'] == progl
#wt &= mtld['SURVEY'] == 'main'
#wt &= mtld['ZDONE'] == 'true'
#mtld = mtld[wt]
#print('there are '+str(len(mtld))+' tiles')


#columns to select from target sample
keys = ['RA', 'DEC', 'BRICKID', 'BRICKNAME','MORPHTYPE','DCHISQ','FLUX_G', 'FLUX_R', 'FLUX_Z','FLUX_W1','FLUX_W2','MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2','FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2','NOBS_G', 'NOBS_R', 'NOBS_Z','PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R',\
       'GALDEPTH_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z',\
       'MASKBITS','WISEMASK_W1','WISEMASK_W2', 'EBV', 'PHOTSYS','TARGETID','DESI_TARGET','BGS_TARGET','SHAPE_R']



#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/main/LSS/'

if not os.path.exists(maindir+'/logs'):
    os.mkdir(maindir+'/logs')
    print('made '+maindir+'/logs')

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
    azf=''
    
    if specrel == 'everest':
        #specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+progl+'-cumulative.fits')
        #zmtlf = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/everest/datcomb_'+progl+'_zmtl_zdone.fits')
        if type[:3] == 'ELG':
            azf = mainp.elgzf
        if type[:3] == 'QSO':
            azf = mainp.qsozf
        dz = ldirspec+'datcomb_'+progl+'_tarspecwdup_zdone.fits' #new
        tlf = ldirspec+'Alltiles_'+progl+'_tilelocs.dat.fits'
    
    if specrel == 'daily':
        dz = ldirspec+'datcomb_'+type+'_tarspecwdup_zdone.fits'
        tlf = ldirspec+type+'_tilelocs.dat.fits'
        if type[:3] == 'ELG':
            azf = '/global/cfs/cdirs/desi/users/raichoor/spectro/daily/main-elg-daily-tiles-cumulative.fits'
    #if specrel == 'daily':
        #specf = Table.read(ldirspec+'datcomb_'+progl+'_spec_zdone.fits')
 
    ftar = fitsio.read(tarf)   

    
    if type == 'BGS_BRIGHT':
        bit = targetmask.bgs_mask[type]
        desitarg='BGS_TARGET'
    else:
        bit = targetmask.desi_mask[type]
        desitarg='DESI_TARGET'
    
    ct.mkfulldat(dz,imbits,ftar,type,bit,dirout+type+notqso+'zdone_full_noveto.dat.fits',tlf,azf=azf,desitarg=desitarg,specver=specrel,notqso=notqso)

if args.apply_veto == 'y':
    print('applying vetos')
    maxp = 3400
    if type[:3] == 'LRG' or notqso == 'notqso':
        maxp = 3200
    if type[:3] == 'BGS':
        maxp = 2100
    fin = dirout+type+notqso+'zdone_full_noveto.dat.fits'
    fout = dirout+type+notqso+'zdone_full.dat.fits'
    ct.apply_veto(fin,fout,ebits=ebits,zmask=False,maxp=maxp)
    print('data veto done, now doing randoms')
    for rn in range(rm,rx):
        fin = dirout+type+notqso+'zdone_'+str(rn)+'_full_noveto.ran.fits'
        fout = dirout+type+notqso+'zdone_'+str(rn)+'_full.ran.fits'
        ct.apply_veto(fin,fout,ebits=ebits,zmask=False,maxp=maxp)
        print('random veto '+str(rn)+' done')
        
    
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
    ct.mkclusdat(dirout+type+notqso+'zdone_',tp=type,dchi2=dchi2,tsnrcut=tsnrcut)#,ntilecut=ntile,ccut=ccut)

if args.fillran == 'y':
    print('filling randoms with imaging properties')
    for ii in range(rm,rx):
        fn = dirout+type+notqso+'zdone_'+str(ii)+'_full.ran.fits'
        ct.addcol_ran(fn,ii)
        print('done with '+str(ii))


if args.imsys == 'y':
    from LSS.imaging import densvar
    regl = ['_DN','_DS','','_N','_S']
    wzm = ''
    fit_maps = ['STARDENS','EBV','GALDEPTH_G', 'GALDEPTH_R','GALDEPTH_Z','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
    use_maps = fit_maps
    if type[:3] == 'ELG':
        zrl = [(0.6,0.8),(0.8,1.1),(1.1,1.5)]
    if type[:3] == 'QSO':
        zrl = [(0.8,1.6),(1.6,2.1),(2.1,3.5)]    
    if type[:3] == 'LRG':
        zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)]    
    if type[:3] == 'BGS':
        zrl = [(0.1,0.5)]    
       
        
    
    for reg in regl:
        for zr in zrl:
            zmin = zr[0]
            zmax = zr[1]
            fb = dirout+type+notqso+'zdone'+wzm+reg
            fcr = fb+'_0_clustering.ran.fits'
            rd = fitsio.read(fcr)
            fcd = fb+'_clustering.dat.fits'
            dd = Table.read(fcd)
            print('getting weights for region '+reg+' and '+str(zmin)+'<z<'+str(zmax))
            wsysl = densvar.get_imweight(dd,rd,zmin,zmax,fit_maps,use_maps,plotr=False)
            sel = wsysl != 1
            dd['WEIGHT_SYS'][sel] = wsysl[sel]
            dd['WEIGHT'][sel] *= wsysl[sel]
            dd.write(fcd,overwrite=True,format='fits')

if mkclusran:
    print('doing clustering randoms')
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
    rcols=['Z','WEIGHT']
    if type[:3] == 'BGS':
        rcols.append('flux_r_dered')

    for ii in range(rm,rx):
        ct.mkclusran(dirout+type+notqso+'zdone_',ii,rcols=rcols,tsnrcut=tsnrcut,tsnrcol=tsnrcol,ebits=ebits)#,ntilecut=ntile,ccut=ccut)

    

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

    regl = ['_DN','_DS','','_N','_S']
    
    for reg in regl:
        fb = dirout+type+notqso+'zdone'+wzm+reg
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
        common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
        common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax)

        