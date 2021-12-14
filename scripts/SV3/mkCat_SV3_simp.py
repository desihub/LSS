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
from desitarget.sv3 import sv3_targetmask

#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.SV3.cattools as ct
from LSS.globals import SV3 
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='everest')
parser.add_argument("--fulld", help="make the 'full' catalog containing info on everything physically reachable by a fiber",default='y')
parser.add_argument("--fullr", help="make the random files associated with the full data files",default='y')
parser.add_argument("--clus", help="make the data clustering files; these are cut to a small subset of columns",default='y')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='y')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=1) 
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')

parser.add_argument("--weightmode",help="what to use as option for completeness weights",default="probobs")
parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
parser.add_argument("--ntile",help="add any constraint on the number of overlapping tiles",default=0,type=int)
parser.add_argument("--rcut",help="add any cut on the rosette radius, use string like rmin,rmax",default=None)
parser.add_argument("--ccut",help="add some extra cut based on target info; should be string that tells cattools what to ",default=None)


#default processes the first of the 18 random files

args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
version = args.version
specrel = args.verspec
ntile = args.ntile
rcut = args.rcut
if rcut is not None:
    rcutstr = rcut.split(',')
    rcut = []
    rcut.append(float(rcutstr[0]))
    rcut.append(float(rcutstr[1]))

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'

print('running catalogs for tracer type '+type)

mkfulld = True #make the 'full' catalog containing info on everything physically reachable by a fiber
if args.fulld == 'n':
    mkfulld = False
    
if mkfulld:
    print('making "full" catalog file for data')    

rm = int(args.minr)
rx = int(args.maxr)
    
mkfullr = True #make the random files associated with the full data files
if args.fullr == 'n':
    mkfullr = False
    
if mkfullr:
    print('making full catalog for randoms, files '+str(rm)+ ' through '+str(rx))
    print('(if running all, consider doing in parallel)')    
    
mkclusdat = True
mkclusran = False
if args.clus == 'n':
    mkclusdat = False
    
if mkclusdat:
    print('making clustering catalog for data')
    
if args.clusran == 'y':
    mkclusran = True
    
if mkclusran:
    print('making clustering catalog for randoms, files '+str(rm)+ ' through '+str(rx))
    print('(if running all, consider doing in parallel)')  
        
mknz = False #get n(z) for type and all subtypes
if args.nz == 'y':
    mknz = True

if mknz:
    print('creating n(z); note, this does so for all tracer types and requires updated randoms to be done properly')


    
    

if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    prog = 'BRIGHT'

else:
    prog = 'DARK'

progl = prog.lower()

# tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')
# #change imaging bits to just what was applied to targeting
# ebits = None
# if type[:3] == 'BGS':
#     imbits = [1,13]
# else:
#     imbits = [1,12,13]
#     if type[:3] == 'LRG' or type[:3] == 'QSO':
#         ebits = [8,9,11]    
#     if type[:3] == 'ELG' or type[:3] == 'BGS':
#         ebits = [11]    

#location of targets
#tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/targets/sv3/resolve/'+progl+'/' 


SV3p = SV3(type,weightmode=args.weightmode)
mdir = SV3p.mdir+progl+'/' #location of ledgers
tdir = SV3p.tdir+progl+'/' #location of targets
mtld = SV3p.mtld
tiles = SV3p.tiles
imbits = SV3p.imbits #mask bits applied to targeting
ebits = SV3p.ebits #extra mask bits we think should be applied

#basedir for your outputs
sv3dir = basedir +'/SV3/LSS/'
if not os.path.exists(basedir +'/SV3'):
    os.mkdir(basedir +'/SV3')
if not os.path.exists(sv3dir):
    os.mkdir(sv3dir)
    print('made '+sv3dir)    
#base directory for previously compiled inputs
indir = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'

if not os.path.exists(sv3dir+'logs'):
    os.mkdir(sv3dir+'logs')
    print('made '+sv3dir+'logs')

ldirspec = sv3dir+specrel+'/'
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

indirspec = indir+specrel+'/'


if specrel == 'everest':
	specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-sv3-'+progl+'-cumulative.fits')
	fbcol = 'COADD_FIBERSTATUS'
	#wt = np.isin(specf['TILEID'],ta['TILEID']) #cut spec file to dark or bright time tiles
	#specf = specf[wt]
if specrel == 'daily':
	specf = Table.read(ldirspec+'datcomb_'+type+'_specwdup_Alltiles.fits')
	fbcol = 'FIBERSTATUS'
        
        
if mkfulld:

    azf=''
    if type[:3] == 'ELG':
        azf = SV3p.elgzf#'/global/cfs/cdirs/desi/users/raichoor/everest/sv3-elg-everest-tiles.fits'
    if type[:3] == 'QSO':
        azf = SV3p.qsozf#'/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/QSO/QSO_catalog_SV3.fits'
    #azf = '/global/homes/r/raichoor/sv3/sv3-elg-daily-thru20210521.fits' #file generated by Anand with OII flux info
    dz = indirspec+'datcomb_'+progl+'_tarspecwdup_Alltiles.fits' #new
    print(dz)
    if type == 'BGS_BRIGHT':
        bit = sv3_targetmask.bgs_mask[type]
        desitarg='SV3_BGS_TARGET'
    else:
        bit = sv3_targetmask.desi_mask[type]
        desitarg='SV3_DESI_TARGET'
    print(desitarg,progl,bit)
    bitweightfile = None
    if progl == 'dark':
        bitweightfile = SV3p.darkbitweightfile
        #bitweightfile='/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run64_2/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
        #bitweightfile='/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run64_2/BitweightsRound2/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
    if progl == 'bright':
        bitweightfile = SV3p.brightbitweightfile
    ct.mkfulldat(specf,dz,imbits,tdir,type,bit,dirout+type+notqso+'_full_noveto.dat.fits',\
    indirspec+'Alltiles_'+progl+'_tilelocs.dat.fits',azf=azf,desitarg=desitarg,\
    specver=specrel,notqso=notqso,bitweightfile=bitweightfile)

if mkfullr:
    if type == 'BGS_BRIGHT':
        bit = sv3_targetmask.bgs_mask[type]
        desitarg='SV3_BGS_TARGET'
    else:
        bit = sv3_targetmask.desi_mask[type]
        desitarg='SV3_DESI_TARGET'

    for ii in range(rm,rx):
        outf = dirout+type+notqso+'_'+str(ii)+'_full_noveto.ran.fits'
        ct.mkfullran(specf,indirspec,ii,imbits,outf,type,progl,bit,desitarg=desitarg,notqso=notqso,fbcol=fbcol)
        
    #logf.write('ran mkfullran\n')
    #print('ran mkfullran\n')

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
        tsnrcut = 800
    ct.mkclusdat(dirout+type+notqso+'_',tp=type,dchi2=dchi2,tsnrcut=tsnrcut,rcut=rcut,ntilecut=ntile,ebits=ebits,weightmd=SV3p.weightmode)
    #logf.write('ran mkclusdat\n')
    #print('ran mkclusdat\n')

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
        tsnrcut = 800
    rcols=['Z','WEIGHT']
    if type[:3] == 'BGS':
        rcols.append('flux_r_dered')

    for ii in range(rm,rx):
        ct.mkclusran(dirout+type+notqso+'_',ii,tsnrcut=tsnrcut,tsnrcol=tsnrcol,rcut=rcut,ntilecut=ntile,ebits=ebits,rcols=rcols)
    #logf.write('ran mkclusran\n')
    #print('ran mkclusran\n')
    
if mknz:
    wzm = ''
#     if zmask:
#         wzm = 'zmask_'
    if rcut is not None:
        wzm += '_rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntile > 0:
        wzm += '_ntileg'+str(ntilecut)+'_'    
    if args.ccut is not None:
        wzm += '_'+args.ccut #you could change this to however you want the file names to turn out

    regl = ['','_N','_S']
    
    for reg in regl:
        fb = dirout+type+wzm+reg
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


        