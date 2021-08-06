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
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='daily')
parser.add_argument("--verfull",help="version for full catalogs",default='1')
parser.add_argument("--cpfull",help="whether to copy over the full catalogs; necessary the first run",default='y')
parser.add_argument("--clus", help="make the data clustering files; these are cut to a small subset of columns",default='y')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='y')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=1) 
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')

parser.add_argument("--ntile",help="add any constraint on the number of overlapping tiles",default=0,type=int)
parser.add_argument("--rcut",help="add any cut on the rosette radius, use string like rmin,rmax",default=None)
parser.add_argument("--ccut",help="add some extra cut based on target info; should be string that tells cattools what to ",default=None)

#default processes the first of the 18 random files

args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
version = args.version
verfull = args.verfull
specrel = args.verspec
ntile = args.ntile
rcut = args.rcut
if rcut is not None:
    rcutstr = rcut.split(',')
    rcut = []
    rcut.append(float(rcutstr[0]))
    rcut.append(float(rcutstr[1]))
cpfull = args.cpfull
ccut = args.ccut
print('using ccut '+ccut)

print('running catalogs for tracer type '+type)


rm = int(args.minr)
rx = int(args.maxr)
    
    
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

tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')
imbits = [1,5,6,7,8,9,11,12,13]

#location of targets
tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/targets/sv3/resolve/'+progl+'/' 
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

indirfull = indir+specrel+'/LSScats/'+verfull+'/'



if cpfull == 'y':
    cpdatcom = 'cp '+indirfull+ type+'Alltiles_full.dat.fits '+dirout
    os.system(cpdatcom)
    for i in range(0,18):       
        cprancom = 'cp '+indirfull+ type+'Alltiles_'+str(i)+'_full.ran.fits '+dirout
        os.system(cprancom)
        

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
    ct.mkclusdat(dirout+type+'Alltiles_',tp=type,dchi2=dchi2,tsnrcut=tsnrcut,rcut=rcut,ntilecut=ntile,ccut=ccut)
    #logf.write('ran mkclusdat\n')
    #print('ran mkclusdat\n')

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
        ct.mkclusran(dirout+type+'Alltiles_',ii,tsnrcut=tsnrcut,tsnrcol=tsnrcol,rcut=rcut,ntilecut=ntile,ccut=ccut)
    #logf.write('ran mkclusran\n')
    #print('ran mkclusran\n')
    
if mknz:
    regl = ['','_N','_S']
    for reg in regl:
        if zma:
            reg = '_zmask'+reg
        fcr = dirout+type+'Alltiles'+reg+'_0_clustering.ran.fits'
        fcd = dirout+type+'Alltiles'+reg+'_clustering.dat.fits'
        fout = dirout+type+reg+'_nz.dat'
        if type == 'QSO':
            zmin = 0.6
            zmax = 4.5
            dz = 0.05
            ct.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
        else:    
            ct.mknz(fcd,fcr,fout,bs=0.02)


        