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
#from desihub
#from desitarget import targetmask
#from regressis, must be installed
#from regressis import DR9Footprint
#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.common_tools as common

from LSS.globals import main
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--clusd", help="make the 'clustering' catalog intended for paircounts",default='n')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--splitGC", help="split into NGC/SGC catalogs",default='y')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18) 
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--addnbar_ran", help="just add nbar/fkp to randoms",default='n')

parser.add_argument("--compmd", help="whether the extra completeness gets added to data or random", choices=['dat','ran'],default='ran')
parser.add_argument("--kemd", help="where to get k+e corrections, from phot or fastspecfit", choices=['phot','spec'],default='phot')
parser.add_argument("--par", help="use parallel processing", default='y')

args = parser.parse_args()
print(args)

tracer = args.tracer
tp = tracer
basedir = args.basedir
version = args.version
specrel = args.verspec
rm = int(args.minr)
rx = int(args.maxr)


print('running catalogs for tracer type '+tp)

    
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


    

if tp[:3] == 'BGS' or tp == 'bright' or tp == 'MWS_ANY':
    prog = 'BRIGHT'

else:
    prog = 'DARK'

progl = prog.lower()

mainp = main(args.tracer,args.verspec,survey=args.survey)
mdir = mainp.mdir+progl+'/' #location of ledgers
tdir = mainp.tdir+progl+'/' #location of targets
#mtld = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied

tsnrcut = mainp.tsnrcut
dchi2 = mainp.dchi2
tsnrcol = mainp.tsnrcol        
zmin = mainp.zmin
zmax = mainp.zmax


#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/'+args.survey+'/LSS/'

#if not os.path.exists(maindir+'/logs'):
#    os.mkdir(maindir+'/logs')
#    print('made '+maindir+'/logs')

ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    print('made '+ldirspec)
    
if not os.path.exists(ldirspec+'LSScats'):
    os.mkdir(ldirspec+'LSScats')
    print('made '+ldirspec+'LSScats')

dirout = ldirspec+'LSScats/'+version+'/'
logfn = dirout+'log.txt'
if os.path.isfile(logfn):
    logf = open(logfn,'a')
else:
    logf = open(logfn,'w')
dirin = dirout

if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)    


tracer_clus = tracer #legacy of previous script

regl = ['_N','_S']    
#needs to happen before randoms so randoms can get z and weights
if mkclusdat:
    ct.mkclusdat(dirout+tracer,tp=tracer,dchi2=dchi2,tsnrcut=tsnrcut,zmin=zmin,zmax=zmax,compmd=args.compmd,kemd=args.kemd)#,ntilecut=ntile)


rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL']#,'WEIGHT_FKP']#,'WEIGHT_RF']
if tp[:3] == 'BGS':
    fcols = ['G','R','Z','W1','W2']
    for col in fcols:
        rcols.append('flux_'+col.lower()+'_dered')
    if args.kemd == 'phot':
        restcols = ['REST_GMR_0P1','REST_GMR_0P0','ABSMAG_RP0','ABSMAG_RP1']
    if args.kemd == 'spec':
        sys.exit('need to add connection to fastspecfit')
    for col in restcols:
        rcols.append(col)


if mkclusran:
    #tsnrcol = 'TSNR2_ELG'
    #if args.tracer[:3] == 'BGS':
    #    tsnrcol = 'TSNR2_BGS'
    ranin = dirin + args.tracer + '_'
    if args.tracer == 'BGS_BRIGHT-21.5':
        ranin = dirin + 'BGS_BRIGHT_'
    clus_arrays = []
    for reg in ['N','S']:
        clus_arrays.append(fitsio.read(dirout + tracer+'_'+reg+'_clustering.dat.fits'))
    def _parfun(rannum):
        ct.mkclusran(ranin, dirout + args.tracer + '_', rannum, rcols=rcols, tsnrcut=tsnrcut, tsnrcol=tsnrcol,clus_arrays=clus_arrays)#, ntilecut=ntile, ccut=ccut)
    nran = args.maxr-args.minr
    inds = np.arange(args.minr,args.maxr)
    if args.par == 'y':    
        from multiprocessing import Pool
        with Pool(processes=nran*2) as pool:
            res = pool.map(_parfun, inds)
    else:
        for ii in inds:
            ct.mkclusran(ranin, dirout + args.tracer + '_', ii, rcols=rcols, tsnrcut=tsnrcut, tsnrcol=tsnrcol,clus_arrays=clus_arrays)
            print(ii,clus_arrays[0].dtype.names)
    

if args.splitGC == 'y':    
    fb = dirout + args.tracer + '_'
    ct.clusNStoGC(fb, args.maxr - args.minr)

if tp == 'QSO':
    #zmin = 0.6
    #zmax = 4.5
    dz = 0.02
    P0 = 6000
    
else:    
    dz = 0.01
    #zmin = 0.01
    #zmax = 1.61

if tp[:3] == 'LRG':
    P0 = 10000
if tp[:3] == 'ELG':
    P0 = 4000
if tp[:3] == 'BGS':
    P0 = 7000

randens = 2500    

if args.nz == 'y':
    nran = args.maxr - args.minr
    gcl = ['_NGC','_SGC']
    for reg in gcl:
        fb = dirout+args.tracer+reg
        fcr = fb+'_0_clustering.ran.fits'
        fcd = fb+'_clustering.dat.fits'
        fout = fb+'_nz.txt'
        common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax,randens=randens)
        common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0,nran=nran)
