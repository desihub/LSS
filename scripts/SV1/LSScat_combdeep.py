'''
combine catalogs from single SV1 deep fields
'''



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

#sys.path.append('../py')

#from this package
import LSS.mkCat_singletile.cattools as ct
import LSS.mkCat_singletile.fa4lsscat as fa
import LSS.mkCat_singletile.xitools as xt


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--release", help="version of the spectroscopic pipeline",default='cascades')
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')

args = parser.parse_args()

type = args.type
basedir = args.basedir
release = args.release
version = args.version




svdir = basedir+'/SV1/LSS/'
dirout = svdir+'LSScats/'+version+'/'

dodata = True
doran = True

if dodata:
    ffds = glob.glob(dirout+type+'8*_deep_full.dat.fits')

    dt = Table.read(ffds[0])
    for i in range(1,len(ffds)):
        dtn = Table.read(ffds[i])
        dt = vstack([dt,dtn])
    dt.write(dirout+type+'alltiles_deep_full.dat.fits',overwrite=True,format='fits')


    fcds = glob.glob(dirout+type+'8*_deep_clustering.dat.fits')
    dt = Table.read(fcds[0])
    for i in range(1,len(fcds)):
        dtn = Table.read(fcds[i])
        dt = vstack([dt,dtn])
    print('number rows in clustering data, number of unique rows')
    print(len(dt),len(np.unique(dt['TARGETID'])))
    datf = dirout+type+'alltiles_deep_clustering.dat.fits'
    dt.write(datf,overwrite=True,format='fits')

rm = 0
rx = 10

if doran:
    for i in range(rm,rx):
        ffrs = glob.glob(dirout+type+'8*_deep_'+str(i)+'_full.ran.fits')
        dt = Table.read(ffrs[0])
        for ii in range(1,len(ffrs)):
            dtn = Table.read(ffrs[ii])
            dt = vstack([dt,dtn])
        dt.write(dirout+type+'alltiles_deep_'+str(i)+'_full.ran.fits',overwrite=True,format='fits')    

        ffrs = glob.glob(dirout+type+'8*_deep_'+str(i)+'_clustering.ran.fits')
        dt = Table.read(ffrs[0])
        for ii in range(1,len(ffrs)):
            dtn = Table.read(ffrs[ii])
            dt = vstack([dt,dtn])
        print('number rows in clustering ran, number of unique rows')
        print(len(dt),len(np.unique(dt['TARGETID'])))

        ranf = dirout+type+'alltiles_deep_'+str(i)+'_clustering.ran.fits'
        dt.write(ranf,overwrite=True,format='fits')
        ct.ran_reassignz(ranf,datf)


subts = ['LRG','ELG','QSO','LRG_IR','LRG_OPT','LRG_SV_OPT','LRG_SV_IR','ELG_SV_GTOT','ELG_SV_GFIB','ELG_FDR_GTOT','ELG_FDR_GFIB','QSO_COLOR_4PASS',\
'QSO_RF_4PASS','QSO_COLOR_8PASS','QSO_RF_8PASS','BGS_ANY']
subtl = []
for subt in subts:
    if subt[:3] == type or subt == type:
        subtl.append(subt)
print(subtl)
ffd = dirout+type+'alltiles_deep_full.dat.fits'
fcd = dirout+type+'alltiles_deep_clustering.dat.fits'
fcr = dirout+type+'alltiles_deep_0_clustering.ran.fits'

for subt in subtl:
    fout = dirout+subt+'alltiles_deep_nz.dat'
    if type[:3] == 'QSO':
        zmin = 0.6
        zmax = 4.5
        dz = 0.05
        ct.mknz(fcd,fcr,subt,fout,bs=dz,zmin=zmin,zmax=zmax)
    else:   
        ct.mknz(fcd,fcr,subt,fout)

doclus = False
if doclus:
    import subprocess
    dirpcadw = os.environ['CSCRATCH']+'/pcadw/'
    dirpc = os.environ['CSCRATCH']+'/paircounts/'
    if not os.path.exists(dirpc):
        os.mkdir(dirpcadw)
    if not os.path.exists(dirpc):
        os.mkdir(dirpc)

    if type[:3] == 'ELG':
        zmin = 1.2
        zmax = 1.6
    if type == 'LRG':
        zmin = .5
        zmax = 1.1
    if type == 'QSO':
        zmin = 1.
        zmax = 2.

    rmax = 10
    gf = xt.createSourcesrd_ad(type,'alltiles','deep',zmin=zmin,zmax=zmax,datadir=dirout)
    subprocess.run(['chmod','+x','dopc'+gf+'.sh'])
    subprocess.run('./dopc'+gf+'.sh')
    for i in range(rm+1,rmax):
        gf = xt.createSourcesrd_ari(type,'alltiles','deep',i,zmin=zmin,zmax=zmax,datadir=dirout)
        subprocess.run(['chmod','+x','dopc'+gf+'.sh'])
        subprocess.run('./dopc'+gf+'.sh')
    #xt.ppxilcalc_LSDfjack_bs(type,'alltiles','deep',zmin=zmin,zmax=zmax,nran=rmax)
    xt.ppxilcalc_LSDfjack_bs(type,'alltiles','deep',zmin=zmin,zmax=zmax,bs=5,nran=rmax)
    #logf.write('computed paircounts\n')
