import numpy as np
import os
import sys
from matplotlib import pyplot as plt


Nmock = 1000

import argparse

import pycorr

parser = argparse.ArgumentParser()
parser.add_argument("--tracer", help="tracer type",default='LRG')
parser.add_argument("--zmin", help="minimum redshift",default=0.4,type=float)
parser.add_argument("--zmax", help="maximum redshift",default=1.1,type=float)

parser.add_argument("--bs", help="bin size in Mpc/h, some integer multiple of 1",default=4,type=int)
parser.add_argument("--dataver", help="data version",default='test')
parser.add_argument("--njack", help="number of jack knife used",default='60')
parser.add_argument("--weight", help="weight type used for xi",default='default')
parser.add_argument("--reg", help="regions used for xi",default='NScomb_')
args = parser.parse_args()

args.rectype = None

zmin = args.zmin
zmax = args.zmax
bs = args.bs

datadir =  '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/'+args.dataver+'/xi/'
covdir = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/2.1/xi/'
try:
    covfn = covdir+'smu/ximonopole_LRG_NScomb_0.4_1.1_'+args.weight+'_lin4_cov_RascalC.txt'
    covth = np.loadtxt(covfn)
except:
    sys.exit('failed to load '+covfn)


ells = 0

def get_xi0cov():
    
    #dirm = '/global/project/projectdirs/desi/users/dvalcin/Mocks/2PCF/'
    dirm = '/global/cfs/cdirs/desi/users/dvalcin/EZMOCKS/'+args.tracer+'/Xi/'
    #fnm = 'xi_lognormal_lrg_sub_'
    fnm = 'xi_ez_'+args.tracer+'_cutsky_seed' #550_z0.6_0.8.npy
    if args.rectype is not None:
        sys.exit('no recon for EZ mocks yet')
    xinpy = dirm+fnm+'1'+'_z'+str(args.zmin)+'_'+str(args.zmax)+'.npy'
    result = pycorr.TwoPointCorrelationFunction.load(xinpy)
    rebinned = result[:(result.shape[0]//bs)*bs:bs]
    xin0 = rebinned(ells=ells)

    nbin = len(xin0)
    print(nbin)
    xiave = np.zeros((nbin))
    cov = np.zeros((nbin,nbin))

    Ntot = 0
    fac = 1.
    for i in range(1,Nmock):
        nr = str(i)
        #xii = np.loadtxt(dirm+fnm+nr+'.txt').transpose()
        xinpy = dirm+fnm+nr+'_z'+str(args.zmin)+'_'+str(args.zmax)+'.npy'
        result = pycorr.TwoPointCorrelationFunction.load(xinpy)
        rebinned = result[:(result.shape[0]//bs)*bs:bs]
        xic = rebinned(ells=ells)

        xiave += xic
        Ntot += 1.
    print( Ntot)        
    xiave = xiave/float(Ntot)
    for i in range(1,Nmock):
        nr = str(i)
        xinpy = dirm+fnm+nr+'_z'+str(args.zmin)+'_'+str(args.zmax)+'.npy'
        result = pycorr.TwoPointCorrelationFunction.load(xinpy)
        rebinned = result[:(result.shape[0]//bs)*bs:bs]
        xic = rebinned(ells=ells)

        #xii = np.loadtxt(dirm+fnm+nr+'.txt').transpose()
        #xic = xii[1]
        for j in range(0,nbin):
            xij = xic[j]#-angfac*xiit[j]
            for k in range(0,nbin):
                xik = xic[k]#-angfac*xiit[k]
                cov[j][k] += (xij-xiave[j])*(xik-xiave[k])

    cov = cov/float(Ntot)                   
        
    return cov


#data = datadir+'xi024LRGDA02_'+str(zmin)+str(zmax)+'2_default_FKPlin'+str(bs)+'.dat'
zw = ''
if zmin == 0.8 and zmax == 2.1:
    zw = 'lowz'

data = datadir +'/smu/allcounts_'+args.tracer+'_'+args.reg+str(zmin)+'_'+str(zmax)+zw+'_'+args.weight+'_lin_njack'+args.njack+'.npy'

result = pycorr.TwoPointCorrelationFunction.load(data)
rebinned = result[:(result.shape[0]//bs)*bs:bs]

s, xiell, cov = rebinned.get_corr(ells=ells, return_sep=True, return_cov=True)
stdjk = np.diag(cov)**0.5

covm = get_xi0cov() #will become covariance matrix to be used with data vector

stdm = np.diag(covm)**.5

stdth = np.diag(covth)**.5

plt.plot(s,s*stdm,'b-',label='EZ mocks')
#plt.plot(rl,rl*d[5],label='jack-knife')
plt.plot(s,s*stdjk,'r-',label='jack-knife')
plt.plot(s,s*stdth,'k-',label='analytic')
plt.xlabel('s (Mpc/h)')
plt.ylabel(r's$\sigma$')
plt.legend()
plt.title(args.tracer+' '+str(args.zmin)+'<z<'+str(args.zmax)+' '+args.weight)
plt.show()


