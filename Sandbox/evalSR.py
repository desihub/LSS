'''
Evaluate SV results wrt SRD
'''

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large',
         'figure.facecolor':'w'} 
plt.rcParams.update(params)

#science requirements L2.2 - L2.5, from use dictionary for each type/tracer

densRD = {'LRG':300,'ELG':1280,'QSOlz':120,'QSOly':50}
zr = {'LRG':(0.4,1.),'ELG':(0.6,1.6),'QSOlz':(0.01,2.1),'QSOly':(2.1,4.5)}
rerr = {'LRG':0.0005,'ELG':0.0005,'QSOlz':0.0025,'QSOly':0.0025}
serr = {'LRG':0.0002,'ELG':0.0002,'QSOlz':0.0004,'QSOly':1} #no sys error requirement for lyman alpha
catfrac = {'LRG':0.05,'ELG':0.05,'QSOlz':0.05,'QSOly':0.02}
catthresh = 1000./3e5


def add_truth(tp,release='blanc',depthfac=2):
    '''
    Add truth redshift values based on deep and return table; table is input to all functions below
    Any target missing a truth value has Z_TRUTH == 0 (probably should choose something better in the future)
    tp is LRG, ELG, QSO, or BGS_ANY (for currently available z files)
    only release that should work for now is blank (needs deep, so wouldn't work on nightly alone)
    depthfac sets the minimum depth to allow to give truth, current default is half the maximum deep depth value (could probably relax this?)
    '''
    f = fitsio.read('/project/projectdirs/desi/users/ajross/catalogs/SV/redshift_comps/'+release+'/v0/'+tp+'/alltiles_'+tp+'zinfo.fits') #fitsio *much* faster than using Table here
    deep = f[f['subset']=='deep'] 
    min_depth = np.max(deep['R_DEPTH'])/depthfac

    #get list of truez in appropriate rows, matching Rongpu's definition for criteria required to allow truth determination
    #rows without truth have z==0
    mzl = np.zeros(len(f))
    tids = np.unique(f['TARGETID'])
    print('number of unique targets is '+str(len(tids)))
    for iid in tids:
        sf = f['TARGETID'] == iid
        fi = f[sf]
        fd = fi[fi['subset']=='deep']
        mask = fd['FIBERSTATUS']==0 # Remove FIBERSTATUS!=0 fibers
        mask &= fd['ZWARN'] & 2**9==0 # Remove "no data" fibers
        mask &= fd['ZWARN']==0
        mask &= fd['R_DEPTH'] > min_depth
        #mask &= fd['DELTACHI2'] > deep_dchi2
        if len(fd[mask]) > 0:
            mzl[sf] = fd['Z'][0]
    tf = Table(f)
    tf['Z_TRUTH'] = mzl
    return tf

def effvsdepth(tf,type,depth='R_DEPTH',nbin=20):
    '''
    input table tf should be created in add_truth
    type should be one of the ones in the above dictionaries
    nbin is number of depth bin
    '''
    #select out data for fair comparison
    masknight = tf['subset'] != 'deep'
    masknight &= tf['subset'] != 'all'
    masknight &= tf['Z_TRUTH'] != 0
    masknight &= tf['FIBERSTATUS']==0
    masknight &= tf['ZWARN'] & 2**9==0
    tcomp = tf[masknight]
    dz = tcomp['Z'] - tcomp['Z_TRUTH']
    gzsel = tcomp['ZWARN'] == 0
    a = plt.hist(tcomp[depth],bins=nbin)
    b = plt.hist(tcomp[gzsel][depth],bins=a[1])
    plt.clf()
    plt.plot(a[1][:-1],b[0]/a[0],'r-',label='zwarn==0')
    zs = zr[type]
    zrsel = gzsel & (tcomp['Z'] > zs[0]) & (tcomp['Z'] < zs[1])
    c = plt.hist(tcomp[zrsel][depth],bins=a[1])
    plt.plot(a[1][:-1],c[0]/a[0],'b--',label='zwarn==0 and '+str(zs[0])+'<z<'+str(zs[1]) )
    bzsel = gzsel & (abs(dz) > catthresh)
    d = plt.hist(tcomp[bzsel][depth],bins=a[1])
    plt.plot(a[1][:-1],d[0]/a[0],'.-',color='purple',label=r'zwarn==0 and $\Delta z >$0.0033' )
    catreq = catfrac[type]*np.ones(len(a[1][:-1]))
    plt.plot(a[1][:-1],catreq,'k:',label='catastrophic failure fraction req.')
    plt.legend()
    plt.show()
    
    

