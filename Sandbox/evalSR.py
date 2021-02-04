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

#science requirements L2.2 - L2.5, from use dictionary for each type/tracer, BGS just copies LRGs where not obvious

densRD = {'LRG':300,'ELG':1280,'QSOlz':120,'QSOly':50,'BGS':300}
zr = {'LRG':(0.4,1.),'ELG':(0.6,1.6),'QSOlz':(0.6,2.1),'QSOly':(2.1,4.5),'BGS':(0.01,0.5),'ELGhiz':(1.2,1.6)}
rerr = {'LRG':0.0005,'ELG':0.0005,'QSOlz':0.0025,'QSOly':0.0025,'BGS':0.0005,'ELGhiz':0.0005}
serr = {'LRG':0.0002,'ELG':0.0002,'QSOlz':0.0004,'QSOly':1,'BGS':0.0002,'ELGhiz':0.0002} #no sys error requirement for lyman alpha
catfrac = {'LRG':0.05,'ELG':0.05,'QSOlz':0.05,'QSOly':0.02,'BGS':0.05,'ELGhiz':0.05}
catthresh = 1000./3e5


def add_truth(tp,release='blanc',depthfac=2,baseline=True,version='v1',bdir = '/global/cscratch1/sd/ajross/SV1/redshift_comps'):
    '''
    Add truth redshift values based on deep and return table; table is input to all functions below
    Any target missing a truth value has Z_TRUTH == 0 (probably should choose something better in the future)
    tp is LRG, ELG, QSO, or BGS_ANY (for currently available z files)
    only release that should work for now is blank (needs deep, so wouldn't work on nightly alone)
    depthfac sets the minimum depth to allow to give truth, current default is half the maximum deep depth value (could probably relax this?)
    '''
    
    f = fitsio.read(bdir+'/'+release+'/'+version+'/'+tp+'/alltiles_'+tp+'zinfo.fits') #fitsio *much* faster than using Table here
    if baseline:
        #from desitarget import targetmask
        #tarbit = targetmask.desi_mask[tp])
        from desitarget.sv1 import sv1_targetmask
        if tp == 'LRG':
            tarbit = sv1_targetmask.desi_mask['LRG_OPT']
            print('Using the SV1 LRG_OPT selection')
        if tp == 'ELG':
            tarbit = sv1_targetmask.desi_mask['ELG_FDR_GTOT']
            print('Using the SV1 ELG_FDR_GTOT selection')
        if tp == 'QSO':
            tarbit = sv1_targetmask.desi_mask['QSO_RF_4PASS']
            print('Using the SV1 QSO_RF_4PASS selection')
        if tp == 'BGS_ANY':
            tarbit = sv1_targetmask.desi_mask['BGS_ANY']
            print('Using the SV1 BGS_ANY selection')
        sel = (f['SV1_DESI_TARGET'] & tarbit) > 0
        print('fraction of targets in nominal selection is '+str(len(f[sel])/len(f)))
        f = f[sel]
    else:
        print('using the SV1 '+tp+' selection')    
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

def effvsdepth(tf,type,depth='R_DEPTH',nbin=10,lplace=(.15,.15)):
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
    zs = zr[type]
    zrsel = gzsel & (tcomp['Z'] > zs[0]) & (tcomp['Z'] < zs[1])
    tzrsel = (tcomp['Z_TRUTH'] > zs[0]) & (tcomp['Z_TRUTH'] < zs[1])
    bzsel = zrsel & (abs(dz) > catthresh*(1+tcomp['Z_TRUTH']))

    a = plt.hist(tcomp[depth],bins=nbin)
    b = plt.hist(tcomp[gzsel][depth],bins=a[1])
    c = plt.hist(tcomp[zrsel][depth],bins=a[1])
    d = plt.hist(tcomp[bzsel][depth],bins=a[1])
    e = plt.hist(tcomp[tzrsel][depth],bins=a[1])
    plt.clf()
    hv = []
    for i in range(0,len(a[1])-1):
        hv.append((a[1][i]+a[1][i+1])/2.)
    plt.plot(hv,c[0]/e[0],'r-',label='spectroscopic completenes')
    plt.plot(hv,e[0]/a[0],'b--',label='targeting completeness' )
    plt.plot(hv,d[0]/c[0],'.-',color='purple',label= 'spectroscopic contamination' )
    catreq = catfrac[type]*np.ones(len(a[1][:-1]))
    plt.plot(hv,catreq,'k:',label='catastrophic failure fraction req.')
    plt.legend(loc='lower left', bbox_to_anchor=lplace)
    plt.grid(alpha=0.5)
    plt.xlabel(depth+' effective exposure time')
    plt.ylabel('fraction')
    plt.ylim(0,1) #to fit label in
    #plt.xlim(0,10000)
    plt.title(type)
    plt.show()
    print('the target redshift range is '+str(zs[0])+'<z<'+str(zs[1]))
    print('spectroscopic completeness is defined as the fraction of redshifts obained within the target range that have no zwarn flag, divided by the number of true redshifts within the target range (it does not exclude catastrophic failures)')
    print('targeting completeness is defined as the fraction of targets with true redshifts within the target range, divided by the total number of targets; variations are only due to variations in target properties in different tiles')
    print('spectroscopic contamination is defined as the fraction of redshifts, within the target range and with no zwarn flag, that are further than 0.0033(1+z_truth) from z_truth (cutting to the target redshift range and to zwarn == 0 in both numerator and denominator)')
    

def repeatvsdchi2(tf,type,nbin=1000,rng=(9,2000),mind=500,maxd=1500,chi2x=50):    
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
    masknight &= (tf['R_DEPTH'] > mind) & (tf['R_DEPTH'] < maxd)
    tcomp = tf[masknight]
    dz = tcomp['Z'] - tcomp['Z_TRUTH']
    gzsel = tcomp['ZWARN'] == 0
    zs = zr[type]
    zrsel = gzsel & (tcomp['Z'] > zs[0]) & (tcomp['Z'] < zs[1])
    tzrsel = (tcomp['Z_TRUTH'] > zs[0]) & (tcomp['Z_TRUTH'] < zs[1])
    bzsel = zrsel & (abs(dz) > catthresh*(1+tcomp['Z_TRUTH']))
    ggzsel = zrsel & (abs(dz) < catthresh*(1+tcomp['Z_TRUTH']))
    a = plt.hist(tcomp[zrsel]['DELTACHI2'],bins=nbin,cumulative=-1,range=rng)
    b = plt.hist(tcomp[ggzsel]['DELTACHI2'],bins=a[1],cumulative=-1)
    plt.clf()
    hv = []
    #print(len(a[1]))
    for i in range(0,nbin):
        hv.append((a[1][i]+a[1][i+1])/2.)

    plt.plot(hv,b[0]/a[0],'k-',label='cumulative fraction not catastrophic')
    plt.xlim(7,chi2x)
    plt.xlabel('DELTACHI2 threshold')
    plt.ylabel(r'fraction with $\Delta z < 0.0033(1+z)$')
    plt.title(type+' using data with '+str(mind)+'< R_DEPTH < '+str(maxd))
    plt.show()
    

