import numpy as np
from scipy import stats
from scipy.stats import norm
import fitsio
import glob
import os
import matplotlib.pyplot as plt
import statistics
import argparse
import astropy
from astropy.table import Table,join
from astropy.time import Time
from astropy.io import fits

import LSS.common_tools as common

parser = argparse.ArgumentParser()
#parser.add_argument("--type", help="tracer type to be selected")
basedir='/global/cfs/cdirs/desi/survey/catalogs'
parser.add_argument("--basedir", help="base directory for input/output",default=basedir)
#parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='DA02')
parser.add_argument("--verspec",help="version for redshifts",default='guadalupe')

args = parser.parse_args()
basedir = args.basedir
#version = args.version
survey  = args.survey
specver = args.verspec



f=basedir+'/'+survey+'/LSS/'+specver+'/ELG_zsuccess.txt'
f1=basedir+'/'+survey+'/LSS/'+specver+'/QSO_zsuccess.txt'
f2=basedir+'/'+survey+'/LSS/'+specver+'/LRG_zsuccess.txt'
f3=basedir+'/'+survey+'/LSS/'+specver+'/BGS_ANY_zsuccess.txt'
  


ELG=Table()
ELG['FIBER']=np.loadtxt(f)[:,0]
ELG['frac_suc']=np.loadtxt(f)[:,1]
ELG['n_suc']=np.loadtxt(f)[:,2]
ELG['n_tot']=np.loadtxt(f)[:,3]

QSO=Table()
QSO['FIBER']=np.loadtxt(f1)[:,0]
QSO['frac_suc']=np.loadtxt(f1)[:,1]
QSO['n_suc']=np.loadtxt(f1)[:,2]
QSO['n_tot']=np.loadtxt(f1)[:,3]

LRG=Table()
LRG['FIBER']=np.loadtxt(f2)[:,0]
LRG['frac_suc']=np.loadtxt(f2)[:,1]
LRG['n_suc']=np.loadtxt(f2)[:,2]
LRG['n_tot']=np.loadtxt(f2)[:,3]


BGS=Table()
BGS['FIBER']=np.loadtxt(f3)[:,0]
BGS['frac_suc']=np.loadtxt(f3)[:,1]
BGS['n_suc']=np.loadtxt(f3)[:,2]
BGS['n_tot']=np.loadtxt(f3)[:,3]


def fse(fiberstats):
    #masknosuc= fiberstats['n_suc']>0
    #print(fiberstats[~masknosuc])
    mask1ntots = fiberstats['n_tot']>1
    fiberstats = fiberstats[mask1ntots]
    fiberstats['frac_suc'] = fiberstats['n_suc']/fiberstats['n_tot']
    mean = np.sum(fiberstats['n_suc'])/np.sum(fiberstats['n_tot'])


    error_floor = True

    n, p = fiberstats['n_tot'].copy(), fiberstats['frac_suc'].copy()
    if error_floor:
        p1 = np.maximum(1-p, 1/n)  # error floor
    else:
        p1 = p
    fiberstats['frac_suc_err'] = np.clip(np.sqrt(n * p * (1-p))/n, np.sqrt(n * p1 * (1-p1))/n, 1)
    
    fiberstats['check'] =(mean - fiberstats['frac_suc'])/fiberstats['frac_suc_err']
    fiberstats.sort('frac_suc')

    
  

    from scipy.stats import binom
   
    bad_stats=Table()
    bad_stats["FIBER"]=fiberstats["FIBER"]
    bad_stats["n_tot"],bad_stats["n_suc"]=fiberstats['n_tot'], fiberstats['n_suc']
    bad_stats['n_fail'] = bad_stats['n_tot']-bad_stats['n_suc']
    bad_stats['frac_suc']= bad_stats['n_suc']/bad_stats['n_tot']
    bad_stats['frac_suc_err']=fiberstats['frac_suc_err'] 
    bad_stats['more_fail_p']=np.zeros(len(bad_stats))
    bad_stats['check']=fiberstats["check"]
    for fiber in fiberstats['FIBER']:
        n = fiberstats['n_tot'][fiberstats['FIBER']==fiber]
        s = fiberstats['n_suc'][fiberstats['FIBER']==fiber]
        p = mean
        bad_stats["more_fail_p"][fiberstats['FIBER']==fiber]= binom.cdf(s-1, n, p)

        
    nsigma=float(input("Enter the req sigma value\n"))
    #mcheck=fiberstats['check']>3
    mfail=bad_stats['more_fail_p']<norm.cdf(-nsigma)
    #print(bad_stats[mfail])
    #print(fiberstats[mcheck])    
    return(bad_stats[mfail], nsigma)#fiberstats,fiberstats[mcheck],)


def combine(fiberstats1,fiberstats2): #LRG+BGS
    full=Table()
    full['FIBER'] = np.arange(5000)
    fiberstats1 = join(fiberstats1,full, keys='FIBER',join_type='outer').filled(0)
    fiberstats2 = join(fiberstats2,full, keys='FIBER',join_type='outer').filled(0)
    fstats_comb = Table()
    fstats_comb['FIBER']=np.arange(5000)
    fstats_comb['n_tot']=np.arange(5000)
    fstats_comb['n_suc']=np.arange(5000)
    for fiber in fstats_comb['FIBER']:
        m1=fiberstats1['FIBER']==fiber
        m2=fiberstats2['FIBER']==fiber
        fstats_comb['n_tot'][fiber] = fiberstats1['n_tot'][m1]+fiberstats2['n_tot'][m2]
        fstats_comb['n_suc'][fiber] = fiberstats1['n_suc'][m1]+fiberstats2['n_suc'][m2]

    mask0= fstats_comb['n_tot']>1
    fstats_comb=fstats_comb[mask0]
    fstats_comb['frac_suc']=fstats_comb['n_suc']/fstats_comb['n_tot']
    return(fstats_comb)

LRGBGS=combine(LRG,BGS)




choice=1
while(choice==1):
    print("\nEnter number\n1 ELG\n2 LRG\n3 QSO\n4 BGS\n5 LRGBGS\n")
    t_t=int(input())
    if(t_t==1):
        bad_t,sig_t =fse(ELG)
        name="ELG"
    elif(t_t==2):
        bad_t,sig_t =fse(LRG)
        name="LRG"
    elif(t_t==3):
        bad_t,sig_t =fse(QSO)
        name="QSO"
    elif(t_t==4):
        bad_t,sig_t =fse(BGS)
        name="BGS"
    elif(t_t==5):
        bad_t,sig_t =fse(LRGBGS)
        name="LRGBGS"

    fn = "/global/homes/s/sidpen90/desicode/LSS/badfibfail/"+name+"_bad_fibers"+str(sig_t)+"_sigma.txt"
    np.savetxt(fn,bad_t['FIBER'],fmt='%i')
    print('saved results to '+fn)
    choice=int(input("\nPress 1 to run again or any other key to exit\n"))
