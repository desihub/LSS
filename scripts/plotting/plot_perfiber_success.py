import numpy as np
#!pip install astropy
#!pip install fitsio
from scipy import stats
from scipy.stats import norm
import fitsio
import glob
import os
import sys
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
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--version",help="version for catalogs",default='test')

args = parser.parse_args()
basedir = args.basedir
survey  = args.survey
specver = args.verspec


tracers = ['QSO','LRG','ELG_LOPnotqso','BGS_BRIGHT']

#bfib = np.loadtxt(basedir+'/'+survey+'/LSS/'+specver+"/"+args.version+'/LRGbad.txt')#pars.badfib
#bfib = np.concatenate((bfib,np.loadtxt(basedir+'/'+survey+'/LSS/'+specver+"/"+args.version+'/BGSbad.txt')))
bfib = np.loadtxt(basedir+'/Y1/LSS/iron/unique_badfibers.txt')

indir = basedir+'/'+survey+'/LSS/'+specver+"/LSScats/"+args.version+'/'

def plot_all_petal(petal):
    for tp in tracers:
        if args.survey != 'SV3':
            from LSS.globals import main
            pars = main(tp,args.verspec)   

        fn = indir+tp+'_zsuccess_fromfull.txt'
        d = np.loadtxt(fn).transpose()
        fibl = d[0]
        f_succ = d[1]
        n_g= d[2]
        n_tot = d[3]
        
        err = np.sqrt(n_g*(1-f_succ))/n_tot

        fmin = petal*500
        fmax = (petal+1)*500
        sel = fibl >= fmin
        sel &= fibl < fmax
        sel_bfib = np.isin(fibl[sel],bfib)
        print(str(petal)+' '+str(np.sum(n_g[sel][sel_bfib]))+' "good" redshifts on masked fibers for tracer type '+tp)

        if tp == 'LRG':
            plt.errorbar(fibl[sel],f_succ[sel],err[sel],fmt='.r',label='LRG')
            plt.plot(fibl[sel][sel_bfib],f_succ[sel][sel_bfib],'kx',label=r'masked in Y1',zorder=1000)
        if tp == 'ELG_LOPnotqso':
            plt.errorbar(fibl[sel]+.25,f_succ[sel],err[sel],fmt='.b',label='ELG')
            plt.plot(fibl[sel][sel_bfib],f_succ[sel][sel_bfib],'kx',zorder=1000)
        if tp == 'QSO':
            plt.errorbar(fibl[sel]-0.25,f_succ[sel],err[sel],fmt='.g',label='QSO')
            plt.plot(fibl[sel][sel_bfib],f_succ[sel][sel_bfib],'kx',zorder=1000)
        if tp == 'BGS_BRIGHT':
            plt.errorbar(fibl[sel]+0.5,f_succ[sel],err[sel],fmt='.',label='BGS',color='brown')
            plt.plot(fibl[sel][sel_bfib],f_succ[sel][sel_bfib],'kx',zorder=1000)
    plt.grid()
    plt.legend()
    plt.xlabel('FIBER')
    plt.ylabel('redshift success rate')
    plt.title(args.verspec+' petal '+str(petal))
    plt.ylim(-0.05,1.05)
    plt.savefig(basedir+'/'+survey+'/LSS/'+specver+'/plots/petal'+str(petal)+'_zsuccess.png')
    plt.clf()
    #plt.show()

def plot_LRGBGS_petal(petal,ymin=0.8,ymax=1.05):
    for tp in ['LRG','BGS_BRIGHT']:
        if args.survey != 'SV3':
            from LSS.globals import main
            pars = main(tp,args.verspec)   

        fn = indir+tp+'_zsuccess_fromfull.txt'
        d = np.loadtxt(fn).transpose()
        fibl = d[0]
        f_succ = d[1]
        n_g= d[2]
        n_tot = d[3]
        
        err = np.sqrt(n_g*(1-f_succ))/n_tot

        fmin = petal*500
        fmax = (petal+1)*500
        sel = fibl >= fmin
        sel &= fibl < fmax
        #bfib = pars.badfib
        sel_bfib = np.isin(fibl[sel],bfib)
        sel_low = f_succ[sel] < ymin
        f_succ[sel][sel_low] = ymin+0.01

        if tp == 'LRG':
            plt.errorbar(fibl[sel],f_succ[sel],err[sel],fmt='.r',label='LRG')
            plt.plot(fibl[sel][sel_low],f_succ[sel][sel_low],'kv',label='true value below min ',zorder=100)
            plt.plot(fibl[sel][sel_bfib],f_succ[sel][sel_bfib],'kx',label=r'masked in Y1',zorder=1000)
        if tp == 'BGS_BRIGHT':
            plt.errorbar(fibl[sel]+0.5,f_succ[sel],err[sel],fmt='.',label='BGS',color='brown')
            plt.plot(fibl[sel][sel_low],f_succ[sel][sel_low],'kv',zorder=100)
            plt.plot(fibl[sel][sel_bfib],f_succ[sel][sel_bfib],'kx',zorder=1000)
    plt.grid()
    plt.legend()
    plt.xlabel('FIBER')
    plt.ylabel('redshift success rate')
    plt.title(args.verspec+' petal '+str(petal))
    plt.ylim(ymin,ymax)
    plt.savefig(basedir+'/'+survey+'/LSS/'+specver+'/plots/petal'+str(petal)+'LRGBGS_zsuccess.png')
    plt.clf()


for i in range(0,10):
    plot_all_petal(i)
    plot_LRGBGS_petal(i)
