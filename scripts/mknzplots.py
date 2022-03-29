import sys,os
import argparse

import numpy as np
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--survey", help="current choices are SV3,DA02,or main",default='SV3')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='everest')


args = parser.parse_args()
print(args)
catdir='/global/cfs/cdirs/desi/survey/catalogs/'
indir = catdir +args.survey+'/LSS/' +args.verspec+'/LSScats/'+args.version+'/'

dirout = indir+'/plots/'

if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)    


types = ['ELG','ELG_LOP','LRG','ELG_LOPnotqso','QSO','BGS_ANY','BGS_BRIGHT']
if args.survey == 'SV3':
    types = ['ELG','ELG_HIP','LRG','LRG_main','ELG_HIPnotqso','QSO','BGS_ANY','BGS_BRIGHT']



for tp in types:
    wzm = ''
    if args.survey != 'SV3':
        wzm = 'zdone'

    regl = ['_N','_S']
    cl = ['-r','-b']
    ll = ['BASS/MzLS','DECaLS']
    p = False
    for reg,c,l in zip(regl,cl,ll):
        fn = indir+tp+wzm+reg+'_nz.dat'
        if os.path.exists(fn):
            p = True
            zdat = np.loadtxt(fn).transpose()
            plt.plot(zdat[0],zdat[3],c,label=l)
            
        else:
            print('did not find '+fn)
    if p:
        if tp[:3] == 'ELG':
            plt.ylim(0,.0013)
            print(tp)
        if tp[:3] == 'BGS':
            plt.ylim(0,.05)
            plt.xlim(0,.6)
            print(tp)

        plt.legend()
        plt.xlabel('z (redshift)')
        plt.ylabel(r'$n(z)~ (h$Mpc$)^3$')
        plt.title(args.survey+' '+tp)

        
        plt.savefig(dirout+'nz'+args.survey+tp+'.png')
        plt.clf()
        