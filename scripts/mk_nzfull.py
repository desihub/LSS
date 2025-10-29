import sys
import os
import numpy as np
from matplotlib import pyplot as plt
import argparse
import LSS.common_tools as common

parser = argparse.ArgumentParser()
parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--inputdir", help="base directory for input, default is location for official catalogs", default='/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2/')


args = parser.parse_args()


zmax = 1.6
zmin = 0.01
bs = 0.01
if args.tracer[:3] == 'QSO':
	zmax = 4
	#bs = 0.02

fb = args.inputdir+args.tracer
fbr = fb
if args.tracer == 'BGS_BRIGHT-21.5':
	fbr = args.inputdir+'BGS_BRIGHT'
fcr = fbr+'_0_full_HPmapcut.ran.fits'
fcd = fb+'_full_HPmapcut.dat.fits'

common.mknz_full(fcd, fcr, args.tracer[:3], bs=bs, zmin=zmin, zmax=zmax,write='y')
nzf = np.loadtxt(fb+'_full_HPmapcut_nz.txt').transpose()
plt.plot(nzf[0],nzf[3],'k-',label='total')
plt.xlabel('redshift')
plt.ylabel('n(z) (h/Mpc)^3')

plt.grid()
if args.tracer[:3] == 'ELG':
    plt.ylim(0,0.0012)
if args.tracer[:3] == 'BGS':
    plt.yscale('log')
    plt.xlim(0,0.6)
    plt.ylim(1e-5,0.15)
if args.tracer == 'BGS_BRIGHT-21.5':
    plt.xlim(0,0.5)
plt.title(args.tracer)

regl = ['N','S']
for reg in regl:
    common.mknz_full(fcd, fcr, args.tracer[:3], bs=bs, zmin=zmin, zmax=zmax,write='y',reg=reg)
    nzf = np.loadtxt(fb+'_full_HPmapcut_'+reg+'_nz.txt').transpose()
    plt.plot(nzf[0],nzf[3],label=reg)

plt.legend()
plt.savefig(args.inputdir+'plots/'+args.tracer+'_nz.png')
