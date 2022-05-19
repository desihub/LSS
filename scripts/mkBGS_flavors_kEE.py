#standard python
import sys
import os
import sys
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import glob
import argparse
from astropy.table import Table,join,unique,vstack

#from kcorr package, needs to be added to path
ke_code_root = '/global/homes/a/ajross/desicode/DESI_ke'
sys.path.append(ke_code_root)
os.environ['CODE_ROOT'] = ke_code_root
from   smith_kcorr     import GAMA_KCorrection
from   rest_gmr        import smith_rest_gmr
from   tmr_ecorr       import tmr_ecorr, tmr_q

#from this package
import LSS.SV3.cattools as ct
import LSS.common_tools as common

from LSS.tabulated_cosmo import TabulatedDESI
cosmo = TabulatedDESI()
dis_dc = cosmo.comoving_radial_distance

parser = argparse.ArgumentParser()
parser.add_argument("--tracer", help="tracer type to be selected; BGS_ANY or BGS_BRIGHT",default='BGS_BRIGHT')
parser.add_argument("--survey", help="e.g., SV3, DA02, main",default='SV3')
parser.add_argument("--verspec",help="version for redshifts",default='fuji')
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--clus", help="make the data clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=18) 

parser.add_argument("--mkcats", help="make the subsampled catalogs ",default='y')
parser.add_argument("--nz", help="get n(z) ",default='y')

args = parser.parse_args()

dirin = args.basedir+'/'+args.survey+ '/LSS/'+args.verspec+'/LSScats/'+args.version+'/'
dirout = dirin +'BGSsubcats/'

zw = ''
#if args.survey == 'DA02':
#    zw = 'zdone'

if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)


kcorr_r   = GAMA_KCorrection(band='R')
kcorr_g   = GAMA_KCorrection(band='G')




def dl(z):   # Luminosity distance from now to z
    return dis_dc(z)*(1.+z)

def dm(z):
    return 5.*np.log10(dl(z)) + 25.


def AbsMag(mag,z,ke=None,ee=None):
    amag = mag - dm(z)
    if ke is not None:
        amag -= ke
    if ee is not None:
        amag -= ee
    return amag

def cut_abr_ct(data,maxr=0,minr=-100,minct=-100,maxct=100,zmin=0.01,zmax=0.5):
    selz = data['Z'] > zmin
    selz &= data['Z'] < zmax
    data = data[selz]
    r_dered = 22.5 - 2.5*np.log10(data['flux_r_dered'])
    g_dered = 22.5 - 2.5*np.log10(data['flux_g_dered'])
    gmr = g_dered-r_dered

    rest_gmr_0p1, rest_gmr_0p1_warn = smith_rest_gmr(data['Z'], gmr)
    KCORR_R0P1 = kcorr_r.k(data['Z'], rest_gmr_0p1)
    KCORR_G0P1 = kcorr_g.k(data['Z'], rest_gmr_0p1)
    KCORR_R0P0 = kcorr_r.k_nonnative_zref(0.0, data['Z'], rest_gmr_0p1)
    KCORR_G0P0 = kcorr_g.k_nonnative_zref(0.0, data['Z'], rest_gmr_0p1)
    REST_GMR_0P0 = gmr - (KCORR_G0P0 - KCORR_R0P0)
    EQ_ALL_0P0   = tmr_ecorr(data['Z'], REST_GMR_0P0, aall=True)
    EQ_ALL_0P1   = tmr_ecorr(data['Z'], rest_gmr_0p1, aall=True)
    abr = r_dered -dm(data['Z'])-KCORR_R0P1-EQ_ALL_0P1 
    abg = g_dered -dm(data['Z'])
    ct = g_dered-r_dered-0.14*(data['Z']-0.1)/0.05
    sel = abr > minr
    sel &= abr < maxr
    sel &= ct > minct
    sel &= ct < maxct
    return data[sel]
    
ctc = 0.7 #rough red/blue cut
abl = [-21.5,-20.5,-19.5]
P0 = 7000
dz = 0.01
zmin = 0.1
zmax = 0.5

regl = ['_N','_S']
for reg in regl:
    if args.mkcats == 'y':
        dat = fitsio.read(dirin+args.tracer+zw+reg+'_clustering.dat.fits')
        for ab in abl:
            dato = cut_abr_ct(dat,maxr=ab)
            outf = dirout+args.tracer+zw+str(ab)+'ke'+reg+'_clustering.dat.fits'
            common.write_LSS(dato,outf)
            dato = cut_abr_ct(dat,maxr=ab,maxct=ctc)
            outf = dirout+args.tracer+zw+str(ab)+'keblue'+reg+'_clustering.dat.fits'
            common.write_LSS(dato,outf)
            dato = cut_abr_ct(dat,maxr=ab,minct=ctc)
            outf = dirout+args.tracer+zw+str(ab)+'kered'+reg+'_clustering.dat.fits'
            common.write_LSS(dato,outf)

        for rann in range(args.minr,args.maxr):
            dat = fitsio.read(dirin+args.tracer+zw+reg+'_'+str(rann)+'_clustering.ran.fits')
            for ab in abl:
                dato = cut_abr_ct(dat,maxr=ab)
                outf = dirout+args.tracer+zw+str(ab)+'ke'+reg+'_'+str(rann)+'_clustering.ran.fits'
                common.write_LSS(dato,outf)
                dato = cut_abr_ct(dat,maxr=ab,maxct=ctc)
                outf = dirout+args.tracer+zw+str(ab)+'keblue'+reg+'_'+str(rann)+'_clustering.ran.fits'
                common.write_LSS(dato,outf)
                dato = cut_abr_ct(dat,maxr=ab,minct=ctc)
                outf = dirout+args.tracer+zw+str(ab)+'kered'+reg+'_'+str(rann)+'_clustering.ran.fits'
                common.write_LSS(dato,outf)
    if args.nz== 'y':
        for ab in abl:
            for cl in ['ke','keblue','kered']:
                fb = dirout+args.tracer+zw+str(ab)+cl+reg
                fcr = dirin+args.tracer+zw+reg+'_0_clustering.ran.fits'
                fcd = fb+'_clustering.dat.fits'
                fout = fb+'_nz.txt'
                common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
                common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0)

