#standard python
'''
Documentation by Sasha Safonova, 28 June 2022.
EXAMPLE USE
===========

To start the blind parameter generation procedure, specify a hashcode. This should be a hexadecimal code that translates
into an integer with ```int(hashcode, 36)```.
For example, the following line will apply BAO blinding using the hashcode '4x1e', which gets decoded as a random seed 229442:

```python main/apply_blinding_main.py --type LRG --survey DA02 --baoblind y --hashcode 4x1e```

The following does the same, but for RSD:

```python main/apply_blinding_main.py --type LRG --survey DA02 --rsdblind y --hashcode 4x1e```


GENERAL NOTES
=============

- If a hashcode is not specified, the random generator uses 1 as its seed.
- The latest values of expected unertainties on w0, wa, and f should be specified with the flags:
--expected_w0_uncertainty
--expected_wa_uncertainty
--expected_f_uncertainty
- By default these uncertainties, are set to be 0.05 for w0, 0.2 for wa, and 0.05 for f.
- The fiducial values are set as -1 for w0, 0 for wa, and 0.8 for f.
- To specify different fiducial values, use the following flags:
--fiducial_w0
--fiducial_wa
--fiducial_f

NOTES FOR TESTING AND VALIDATION
================================

To override the random blind parameter generation, use any combination of the following flags:
--specified_w0
--specified_wa
--specified_f

For example, to set blind w0 and wa values, but leave f as a randomly generated blind variable, run the following:

```python main/apply_blinding_main.py --type LRG --survey DA02 --rsdblind y --hashcode 4x1e --specified_w0 1.0 --specified_wa 0.0```
'''

import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
from numpy.random import MT19937
from numpy.random import RandomState, SeedSequence
import fitsio
import glob
import argparse
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt

import LSS.main.cattools as ct
from LSS.globals import main
import LSS.blinding_tools as blind

from cosmoprimo.fiducial import DESI
from cosmoprimo.utils import DistanceToRedshift
from cosmoprimo import Cosmology


if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir_in", help="base directory for input, default is location for official catalogs",default='/global/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--basedir_out", help="base directory for output, default is C(P)SCRATCH",default=os.environ[scratch])
parser.add_argument("--version", help="catalog version",default='EDAbeta')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='DA02')
parser.add_argument("--verspec",help="version for redshifts",default='guadalupe')
parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
parser.add_argument("--reg_md",help="whether to runn on split N/S or NGC/SGC",default='NS')

parser.add_argument("--split_GC",help="whether to make the split NGC/SGC",default='n')

parser.add_argument("--baoblind",help="if y, do the bao blinding shift",default='n')
parser.add_argument("--mkclusdat",help="if y, make the clustering data files after the BAO blinding (needed for RSD blinding)",default='n')
parser.add_argument("--mkclusran",help="if y, make the clustering random files after the BAO blinding (needed for RSD blinding)",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)#use 1 for abacus mocks
parser.add_argument("--maxr", help="maximum for random files, default is 1",default=1,type=int) #use 2 for abacus mocks
parser.add_argument("--dorecon",help="if y, run the recon needed for RSD blinding",default='n')
parser.add_argument("--rsdblind",help="if y, do the RSD blinding shift",default='n')
parser.add_argument("--hashcode", help="Code for the blinding procedure", default='0x1')
parser.add_argument("--fiducial_w0", help="Value for w0 in the DESI fiducial cosmology", default=-1)
parser.add_argument("--fiducial_wa", help="Value for wa in the DESI fiducial cosmology", default=0)
parser.add_argument("--fiducial_f", help="Value for the RSD parameter in the DESI fiducial cosmology", default=0.8)
parser.add_argument("--expected_w0_uncertainty", help="Expected uncertainty for w0", default=0.05)
parser.add_argument("--expected_wa_uncertainty", help="Expected uncertainty for wa", default=0.2)
parser.add_argument("--expected_f_uncertainty", help="Expected uncertainty for RSD f", default=0.05)
parser.add_argument("--specified_w0",
                    help="Specify a blind w0 value to overwrite the random blinding procedure",
                    default=None)
parser.add_argument("--specified_wa",
                    help="Specify a blind wa value to overwrite the random blinding procedure",
                    default=None)
parser.add_argument("--specified_f",
                    help="Specify a blind f value to overwrite the random blinding procedure",
                    default=None)
parser.add_argument("--fix_monopole",help="whether to choose f such that the amplitude of the monopole is fixed",default='n')

def make_parameter_blind(expected_value,
                         expected_error,
                         random_state):
    blind_offset = random_state.uniform(-5*expected_error, 5*expected_error)
    return expected_value + blind_offset

def translate_hashcode_to_seed(hashcode):
    decoded_seed = int(hashcode, 36)
    return decoded_seed


args = parser.parse_args()
print(args)

type = args.type
#basedir = args.basedir
version = args.version
specrel = args.verspec

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'

print('blinding catalogs for tracer type '+type+notqso)


if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    prog = 'BRIGHT'

else:
    prog = 'DARK'

progl = prog.lower()

mainp = main(args.type)
zmin = mainp.zmin
zmax = mainp.zmax
tsnrcol = mainp.tsnrcol  


#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
if 'mock' not in args.verspec:
    maindir = args.basedir_in +'/'+args.survey+'/LSS/'

    ldirspec = maindir+specrel+'/'

    dirin = ldirspec+'LSScats/'+version+'/'
    tsnrcut = mainp.tsnrcut
    dchi2 = mainp.dchi2
          


elif 'Y1/mock' in args.verspec: #e.g., use 'mocks/FirstGenMocks/AbacusSummit/Y1/mock1' to get the 1st mock with fiberassign
    dirin = args.basedir_in +'/'+args.survey+'/'+args.verspec+'/LSScats/'+version+'/'
    dchi2=None
    tsnrcut=0

else:
    sys.exit('verspec '+args.verspec+' not supported')
    

dirout = args.basedir_out+'/LSScats/'+version+'/blinded/'

# if not os.path.exists(args.basedir_out+'/LSScats/'):
#     os.makedirs(args.basedir_out+'/LSScats/')
#     print('made '+args.basedir_out+'/LSScats/')
# 
# if not os.path.exists(args.basedir_out+'/LSScats/'+version):
#     os.makedirs(args.basedir_out+'/LSScats/'+version)
#     print('made '+args.basedir_out+'/LSScats/'+version)

if not os.path.exists(dirout):
    os.makedirs(dirout)
    print('made '+dirout)


tp2z = {'LRG':0.8,'ELG':1.1,'QSO':1.6}
tp2bias = {'LRG':2.,'ELG':1.3,'QSO':2.3}
ztp = tp2z[args.type]
bias = tp2bias[args.type]

# Generate the blinded parameters
#Needs to be change to randomly read from file instead
rs = RandomState(MT19937(SeedSequence(translate_hashcode_to_seed(args.hashcode))))
w0_blind = make_parameter_blind(args.fiducial_w0,
                                args.expected_w0_uncertainty, rs)
wa_blind = make_parameter_blind(args.fiducial_wa,
                                args.expected_wa_uncertainty, rs)

# If blinded values have been specified, overwrite the random procedure here:
if args.specified_w0 is not None:
    w0_blind = float(args.specified_w0)

if args.specified_wa is not None:
    wa_blind = float(args.specified_wa)

if args.specified_f is not None:
    fgrowth_blind = float(args.specified_f)
else:
    if args.fix_monopole == 'n':
        fgrowth_blind = make_parameter_blind(args.fiducial_f,
                                             args.expected_f_uncertainty, rs)
    elif args.fix_monopole == 'y':
        #choose f_shift to compensate shift in monopole amplitude
        cosmo_fid = DESI()
        cosmo_shift = cosmo_fid.clone(w0_fld=w0_blind, wa_fld=wa_blind)

        DM_fid = cosmo_fid.comoving_angular_distance(ztp)
        DH_fid = 1./cosmo_fid.hubble_function(ztp)

        DM_shift = cosmo_shift.comoving_angular_distance(ztp)
        DH_shift = 1./cosmo_shift.hubble_function(ztp)


        vol_fac =  (DM_shift**2*DH_shift)/(DM_fid**2*DH_fid)

        #a, b, c for quadratic formula
        a = 0.2/bias**2.
        b = 2/(3*bias)
        c = 1-(1+0.2*(args.fiducial_f/bias)**2.+2/3*args.fiducial_f/bias)/vol_fac

        f_shift = (-b+np.sqrt(b**2.-4.*a*c))/(2*a)

        dfper = (f_shift-args.fiducial_f)/args.fiducial_f

        maxfper = 0.1
        if abs(dfper) > maxfper:
            dfper = maxfper*dfper/abs(dfper)
            print('f computed is '+str(round(f_shift,3))+', which is more than '+str(maxfper*100)+'% different than the fiducial value. Setting f to '+str(round((1+dfper)*args.fiducial_f,3)))
            f_shift = (1+dfper)*args.fiducial_f

        fgrowth_blind = f_shift

        print('f is '+str(fgrowth_blind))

# Write out the blind parameter values
to_write = [['w0', 'wa', 'f'],
            [f"{w0_blind}", f"{wa_blind}", f"{fgrowth_blind}"]]
np.savetxt(dirout + "blinded_parameters_"+type+".csv",
           to_write,
           delimiter=", ",
           fmt="%s")

#if args.reg_md == 'NS':
regl = ['_S','_N']
#if args.reg_md == 'GC':
gcl = ['_SGC','_NGC']

if args.baoblind == 'y':
    data = Table(fitsio.read(dirin+type+notqso+'_full.dat.fits'))
    outf = dirout + type+notqso+'_full.dat.fits'
    blind.apply_zshift_DE(data,outf,w0=w0_blind,wa=wa_blind,zcol='Z_not4clus')



if args.mkclusdat == 'y':
    ct.mkclusdat(dirout+type+notqso,tp=type,dchi2=dchi2,tsnrcut=tsnrcut,zmin=zmin,zmax=zmax)


if args.mkclusran == 'y':
    rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL']
    tsnrcol = 'TSNR2_ELG'
    if args.type[:3] == 'BGS':
        tsnrcol = 'TSNR2_BGS'
    for rannum in range(args.minr,args.maxr):
        ct.mkclusran(dirin+args.type+notqso+'_',dirout+args.type+notqso+'_',rannum,rcols=rcols,tsnrcut=tsnrcut,tsnrcol=tsnrcol)#,ntilecut=ntile,ccut=ccut)
        #for clustering, make rannum start from 0
        if 'Y1/mock' in args.verspec:
            for reg in regl:
                ranf = dirout+args.type+notqso+reg+'_'+str(rannum)+'_clustering.ran.fits'
                ranfm = dirout+args.type+notqso+reg+'_'+str(rannum-1)+'_clustering.ran.fits'
                os.system('mv '+ranf+' '+ranfm)

reg_md = args.reg_md

if args.split_GC == 'y':
    fb = dirout+args.type+notqso+'_'                
    ct.clusNStoGC(fb,args.maxr-args.minr)

sys.stdout.flush()

if args.dorecon == 'y':
    nran = args.maxr-args.minr
    if reg_md == 'NS':
        os.system('python recon.py --tracer '+args.type+' --prepare_blinding True --indir '+dirout+' --outdir '+dirout+' --nran '+str(nran))
    else:
        for gc in gcl:
            os.system('python recon.py --tracer '+args.type+' --prepare_blinding True --indir '+dirout+' --outdir '+dirout+' --nran '+str(nran)+' --region '+gc.strip('_'))

if args.rsdblind == 'y':
    if reg_md == 'NS':
        cl = regl
    if reg_md == 'GC':
        cl = gcl
    for reg in cl:
        fnd = dirout+type+notqso+reg+'_clustering.dat.fits'
        fndr = dirout+type+notqso+reg+'_clustering.MGrsd.dat.fits'
        data = Table(fitsio.read(fnd))
        data_real = Table(fitsio.read(fndr))

        out_file = fnd
        blind.apply_zshift_RSD(data,data_real,out_file,
                               fgrowth_fid=args.fiducial_f,
                               fgrowth_blind=fgrowth_blind,
                               comments=f"f_blind: {fgrowth_blind}, w0_blind: {w0_blind}, wa_blind: {wa_blind}")

