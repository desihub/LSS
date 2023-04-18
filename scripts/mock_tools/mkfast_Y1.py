from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import fitsio
import numpy as np
import glob
import os
import healpy as hp
import argparse
import sys


def radec2thphi(ra,dec):
    return (-dec+90.)*np.pi/180.,ra*np.pi/180.


if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--tracer", help="tracer type to use",choices=['LRG','ELG_LOPnotqso','QSO','BGS_BRIGHT'],default='LRG')
parser.add_argument("--mocktype", help="type of mock to use",choices=['Ab','EZ','EZ2gpc'],default='EZ')
parser.add_argument("--mockversion", help="version of mock type ",choices=['1stgen'],default='1stgen')
#parser.add_argument("--mockpath", help="Location of mock file(s)",default='/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/')
#parser.add_argument("--mockfile", help="formattable name of mock file(s). e.g. cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits. TYPE will be replaced with tracer type. PH will be replaced with realization number for simulation of mock.",default='cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits')
parser.add_argument("--real", help="number for the realization",default=1,type=int)
parser.add_argument("--survey", help="points to set of tiles",default='Y1')
parser.add_argument("--reg", help="region",choices=['NGC','SGC'],default='NGC')
parser.add_argument("--dataver", help="points to set of tiles",default='v0.1')
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/main/mocks/')

args = parser.parse_args()

if args.tracer[3] == 'BGS':
    prog = 'bright'
else:
    prog = 'dark'
#select mock data

if args.mockversion == '1stgen':
    zs = {'ELG':'z1.100','LRG':'z0.800','QSO':'z1.400'}
    if args.mocktype == 'Ab':
        mockpath = '/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/'    
        file_name = 'cutsky_'+args.tracer[3]+'_'+zs[args.tracer[3]]+'_AbacusSummit_base_c000_ph'+str(args.real)+'.fits'
    if args.mocktype == 'EZ':
        mockpath = '/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/'+args.tracer[:3]+'/'+zs[args.tracer[3]]+'/'
        if args.tracer == 'LRG':
            file_name = 'cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed'+str(args.real)+'_'+args.reg+'.fits'
    def mask(main=0, nz=0, Y5=0, sv3=0):
        return main * (2**3) + sv3 * (2**2) + Y5 * (2**1) + nz * (2**0)
    thepath = os.path.join(mockpath,file_name)
    data = fitsio.read(thepath,columns=['RA','DEC','Z','Z_COSMO','STATUS'])#f[1].data
    status = data['STATUS'][()]
    idx = np.arange(len(status))
    mask_main = mask(main=0, nz=1, Y5=1, sv3=0)
    if args.tracer == 'LRG':
        mask_main = mask(main=1, nz=1, Y5=1, sv3=0)
    idx_main = idx[(status & (mask_main))==mask_main]
    data = data[idx_main]

#load geometric mask

healpix_mask = fitiso.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/'+args.version+'/healpix_map_ran_comp_'+args.tracer+'.fits')

th,phi = radec2thphi(data['RA'],data['DEC'])
dpix = hp.ang2pix(1024,th,phi)

### ADD LINES TO SUBSAMPLE ###

### ADD LINES TO GET RANDOMS AND SUBSAMPLE ###

### PUT IN COMPLETENESS AS FUNCTION OF NTILE FROM https://desi.lbl.gov/trac/wiki/keyprojects/y1kp3/Y1details#Completenessfootprintstatistics

ntile_map = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/'+args.version+'/healpix_map_ntile_'+prog+'.fits')

### ADD LINES TO SUBSAMPLE AS FUNCTION OF NTILE AND ADD WEIGHT THAT COMPENSATES