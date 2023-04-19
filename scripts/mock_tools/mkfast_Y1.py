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
#parser.add_argument("--survey", help="points to set of tiles",default='Y1')
parser.add_argument("--reg", help="region",choices=['NGC','SGC'],default='NGC')
parser.add_argument("--dataver", help="points to set of tiles",default='v0.1')
parser.add_argument("--fastver", help="version for output",default='test')
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/fast/')

args = parser.parse_args()

if args.tracer[:3] == 'BGS':
    prog = 'bright'
else:
    prog = 'dark'
#select mock data

outdir = args.base_output+args.mockversion+'/'+args.mocktype+'/'+args.fastver+'/'+args.tracer+'/'
if not os.path.exists(outdir):
    os.makedirs(outdir)
foutname = outdir + args.tracer+'_'+str(args.real)+'_'+args.reg+'_clustering.dat.fits'

print('output will be written to '+foutname)

if args.mockversion == '1stgen':
    zs = {'ELG':'z1.100','LRG':'z0.800','QSO':'z1.400'}
    if args.mocktype == 'Ab':
        mockpath = '/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/'    
        file_name = 'cutsky_'+args.tracer[:3]+'_'+zs[args.tracer[:3]]+'_AbacusSummit_base_c000_ph'+str(args.real)+'.fits'
    if args.mocktype == 'EZ':
        mockpath = '/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/'+args.tracer[:3]+'/'+zs[args.tracer[:3]]+'/'
        if args.tracer == 'LRG':
            file_name = 'cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed'+str(args.real)+'_'+args.reg+'.fits'
        if args.tracer[:3] == 'ELG':
            file_name = 'cutsky_ELG_z1.100_EZmock_B6000G1536Z1.1N648012690_b0.345d1.45r40c0.05_seed'+str(args.real)+'_'+args.reg+'.fits'
        if args.tracer == 'QSO':
            file_name = 'cutsky_QSO_z1.400_EZmock_B6000G1536Z1.4N27395172_b0.053d1.13r0c0.6_seed'+str(args.real)+'_'+args.reg+'.fits'
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

healpix_mask = hp.read_map('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/'+args.dataver+'/healpix_map_ran_comp_'+args.tracer+'.fits')

### SUBSAMPLE ###

th,phi = radec2thphi(data['RA'],data['DEC'])
dpix = hp.ang2pix(1024,th,phi,nest=True)
mask_comp = np.zeros(len(data))
mask_comp = healpix_mask[dpix]
rans = np.random.random(len(data))
mask_keep = (rans < mask_comp)

print('mask will keep '+str(np.sum(mask_keep))+' for Y1 out of '+str(len(data)))

### ADD LINES TO GET RANDOMS AND SUBSAMPLE ###

### COMPLETENESS AS FUNCTION OF NTILE FROM https://desi.lbl.gov/trac/wiki/keyprojects/y1kp3/Y1details#Completenessfootprintstatistics

ntile_dic = {'BGS_BRIGHT':{0:0,1:0.512,2:755,3:0.889,4:0.954},'ELG_LOPnotqso':{0:0,1:0.308,2:0.400,3:0.528,4:0.658,5:0.767,6:0.853,7:0.917},\
'LRG':{0:0,1:0.585,2:0.741,3:0.869,4:0.936,5:0.968,6:0.985,7:0.995},'QSO':{0:0,1:0.794,2:0.957,3:0.989,4:0.994,5:0.996,6:0.998,7:0.999}}

#load ntile map

ntile_map = hp.read_map('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/'+args.dataver+'/healpix_map_ntile_'+prog+'.fits')

ntile = np.zeros(len(data))
ntile = ntile_map[dpix]
ntile_comp = np.zeros(len(data))
for i in range(0,len(ntile)):
    nt = ntile[i]
    nti = int(nt)
    if nti == nt:
        cp = ntile_dic[args.tracer][nti]
    else:
        cp_low = ntile_dic[args.tracer][nti]
        cp_high = ntile_dic[args.tracer][nti+1]
        cp = cp_low + (nt-nti)*(cp_high-cp_low) #just linearly interpolate
    ntile_comp[i] = cp

wts = 1/ntile_comp

rans = np.random.random(len(data))
comp_keep = (rans < ntile_comp)

print('we will keep '+str(np.sum(mask_keep&comp_keep))+' for Y1 out of '+str(len(data)))

data = Table(data)

data['WEIGHT'] = wts
data = data[mask_keep&comp_keep]

### NEED TO ADD SOME KIND OF n(z) sub-sampling ###


data.write(foutname,overwrite=True,format='fits')




### ADD LINES TO SUBSAMPLE AS FUNCTION OF NTILE AND ADD WEIGHT THAT COMPENSATES