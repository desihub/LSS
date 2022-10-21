import os
import sys
import time
import numpy as np
import argparse
import itertools
import astropy.io.fits         as     fits 

from   cosmo                   import cosmo, volcom
from   scipy.interpolate       import interp1d
from   astropy.table           import Table, vstack
from   cartesian               import cartesian, rotate
from   runtime                 import calc_runtime
from   desi_randoms            import desi_randoms
from   findfile                import fetch_fields, findfile, overwrite_check, call_signature
from   gama_limits             import gama_limits, gama_field
from   scipy.spatial.transform import Rotation as R
from   ros_tools               import roscen, ros_limits
from   config                  import Configuration

np.random.seed(314)

def rotate2rosette(ros_ra, ros_dec, pos):
    pos        = np.array(pos, copy=True)
    
    rot        = R.from_rotvec(-np.radians(90. - ros_dec) * np.array([1, 0, 0]))
    res        = rot.apply(pos)
    
    rot        = R.from_rotvec(np.radians(ros_ra - 90.) * np.array([0, 0, 1]))
    
    resres     = rot.apply(res)
    
    return  resres

parser  = argparse.ArgumentParser(description='Calculate a set of boundary points')
parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')
parser.add_argument('-f', '--field',  type=str, help='select GAMA field [G9, G12, G15] or DESI rosette [R1...]', required=True)
parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
parser.add_argument('-s', '--survey', help='Survey, e.g. GAMA, DESI, etc.', type=str, default='gama')
parser.add_argument('--sampling',     help='Sampling rate', default=90000, type=int)
parser.add_argument('--prefix',       help='filename prefix', default='randoms')
parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
parser.add_argument('--config',       help='Path to configuration file', type=str, default=findfile('config'))
# Defaults to GAMA Gold limits. 
parser.add_argument('--zmin', type=float, help='Minimum redshift limit', default=0.039)
parser.add_argument('--zmax', type=float, help='Maximum redshift limit', default=0.263)


args     = parser.parse_args()
log      = args.log
field    = args.field.upper()
dryrun   = args.dryrun
survey   = args.survey.lower()
zmin     = args.zmin
zmax     = args.zmax
prefix   = args.prefix 
sampling = args.sampling
realz    = 0

start    = time.time()

fields   = fetch_fields(survey)

assert  field in fields, f'Provided {field} field is not compatible with those available for {survey} survey ({fields})'

opath    = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=prefix, realz=realz)

if log:
    logfile = findfile(ftype='boundary', dryrun=False, field=field, survey=survey, prefix=prefix, realz=realz, log=True)

    print(f'Logging to {logfile}')

    sys.stdout = open(logfile, 'w')
'''
config = Configuration(args.config)
config.update_attributes('boundary', args)
config.write()
'''
if args.nooverwrite:
    overwrite_check(opath, ext='BOUNDARY')
    
if args.dryrun:
    sampling   = 1000

call_signature(dryrun, sys.argv)

##  ras and decs.                                                                                                                                                              
if survey == 'gama':    
    area       = 60. 

    ra_min     = gama_limits[field]['ra_min']
    ra_max     = gama_limits[field]['ra_max']

    dec_min    = gama_limits[field]['dec_min']
    dec_max    = gama_limits[field]['dec_max']

    pairs      = {'RA': (ra_min, ra_max), 'DEC': (dec_min, dec_max), 'Z': (zmin, zmax)}
    names      = list(pairs.keys())

    randoms    = []

    for key0 in names:
        keys         = list(pairs.keys())
        keys.remove(key0)
        
        key1         = keys[0]
        key2         = keys[1]

        print('Solving for {} boundary ({}, {})'.format(key0, key1, key2))
        
        pair0        = pairs[key0]
        pair1        = pairs[key1]
        pair2        = pairs[key2]

        continuous   = np.linspace(pair1[0], pair1[1], sampling)
        continuous   = np.tile(continuous, 2)

        np.random.shuffle(continuous)
         
        continuous2  = np.linspace(pair2[0], pair2[1], sampling)
        continuous2  = np.tile(continuous2, 2)

        np.random.shuffle(continuous2)

        discrete      = pair0[0] * np.ones_like(continuous)
        discrete[::2] = pair0[1]

        np.random.shuffle(discrete)

        to_add        = Table(np.c_[discrete, continuous, continuous2], names=['BOUND_{}'.format(key0), 'BOUND_{}'.format(key1), 'BOUND_{}'.format(key2)])
        to_add        = to_add['BOUND_RA', 'BOUND_DEC', 'BOUND_Z']

        randoms.append(to_add)

    randoms = vstack(randoms)
    randoms.rename_column('BOUND_Z', 'Z')

elif survey == 'desi':
    # No requirement on NERSC HOST for boundary.
    inner = 0.20  # deg.                                                                                                                                                                            
    outer = 1.75  # deg.                                                                                                                                                                                

    # TODO/HACK?
    area  = np.pi * (outer**2. - inner**2.)
        
    ras   = np.arange(0., 360., 1.e-3)
    
    idecs = (90. - inner) * np.ones_like(ras)
    odecs = (90. - outer) * np.ones_like(ras)
        
    np.random.shuffle(ras)

    randoms = np.c_[ras, idecs]
    randoms = np.vstack((randoms, np.c_[ras, odecs]))

    randoms = Table(randoms, names=['BOUND_RA', 'BOUND_DEC'])
    randoms['Z'] = 0.2
        
    chis    = np.ones_like(randoms['Z'])
    chis   *= cosmo.comoving_distance(0.2).value # Mpc/h
    
    xyz     = cartesian(randoms['BOUND_RA'], randoms['BOUND_DEC'], randoms['Z'], rotate=False)

    rr      = int(field[1:])
    rr      = roscen[rr]

    ros_xyz = rotate2rosette(rr[0], rr[1], xyz)

    ras     = np.degrees(np.arctan2(ros_xyz[:,1], ros_xyz[:,0]))
      
    thetas  = np.degrees(np.arccos(ros_xyz[:,2] / chis))
    decs    = 90. - thetas

    to_wrap = ras < 0.0
    ras[to_wrap] += 360.

    randoms = Table(np.c_[ras, decs], names=['BOUND_RA', 'BOUND_DEC'])
    randoms['Z'] = np.random.uniform(zmin, zmax, len(randoms))

else:
    raise  NotImplementedError(f'No implementation for survey: {survey}')

if dryrun:
    nrand = 500

else:
    nrand = len(randoms)

print('Solved {:d} for field {}'.format(nrand, field))

randoms.pprint()

randoms['V']          = volcom(randoms['Z'].data, area=area) - volcom(zmin, area=area)
randoms['BOUNDID']    = np.arange(len(randoms))

randoms['FIELD']      = field
randoms['GAMA_FIELD'] = gama_field(randoms['BOUND_RA'], randoms['BOUND_DEC'])

xyz                    = cartesian(randoms['BOUND_RA'], randoms['BOUND_DEC'], randoms['Z'])

randoms['CARTESIAN_X'] = xyz[:,0]
randoms['CARTESIAN_Y'] = xyz[:,1]
randoms['CARTESIAN_Z'] = xyz[:,2]

xyz                    = rotate(randoms['BOUND_RA'], randoms['BOUND_DEC'], xyz)

randoms['ROTCARTESIAN_X'] = xyz[:,0]
randoms['ROTCARTESIAN_Y'] = xyz[:,1]
randoms['ROTCARTESIAN_Z'] = xyz[:,2]

randoms.meta = {'ZMIN': zmin,\
                'ZMAX': zmax,\
                'NBOUND': nrand,\
                'FIELD': field,\
                'SAMPLING': sampling,\
                'AREA': area}

print(randoms.meta)

randoms.meta['EXTNAME'] = 'BOUNDARY'

if os.path.isfile(opath):
    runtime = calc_runtime(start, f'Appending BOUNDARY extension to {opath}', xx=randoms)

else:
    raise  RuntimeError(f'Failed to find {opath} needed to append.')

boundary = Table(randoms, copy=True)

# https://github.com/desihub/redrock/blob/7952a4d8e2692a4a4f07b85286c4346579e447ce/py/redrock/external/desi.py#L64
randoms  = Table.read(opath)
randoms.meta['EXTNAME'] = 'RANDOMS'

header   = fits.Header()

hx       = fits.HDUList()
hx.append(fits.PrimaryHDU(header=header))
hx.append(fits.convenience.table_to_hdu(randoms))
hx.append(fits.convenience.table_to_hdu(boundary))

hx.writeto(opath, overwrite=True) 

runtime = calc_runtime(start, 'Finished'.format(opath))

if log:
    sys.stdout.close()
