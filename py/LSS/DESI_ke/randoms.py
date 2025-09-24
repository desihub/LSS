import os
import sys
import time
import numpy as np
import argparse

from   cosmo             import cosmo, volcom
from   scipy.interpolate import interp1d
from   astropy.table     import Table
from   cartesian         import cartesian, rotate
from   runtime           import calc_runtime
from   desi_randoms      import desi_randoms
from   findfile          import fetch_fields, findfile, overwrite_check, call_signature
from   gama_limits       import gama_limits, gama_field
from   bitmask           import lumfn_mask, consv_mask
from   config            import Configuration


def randoms(field='G9', survey='gama', density=1., zmin=0.039, zmax=0.263, dryrun=False, prefix='', seed=314, oversample=8, realz=0):
    start   = time.time()

    fields  = fetch_fields(survey)

    assert  field in fields, f'Provided {field} field is not compatible with those available for {survey} survey ({fields})'

    opath   = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=prefix, realz=realz, oversample=oversample)

    if args.nooverwrite:
        overwrite_check(opath)

    seed    = seed + realz + 50 * oversample

    np.random.seed(seed)

    call_signature(dryrun, sys.argv)

    ##  ras and decs.                                                                                                                                                              
    if survey == 'gama':    
        Area    = 60.

        ra_min  = gama_limits[field]['ra_min']
        ra_max  = gama_limits[field]['ra_max']

        dec_min = gama_limits[field]['dec_min']
        dec_max = gama_limits[field]['dec_max']

        ctheta_min = np.cos(np.pi/2. - np.radians(dec_min))
        ctheta_max = np.cos(np.pi/2  - np.radians(dec_max))

        vol        = volcom(zmax, Area) - volcom(zmin, Area)

        nrand     = int(np.ceil(vol * density * oversample) / 2.0)
        
        cos_theta = np.random.uniform(ctheta_min, ctheta_max, nrand)
        theta     = np.arccos(cos_theta)
        decs      = np.pi/2. - theta
        decs      = np.degrees(decs)

        ras       = np.random.uniform(ra_min, ra_max, nrand)

        randoms   = Table(np.c_[ras, decs], names=['RANDOM_RA', 'RANDOM_DEC'])
        nrand     = len(randoms)
        
        if dryrun:
            # Dryrun:  2x2 sq. patch of sky.  
            # G12
            delta_deg = 0.5
 
            isin    = (randoms['RANDOM_RA'] > 180. - delta_deg) & (randoms['RANDOM_RA'] < 180. + delta_deg)
            isin   &= (randoms['RANDOM_DEC'] > 0. - delta_deg) & (randoms['RANDOM_DEC'] < 0. + delta_deg)

            allin  = isin

            # G9                                                                                                                                                                                       
            isin   = (randoms['RANDOM_RA']  > 135. - delta_deg) & (randoms['RANDOM_RA']  < 135. + delta_deg)
            isin  &= (randoms['RANDOM_DEC'] > 0. - delta_deg) & (randoms['RANDOM_DEC'] < 0. + delta_deg)

            allin |= isin

            # G15                                                                                                                                                                                         
            isin   = (randoms['RANDOM_RA']  > 217. - delta_deg) & (randoms['RANDOM_RA']  < 217. + delta_deg)
            isin  &= (randoms['RANDOM_DEC'] > 0.0 - delta_deg) & (randoms['RANDOM_DEC'] < 0.0 + delta_deg)
            
            allin |= isin

            randoms = randoms[allin]
            nrand   = len(randoms)
            
    elif survey == 'desi':
        if 'NERSC_HOST' in os.environ.keys():
            # Support to run on nersc only.
            randoms = desi_randoms(int(field[1:]), oversample=oversample, dryrun=dryrun)

            nrand   = randoms.meta['NRAND']
            Area    = randoms.meta['AREA']
            
        elif 'ddp1' in prefix:
            rpath   = findfile(ftype='randoms', dryrun=dryrun, field=field, survey=survey, prefix=None, realz=realz, oversample=oversample)
            randoms = Table.read(rpath)

            nrand   = randoms.meta['NRAND']
            Area    = randoms.meta['AREA']
            
        else:
            print(f'As you are not running on nersc, the output of this script is assumed to be present at {opath} for dryrun: {dryrun}.')
            return 0

    else:
        raise  NotImplementedError(f'No implementation for survey: {survey}')
    
    ##  Vs and zs.
    dz      = 1.e-4

    Vmin    = volcom(zmin, Area)
    Vmax    = volcom(zmax, Area)

    vol     = Vmax - Vmin

    density = nrand / vol

    rand_dir = os.path.dirname(opath)
        
    if not os.path.isdir(rand_dir):
        print('Creating {}'.format(rand_dir))

        os.makedirs(rand_dir)

    print('Volume [1e6]: {:.2f}; oversample: {:.2f};  density: {:.2e}; nrand [1e6]: {:.2f}'.format(vol/1.e6, oversample, density, nrand / 1.e6))

    zs      = np.arange(0.0, zmax+dz, dz)
    Vs      = volcom(zs, Area) 

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
    Vz      = interp1d(Vs, zs, kind='linear', copy=True, bounds_error=True, fill_value=np.NaN, assume_sorted=False)

    Vdraws  = np.random.uniform(0., 1., nrand)
    Vdraws  = Vmin + Vdraws * (Vmax - Vmin)

    zs      = Vz(Vdraws)

    print('Solved {:d} for field {}'.format(nrand, field))

    print('Applying rotation.')

    ras      = randoms['RANDOM_RA']
    decs     = randoms['RANDOM_DEC']

    randoms['Z'] = zs
    randoms['V'] = Vdraws
    randoms['RANDID'] = np.arange(len(randoms))

    randoms['FIELD']      = field

    # TODO/HACK/RESTORE
    randoms['GAMA_FIELD'] = gama_field(ras, decs)

    xyz      = cartesian(ras, decs, zs)

    randoms['CARTESIAN_X'] = xyz[:,0]
    randoms['CARTESIAN_Y'] = xyz[:,1]
    randoms['CARTESIAN_Z'] = xyz[:,2]

    xyz = rotate(randoms['RANDOM_RA'], randoms['RANDOM_DEC'], xyz)

    randoms['ROTCARTESIAN_X'] = xyz[:,0]
    randoms['ROTCARTESIAN_Y'] = xyz[:,1]
    randoms['ROTCARTESIAN_Z'] = xyz[:,2]

    '''
    elif survey == 'desi':    
        randoms['IS_BOUNDARY'][randoms['ROS_DIST']   > np.percentile(randoms['ROS_DIST'],   100. - boundary_percent)] = 1
        randoms['IS_BOUNDARY'][randoms['ROS_DIST']   < np.percentile(randoms['ROS_DIST'],   boundary_percent)]        = 1
    '''

    randoms['ZSURV']          = randoms['Z']
    randoms['CONSERVATIVE']   = np.zeros_like(randoms['FIELD'], dtype=int)

    if 'IN_D8LUMFN' not in randoms.dtype.names:
        randoms['IN_D8LUMFN'] = np.zeros_like(randoms['FIELD'], dtype=int)

    updates      = {'ZMIN':   zmin,\
                    'ZMAX':   zmax,\
                    'DZ':       dz,\
                    'NRAND': nrand,\
                    'FIELD': field,\
                    'AREA':   Area,\
                    'VOL':     vol,\
                    'RAND_DENS': density,\
                    'VOL8': (4./3.)*np.pi*(8.**3.),\
                    'OVERSAMPLE': oversample,\
                    'SEED': seed,\
                    'PREFIX': prefix,\
                    'REALZ': realz,\
                    'FPATH': opath}

    randoms.meta.update(updates)

    randoms.meta['NRAND8']      = randoms.meta['VOL8'] * randoms.meta['RAND_DENS']
    randoms.meta['NRAND8_PERR'] = np.sqrt(randoms.meta['NRAND8'])

    print(randoms.meta)

    runtime = calc_runtime(start, 'Writing {}'.format(opath), xx=randoms)

    randoms.write(opath, format='fits', overwrite=True)

    runtime = calc_runtime(start, 'Finished'.format(opath))


if __name__ == '__main__':
    parser  = argparse.ArgumentParser(description='Select GAMA field.')
    parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')
    parser.add_argument('-f', '--field',  type=str, help='select GAMA field [G9, G12, G15] or DESI rosette [R1...]', default='G9')
    parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
    parser.add_argument('-s', '--survey', help='Survey, e.g. GAMA, DESI, etc.', type=str, default='gama')
    parser.add_argument('--realz',        help='Realization', default=0, type=int)
    parser.add_argument('--prefix',       help='filename prefix', default='randoms')
    parser.add_argument('--config',       help='Path to configuration file', type=str, default=findfile('config'))
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
    parser.add_argument('--density',      help='Random density per (Mpc/h)^3', default=1.0, type=float)
    parser.add_argument('--oversample',   help='Oversampling factor for fillfactor counting.', default=8, type=int)
    parser.add_argument('--seed',         help='Random seed.', default=314, type=int)
    
    # Defaults to GAMA Gold limits. 
    parser.add_argument('--zmin', type=float, help='Minimum redshift limit', default=0.039)
    parser.add_argument('--zmax', type=float, help='Maximum redshift limit', default=0.263)

    args    = parser.parse_args()
    log     = args.log
    field   = args.field.upper()
    dryrun  = args.dryrun
    survey  = args.survey.lower()
    zmin    = args.zmin
    zmax    = args.zmax
    prefix  = args.prefix 
    realz   = args.realz
    seed    = args.seed

    density    = args.density
    oversample = args.oversample

    assert oversample in np.arange(1, 21, 1)

    if log:
        logfile = findfile(ftype='randoms', dryrun=False, field=field, survey=survey, prefix=prefix, realz=realz, log=True)

        print(f'Logging to {logfile}')

        sys.stdout = open(logfile, 'w')
    '''
    config = Configuration(args.config)
    config.update_attributes('randoms', args)
    config.write()
    '''
    for xx in [1, oversample]:        
        randoms(field=field, survey=survey, density=density, zmin=zmin, zmax=zmax, dryrun=dryrun, prefix=prefix, seed=seed, oversample=xx, realz=realz)

    if log:
        sys.stdout.close()

