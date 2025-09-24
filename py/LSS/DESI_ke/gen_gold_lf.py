import os
import sys
import yaml
import runtime
import argparse
import pylab as pl
import numpy as np
import astropy.io.fits as fits

from   astropy.table    import Table, vstack
from   vmaxer           import vmaxer, vmaxer_rand
from   lumfn            import lumfn
from   lumfn_stepwise   import lumfn_stepwise
from   schechter        import schechter, named_schechter
from   renormalise_d8LF import renormalise_d8LF
from   delta8_limits    import d8_limits
from   config           import Configuration
from   findfile         import findfile, fetch_fields, overwrite_check, gather_cat, call_signature, write_desitable
from   jackknife_limits import solve_jackknife, set_jackknife


def jackknife_mean(fpath):           
    print('Appending JK mean and error to lumfn. extension.')

    with fits.open(fpath, mode='update') as hdulist:
        nphi =  0
        phis = []

        for i, hdu in enumerate(hdulist):
            # skip primary.
            if i > 0:
                phis.append(hdu.data['PHI_IVMAX'])

                nphi += 1 

        phis  = np.array(phis)

        mean  = np.mean(phis, axis=0)

        err   =  np.std(phis, axis=0)

        hdr   = hdulist['LUMFN'].header

        lumfn = hdulist['LUMFN'].data 
        lumfn = Table(lumfn, names=lumfn.names)

        lumfn['PHI_IVMAX_JK']       = mean
        lumfn['PHI_IVMAX_ERROR_JK'] = err 

        lumfn.pprint()

        lumfn = fits.BinTableHDU(lumfn, name='LUMFN', header=hdr)

        hdulist[1] = lumfn

        hdulist.flush()
        hdulist.close()

def process_cat(fpath, vmax_opath, field=None, survey='gama', rand_paths=[], extra_cols=[], bitmasks=[], fillfactor=False, conservative=False, stepwise=False, version='GAMA4'):        
    assert 'vmax' in vmax_opath

    opath = vmax_opath

    if not os.path.isfile(fpath):
        # Do not crash and burn, but proceed on gracefully. 
        print('WARNING:  Failed to find {}'.format(fpath))
        return  1

    zmax = Table.read(fpath)

    if len(zmax) == 0:
        print('Zero length catalogue, nothing to be done; Exiting.') 
        return 0
         
    found_fields = np.unique(zmax['FIELD'].data)
        
    print('Found fields: {}'.format(found_fields))
    
    minz = zmax['ZSURV'].min()
    maxz = zmax['ZSURV'].max()
    
    print('Found redshift limits: {:.3f} < z < {:.3f}'.format(minz, maxz))

    if field != None:
        assert  len(found_fields) == 1, 'ERROR: expected single-field restricted input, e.g. G9.'

    vmax  = vmaxer(zmax, minz, maxz, fillfactor=fillfactor, conservative=conservative, extra_cols=extra_cols)

    print('WARNING:  Found {:.3f}% with zmax < 0.0'.format(100. * np.mean(vmax['ZMAX'] <= 0.0)))

    vmax.meta['EXTNAME'] = 'VMAX'
    # vmax.meta['INPUT_CAT'] = fpath.replace(os.environ['GOLD_DIR'], '$GOLD_DIR')
        
    print('Writing {}.'.format(opath))

    write_desitable(opath, vmax)
    
    ##  Luminosity fn.
    opath  = opath.replace('vmax', 'lumfn')
    result = lumfn(vmax, bitmask='IN_D8LUMFN')

    result.meta['EXTNAME'] = 'LUMFN'
    # result.meta['INPUT_CAT'] = fpath.replace(os.environ['GOLD_DIR'], '$GOLD_DIR')

    write_desitable(opath, result)
    
    return  0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Gold luminosity function.')
    parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')
    parser.add_argument('--field', type=str, help='Select equatorial GAMA field: G9, G12, G15', default='G9')
    parser.add_argument('--survey', help='Select survey', default='gama')
    parser.add_argument('--density_split', help='Trigger density split luminosity function.', action='store_true')
    parser.add_argument('--dryrun', action='store_true', help='dryrun.')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
    parser.add_argument('--selfcount_volfracs', help='Apply volfrac corrections based on randoms counting themselves as ddps.', action='store_true')
    parser.add_argument('--version', help='Version', default='GAMA4')
    parser.add_argument('--conservative', help='Conservative analysis choices', action='store_true')
    
    args          = parser.parse_args()

    log           = args.log
    field         = args.field.upper()
    dryrun        = args.dryrun
    survey        = args.survey
    density_split = args.density_split
    self_count    = args.selfcount_volfracs
    version       = args.version
    conservative  = args.conservative
    
    if not density_split:
        if log:
            logfile = findfile(ftype='lumfn', dryrun=False, survey=survey, log=True)
            
            print(f'Logging to {logfile}')
                
            sys.stdout = open(logfile, 'w')

        print('Generating Gold reference LF.')

        call_signature(dryrun, sys.argv)

        # Bounded by gama gold, reference schechter limits:  
        # 0.039 < z < 0.263.
        # Note: not split by field. 

        prefix = 'randoms'
        
        # MJW/HACK:  repeated calls in this script to specify version == GAMA4? 
        fpath  = findfile(ftype='ddp',  dryrun=dryrun, survey=survey, prefix=prefix, version=version)
        opath  = findfile(ftype='vmax', dryrun=dryrun, survey=survey, prefix=prefix, version=version)

        if args.nooverwrite:
            overwrite_check(opath)

        print(f'Reading: {fpath}')
        print(f'Writing: {opath}')

        process_cat(fpath, opath, survey=survey, fillfactor=False)

        vmax                           = Table.read(opath)
        rand_vmax                      = vmaxer_rand(survey=survey, ftype='randoms_bd_ddp_n8', dryrun=dryrun, prefix=prefix, conservative=conservative)

        # Solve for jack knife limits.
        njack, jk_volfrac, limits, jks = solve_jackknife(rand_vmax)

        rand_vmax['JK']                = jks
        rand_vmax.meta['NJACK']        = njack
        rand_vmax.meta['JK_VOLFRAC']   = jk_volfrac

        # Set jack knife limits to data.
        vmax['JK']                     = set_jackknife(vmax['RA'], vmax['DEC'], limits=limits, debug=False)
        vmax.meta['NJACK']             = njack
        vmax.meta['JK_VOLFRAC']        = jk_volfrac

        # Save jack knife limits.
        jpath                          = findfile(ftype='jackknife', prefix=prefix, dryrun=dryrun)

        with open(jpath, 'w') as ofile:
            yaml.dump(dict(limits), ofile, default_flow_style=False)

        print(f'Writing: {jpath}')

        lpath                          = findfile(ftype='lumfn', dryrun=dryrun, survey=survey, prefix=prefix, version=version)
        jackknife                      = np.arange(njack).astype(int)

        lumfn(vmax, jackknife=jackknife, opath=lpath)

        print(f'Written {lpath}')

        jackknife_mean(lpath)

        print('Done.')

        if log:
            sys.stdout.close()

    else:
        if log:
            # TODO NOTE: Do not support version.
            logfile = findfile(ftype='ddp_n8_d0_vmax', dryrun=False, field=field, survey=survey, log=True).replace('vmax', 'lumfn').replace('_{utier}', '')
                        
            print(f'Logging to {logfile}')
        
            sys.stdout = open(logfile, 'w')

        print('Generating Gold density-split LF.')

        call_signature(dryrun, sys.argv)

        assert  field != None

        prefix = 'randoms_ddp1'

        rpath = findfile(ftype='randoms_bd_ddp_n8', dryrun=dryrun, field=field, survey=survey, prefix=prefix, version=version)
        
        if dryrun:
            # A few galaxies have a high probability to be in highest density only. 
            utiers = np.array([8])

        else:
            utiers = np.arange(len(d8_limits))
                    
        for idx in utiers:
            print(f'\n\n\n\n----------------  Solving for density tier {idx}  ----------------\n\n')

            ddp_idx   = idx + 1

            # Bounded by DDP1 z limits. 
            ddp_fpath = findfile(ftype='ddp_n8_d0', dryrun=dryrun, field=field, survey=survey, utier=idx, prefix=prefix, version=version)
            ddp_opath = findfile(ftype='ddp_n8_d0_vmax', dryrun=dryrun, field=field, survey=survey, utier=idx, prefix=prefix, version=version)
    
            print()
            print('Reading: {}'.format(ddp_fpath))

            try:
                failure   = process_cat(ddp_fpath, ddp_opath, field=field, rand_paths=[rpath], extra_cols=['MCOLOR_0P0', 'FIELD'], fillfactor=True, stepwise=False)

            except Exception as E:
                print('Error: Failed gen_gold_lf --density_split on d0 tier {:d} with Exception:'.format(idx))
                print(E)
                print('skipping.')
                
                continue 
        
            print('LF process cat. complete.')

            lpath                          = findfile(ftype='ddp_n8_d0_lumfn', field=field, dryrun=dryrun, survey=survey, utier=idx, prefix=prefix, version=version)

            result                         = Table.read(lpath)
            # result.pprint()                                                                                                                                                                          

            # Single-field values.                                                                                                                                                                       
            print('Calculating single-field volume fractions.')
            
            rand                           = Table.read(findfile(ftype='randoms_bd_ddp_n8', dryrun=dryrun, field=field, survey=survey, prefix=prefix))

            fdelta                         = float(rand.meta['DDP1_d{}_VOLFRAC'.format(idx)])
            fdelta_zp                      = float(rand.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(idx)])

            result.meta['DDP1_d{}_VOLFRAC'.format(idx)]   = '{:.6e}'.format(fdelta)
            result.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(idx)]   = '{:.6e}'.format(fdelta_zp)
                    
            # MJW:  Load three-field randoms/meta directly, for e.g. volume fractions. 
            print('Calculating multi-field volume fractions.')

            rand_vmax                      = vmaxer_rand(survey=survey, ftype='randoms_bd_ddp_n8', dryrun=dryrun, prefix=prefix, conservative=conservative)            
            rand_vmax                      = rand_vmax[rand_vmax['DDP1_DELTA8_TIER'] == idx]

            fdelta                         = float(rand_vmax.meta['DDP1_d{}_VOLFRAC'.format(idx)])
            fdelta_zp                      = float(rand_vmax.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(idx)])

            d8                             = float(rand_vmax.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(idx)])
            d8_zp                          = float(rand_vmax.meta['DDP1_d{}_TIERMEDd8'.format(idx)])

            ##  Aside for jack knife.
            print('Solving for jack knife limits.')
            
            njack, jk_volfrac, limits, jks = solve_jackknife(rand_vmax)
                
            rand_vmax['JK']                = jks
            rand_vmax.meta['NJACK']        = njack
            rand_vmax.meta['JK_VOLFRAC']   = jk_volfrac

            print('Setting data jack knife limits.')
            
            vmax_path                      = findfile(ftype='ddp_n8_d0_vmax', dryrun=False, field=field, utier=idx, survey=survey)
            vmax                           = Table.read(vmax_path, format='fits')
            
            vmax['JK']                     = set_jackknife(vmax['RA'], vmax['DEC'], limits=limits, debug=False)
            vmax.meta['NJACK']             = njack
            vmax.meta['JK_VOLFRAC']        = jk_volfrac
        
            print('Writing jack knife limits yaml')

            jpath                          = findfile(ftype='jackknife', prefix=prefix, dryrun=dryrun)
            
            with open(jpath, 'w') as jfile:
                yaml.dump(dict(limits), jfile, default_flow_style=False)

            jackknife                      = np.arange(njack).astype(int)

            print('Solving for jacked up luminosity functions.')

            lumfn(vmax, jackknife=jackknife, opath=lpath)

            print('Solving for jacked up luminosity function mean.')
            
            jackknife_mean(lpath)

            # Reload result with JK columns.
            result                         = Table.read(lpath)

            print('Renormalising LUMFN.')

            if (fdelta > 0.0) & (fdelta_zp > 0.0):
                result = renormalise_d8LF(idx, result, fdelta, fdelta_zp, self_count)
            
            else:
                assert dryrun, 'ERROR:  lf renormalisation has failed.'

            print('Solving for reference Schechter.')

            result['REF_SCHECHTER']  = named_schechter(result['MEDIAN_M'], named_type='TMR')
            result['REF_SCHECHTER'] *= (1. + d8) / (1. + 0.007)

            result['REF_RATIO']      = result['PHI_IVMAX'] / result['REF_SCHECHTER']

            print('LF renormalization and ref. schechter complete.')
            
            result.pprint()

            # Reference Schechter - finer binning
            sch_Ms = np.arange(-23., -15., 1.e-3)

            sch    = named_schechter(sch_Ms, named_type='TMR')
            sch   *= (1. + d8) / (1. + 0.007)

            ##
            ref_result = Table(np.c_[sch_Ms, sch], names=['MS', 'REFSCHECHTER'])            
            ref_result.meta['DDP1_d{}_VOLFRAC'.format(idx)]   = '{:.6e}'.format(fdelta)
            ref_result.meta['DDP1_d{}_TIERMEDd8'.format(idx)] = '{:.6e}'.format(d8)
            ref_result.meta['DDP1_d{}_ZEROPOINT_VOLFRAC'.format(idx)]   = '{:.6e}'.format(fdelta_zp)
            ref_result.meta['DDP1_d{}_ZEROPOINT_TIERMEDd8'.format(idx)] = '{:.6e}'.format(d8_zp)
            
            ##  
            keys           = sorted(result.meta.keys())
            
            header         = {}
            
            for key in keys:
                header[key] = str(result.meta[key])

            primary_hdu    = fits.PrimaryHDU()
            hdr            = fits.Header(header)
            result_hdu     = fits.BinTableHDU(result, name='LUMFN', header=hdr)
            ref_result_hdu = fits.BinTableHDU(ref_result, name='REFERENCE')
            
            # hdul         = fits.HDUList([primary_hdu, result_hdu, ref_result_hdu])

            print('Writing {}'.format(lpath))

            with fits.open(lpath, mode='update') as hdulist:
                assert  hdulist[1].header['EXTNAME'] == 'LUMFN'

                hdulist[1] = result_hdu

                for i, hdu in enumerate(hdulist):
                    hdr     = hdu.header

                    if 'EXTNAME' not in hdu.header:
                        continue

                    if 'JK' in hdu.header['EXTNAME']:
                        extname = hdu.header['EXTNAME']

                        print(f'Updating {extname}')

                        result_jk = Table(hdu.data, names=hdu.data.names)
                        result_jk = renormalise_d8LF(idx, result_jk, fdelta, fdelta_zp, self_count)
                        result_jk = fits.BinTableHDU(result_jk, name=extname, header=hdr)

                        hdulist[i] = result_jk

                hdulist.append(ref_result_hdu)
                hdulist.flush()
                hdulist.close()
            
        print('Done.')

        if log:
            sys.stdout.close()

