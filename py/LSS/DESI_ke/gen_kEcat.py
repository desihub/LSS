import os
import sys
import argparse
import runtime
import numpy as np

from   astropy.table   import Table, vstack
from   smith_kcorr     import GAMA_KCorrection
from   rest_gmr        import smith_rest_gmr
from   tmr_ecorr       import tmr_ecorr, tmr_q
from   abs_mag         import abs_mag
from   findfile        import findfile, fetch_fields, overwrite_check, write_desitable
from   multiprocessing import Pool
from   functools       import partial
from   config          import Configuration

np.random.seed(314)

def sub_kE(dat, kcorr_r, kcorr_g):
    rest_gmr_0p1, rest_gmr_0p1_warn = smith_rest_gmr(dat['ZSURV'], dat['GMR'])

    dat['REST_GMR_0P1']      = rest_gmr_0p1
    dat['REST_GMR_0P1_WARN'] = rest_gmr_0p1_warn.astype(np.int32)

    dat['REST_GMR_0P1_INDEX'] = kcorr_r.rest_gmr_index(dat['REST_GMR_0P1'], kcoeff=False)

    dat['KCORR_R0P1'] = kcorr_r.k(dat['ZSURV'], dat['REST_GMR_0P1'])
    dat['KCORR_G0P1'] = kcorr_g.k(dat['ZSURV'], dat['REST_GMR_0P1'])

    dat['KCORR_R0P0'] = kcorr_r.k_nonnative_zref(0.0, dat['ZSURV'], dat['REST_GMR_0P1'])
    dat['KCORR_G0P0'] = kcorr_g.k_nonnative_zref(0.0, dat['ZSURV'], dat['REST_GMR_0P1'])

    dat['REST_GMR_0P0'] = dat['GMR'] - (dat['KCORR_G0P0'] - dat['KCORR_R0P0'])

    dat['Q_COLOR_0P0'] = tmr_q(dat['REST_GMR_0P0'], aall=False)

    dat['EQ_ALL_0P0']   = tmr_ecorr(dat['ZSURV'], dat['REST_GMR_0P0'], aall=True)
    dat['EQ_COLOR_0P0'] = tmr_ecorr(dat['ZSURV'], dat['REST_GMR_0P0'], aall=False)

    dat['MALL_0P0']     = abs_mag(dat['DETMAG'], dat['DISTMOD'], dat['KCORR_R0P0'], dat['EQ_ALL_0P0'])
    dat['MCOLOR_0P0']   = abs_mag(dat['DETMAG'], dat['DISTMOD'], dat['KCORR_R0P0'], dat['EQ_COLOR_0P0'])
    dat['MQZERO_0P0']   = abs_mag(dat['DETMAG'], dat['DISTMOD'], dat['KCORR_R0P0'], np.zeros_like(dat['EQ_ALL_0P0']))

    dat['Z_THETA_QALL']   = dat['DISTMOD'] + dat['KCORR_R0P0'] + dat['EQ_ALL_0P0']
    dat['Z_THETA_QZERO']  = dat['DISTMOD'] + dat['KCORR_R0P0'] + np.zeros_like(dat['EQ_ALL_0P0'])
    dat['Z_THETA_QCOLOR'] = dat['DISTMOD'] + dat['KCORR_R0P0'] + dat['EQ_COLOR_0P0']

    ##  ----  DDP  ----                                                                                                                                                                                    
    ##  Note:  assumes median rest-frame colour and QALL.                                                                                                                                                  
    dat['DDPKCORR_R0P1'] = kcorr_r.k(dat['ZSURV'], dat['REST_GMR_0P1'], median=True)
    dat['DDPKCORR_G0P1'] = kcorr_g.k(dat['ZSURV'], dat['REST_GMR_0P1'], median=True)

    dat['DDPKCORR_R0P0'] = kcorr_r.k_nonnative_zref(0.0, dat['ZSURV'], dat['REST_GMR_0P1'], median=True)
    dat['DDPKCORR_G0P0'] = kcorr_g.k_nonnative_zref(0.0, dat['ZSURV'], dat['REST_GMR_0P1'], median=True)

    dat['DDPMALL_0P0']   = abs_mag(dat['DETMAG'], dat['DISTMOD'], dat['DDPKCORR_R0P0'], dat['EQ_ALL_0P0'])

    return dat

def gen_kE(log, dryrun, survey, nooverwrite, nproc=12):
    root      = os.environ['GOLD_DIR']

    fpath     = findfile(ftype='gold', dryrun=dryrun, survey=survey)
    opath     = findfile(ftype='kE',   dryrun=dryrun, survey=survey)

    if log:
        logfile = findfile(ftype='kE', dryrun=False, survey=survey, log=True)

        print(f'Logging to {logfile}')

        sys.stdout = open(logfile, 'w')

    if args.nooverwrite:
        overwrite_check(opath)

    fields    = fetch_fields(survey)

    print(f'Reading {fpath}')
    print(f'Writing {opath}')
  
    dat       = Table.read(fpath)
    dat.pprint()

    kcorr_r   = GAMA_KCorrection(band='R')
    kcorr_g   = GAMA_KCorrection(band='G')

    split_idx = np.arange(len(dat))
    split_idx = np.array_split(split_idx, 10)
    dat       = [dat[x] for x in split_idx]

    # --  Serial  --
    # dat     = [sub_kE(x, kcorr_r, kcorr_g) for x in dat]
    # dat     = vstack(dat) 

    # --  Multi-processing  --
    with Pool(nproc) as pool:
        dat = vstack(pool.map(partial(sub_kE, kcorr_r=kcorr_r, kcorr_g=kcorr_g), dat))

    ##  Stack
    dat.meta['IMMUTABLE'] = 'False'

    dat.pprint()

    print('Writing {}.'.format(opath))

    write_desitable(opath, dat)

    if log:
        sys.stdout.close()


if __name__ == '__main__':
    parser  = argparse.ArgumentParser(description='Gen kE cat.')
    parser.add_argument('--log', help='Create a log file of stdout.', action='store_true')
    parser.add_argument('-d', '--dryrun', help='Dryrun.', action='store_true')
    parser.add_argument('-s', '--survey', help='Select survey', default='gama')
    parser.add_argument('--config',       help='Path to configuration file', type=str, default=findfile('config'))
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
  
    args        = parser.parse_args()
    log         = args.log
    dryrun      = args.dryrun
    survey      = args.survey.lower()
    nooverwrite = args.nooverwrite

    config = Configuration(args.config)
    config.update_attributes('kE', args)
    config.write()

    gen_kE(log, dryrun, survey, nooverwrite)
