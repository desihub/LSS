'''
Take as input clustering catalog and save blinded clustering catalog.
'''

import os
import sys
import logging

import numpy as np


logger = logging.getLogger('Blinding')
# disable jax warning:
logging.getLogger("jax._src.lib.xla_bridge").setLevel(logging.ERROR)


def build_blinding(tpe='LRG', seed=751):
    from cosmoprimo.fiducial import DESI
    from mockfactory.blinding import get_cosmo_blind, CutskyCatalogBlinding

    zeff, bias = {'LRG':0.8,'ELG':1.1,'QSO':1.6}[tpe], {'LRG':2.,'ELG':1.3,'QSO':2.3}[tpe]

    # load pre-generated w0, wa blinded value. 
    # fix the seed to be same for all the tracer
    np.random.seed(seed)
    w0wa = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/w0wa_initvalues_zeffcombined_1000realisations.txt')
    w0_blind, wa_blind = w0wa[int(np.random.random() * 1000)]
    fnl_blind = np.random.uniform(low=-15, high=15, size=1)[0]

    #choose f_shift to compensate shift in monopole amplitude
    cosmo_fid = DESI()
    cosmo_shift = cosmo_fid.clone(w0_fld=w0_blind, wa_fld=wa_blind)
    fiducial_f = cosmo_fid.get_fourier().sigma8_z(zeff, of='theta_cb') / cosmo_fid.get_fourier().sigma8_z(zeff, of='delta_cb')

    #a, b, c for quadratic formula
    DM_fid, DH_fid = cosmo_fid.comoving_angular_distance(zeff), 1. / cosmo_fid.hubble_function(zeff)
    DM_shift, DH_shift = cosmo_shift.comoving_angular_distance(zeff), 1. / cosmo_shift.hubble_function(zeff)
    vol_fac =  (DM_shift**2 * DH_shift) / (DM_fid**2 * DH_fid)
    a = 0.2 / bias**2.
    b = 2 / (3 * bias)
    c = 1-(1 + 0.2 * (fiducial_f / bias)**2. + 2 / 3 * fiducial_f / bias) / vol_fac

    f_blind = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
    dfper = (f_blind - fiducial_f) / fiducial_f
    maxfper = 0.1
    if abs(dfper) > maxfper:
        dfper = maxfper * dfper / abs(dfper)
        f_blind = (1 + dfper) * fiducial_f

    cosmo_blind = get_cosmo_blind('DESI', z=zeff)
    cosmo_blind.params['w0_fld'] = w0_blind
    cosmo_blind.params['wa_fld'] = wa_blind
    cosmo_blind._derived['f'] = f_blind
    cosmo_blind._derived['fnl'] = fnl_blind   # on fixe la valeur pour de bon
    blinding = CutskyCatalogBlinding(cosmo_fid='DESI', cosmo_blind=cosmo_blind, bias=bias, z=zeff, position_type='rdz', mpicomm=mpicomm)

    return blinding


def collect_argparser():
    import argparse

    parser = argparse.ArgumentParser(description="Blindind scheme for DESI clustering catalogs.")

    parser.add_argument("--basedir", type=str, default='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/',
                        help="base directory where the catalogs clustering are saved. The mocks are here: /global/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit/Y1/mock1/LSScats")

    parser.add_argument("--type", type=str, default='LRG', help="tracer type to be selected")
    parser.add_argument("--notqso", type=str, default='n', help="if y, do not include any qso targets")
    parser.add_argument("--split_GC", type=str, default='n', help="whether to make the split NGC/SGC")
    parser.add_argument("--minr", type=int, default=0, help="minimum number for random files; use 1 for abacus mocks")
    parser.add_argument("--maxr", type=int, default=1, help="maximum for random files, default is 1; use 2 for abacus mocks") 

    parser.add_argument("--seed", type=int, required=False, default=741, help="Fix the seed, when draw the blinded parameters.")

    parser.add_argument("--apply_rsd_blinding", type=str, required=False, default='False', help="if True, do the rsd blinding.")
    parser.add_argument("--apply_bao_blinding", type=str, required=False, default='False', help="if True, do the bao blinding.")
    parser.add_argument("--apply_fnl_blinding", type=str, required=False, default='True', help="if True, do the fnl blinding.")

    parser.add_argument("--suff_output", type=str, required=False, default='-test',
                        help="If you do not want to overwrite the input catalog, set --suff_output != '' ")

    return parser.parse_args()


if __name__ == '__main__':
    from mockfactory import CutskyCatalog, setup_logging
    
    from mpi4py import MPI
    mpicomm = MPI.COMM_WORLD

    setup_logging(level=(logging.INFO if mpicomm.rank == 0 else logging.ERROR))

    args = collect_argparser()
    logger.info(args)

    # Load blinded cosmo:
    blinding = build_blinding(tpe=args.type, seed=args.seed)

    # # Load data / randoms:
    logger.info(args.basedir)

    if args.split_GC == 'y':
        regions = ['_SGC','_NGC']
    else:
        regions = ['_S','_N']

    for region in regions:
        notqso = 'notqso' if args.notqso == 'y' else ''
        logger.info(f'blinding catalogs for tracer type : {args.type}{notqso}{region}')  

        path_data = os.path.join(args.basedir, f"{args.type}{notqso}{region}_clustering.dat.fits")
        path_data_output = os.path.join(args.basedir, f"{args.type}{notqso}{region}_clustering{args.suff_output}.dat.fits")
        path_randoms = [os.path.join(args.basedir, f"{args.type}{notqso}{region}_{i}_clustering.ran.fits") for i in range(args.minr, args.maxr)]

        randoms = CutskyCatalog.read(path_randoms)
        data = CutskyCatalog.read(path_data)

        # Apply RSD blinding:
        if args.apply_rsd_blinding == 'True':
            data['RA'], data['DEC'], data['Z'] = blinding.rsd([data['RA'], data['DEC'], data['Z']], data_weights=data['WEIGHT'],
                                                            randoms_positions=[randoms['RA'], randoms['DEC'], randoms['Z']], randoms_weights=randoms['WEIGHT'])

        # Apply AP blinding:
        if args.apply_bao_blinding == 'True':
            data['RA'], data['DEC'], data['Z'] = blinding.ap([data['RA'], data['DEC'], data['Z']])

        # Apply fNL blinding: 
        if args.apply_fnl_blinding == 'True':
            data['WEIGHT'] = blinding.png([data['RA'], data['DEC'], data['Z']], data_weights=data['WEIGHT'],
                                        randoms_positions=[randoms['RA'], randoms['DEC'], randoms['Z']], randoms_weights=randoms['WEIGHT'],
                                        method='data_weights', shotnoise_correction=True)
        
        # Save blinded data:
        data.save(path_data_output)
