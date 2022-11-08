import  os
import  sys
import  time
import  tqdm
import  argparse
import  numpy           as     np

from    runtime         import calc_runtime
from    functools       import partial
from    multiprocessing import Pool
from    astropy.table   import Table
from    findfile        import findfile, overwrite_check, write_desitable
from    schechter       import named_schechter
from    ddp             import initialise_ddplimits


def lum_binner(x, dM=0.1):
    '''
    Eqn. 2.10a, W(x), of Efstathiou, Ellis & Peterson.
    '''
    return  np.abs(x) <= (dM / 2.)

def lum_visible(x, dM=0.1):
    '''                                                                                                                          
    Eqn. 2.10b, H(x), of Efstathiou, Ellis & Peterson.                                                                            
    '''

    result = -x/dM + 1./2.

    result[x >=  (dM / 2.)] = 0.0
    result[x <= -(dM / 2.)] = 1.0
    
    return  result

def process_one(split, Mmins, Mmaxs, phi_Ms, phis):
    weights     = []

    Mmins       = np.array(Mmins[split], copy=True)
    Mmaxs       = np.array(Mmaxs[split], copy=True)

    dM          = np.abs(np.diff(phi_Ms)[0])

    for i in np.arange(len(Mmins)):
        Mmin    = Mmins[i]
        Mmax    = Mmaxs[i]
            
        isin    = (phi_Ms > Mmin) & (phi_Ms < Mmax)

        assert  np.count_nonzero(isin)

        weight  = 1. / np.sum(phis[isin])
        weight /= dM

        weights.append(weight)

        # print(Mmin, Mmax, weight)

    weights = np.array(weights)

    return  weights.tolist()

def lumfn_stepwise_eval(vmax, phi_M, phi, phis, phi_Ms, dM, Mcol='MALL_0P0', survey='gama', nproc=12):
    '''
    Eqn. 2.12, of Efstathiou, Ellis & Peterson.   
    '''
    bright_curve, bright_curve_r, faint_curve, faint_curve_r = initialise_ddplimits(survey=survey)

    # HACK MALL, MQALL?
    assert  Mcol == 'MALL_0P0' 

    zmin      = 0.0 # bright_curve(phi_M) 
    zmax      = faint_curve(phi_M)

    # TODO: switch to ZSURV.
    try:
        zcol      = 'Z{}'.format(survey.upper())
        vol_lim   = vmax[(vmax[zcol] > zmin) & (vmax[zcol] < zmax)]
    except:
        zcol      = 'ZSURV'
        vol_lim   = vmax[(vmax[zcol] > zmin) & (vmax[zcol] < zmax)]

        
    Mmins     = bright_curve_r(vol_lim[zcol].data)
    Mmaxs     =  faint_curve_r(vol_lim[zcol].data)

    if len(vol_lim) == 0:
        return  np.ones_like(vol_lim[zcol].data)

    print('{:.4f}\t{:.3f}\t{:.3f}\t{:.2f}\t{:.2f}\t{:d}\t{:d}'.format(phi_M, zmin, zmax, Mmins.min(), Mmaxs.max(), len(vol_lim), len(vmax)))

    split_idx = np.arange(len(vol_lim))
    splits    = np.array_split(split_idx, 10 * nproc)

    results   = []

    with Pool(nproc) as pool:
        for result in tqdm.tqdm(pool.imap(partial(process_one, Mmins=Mmins, Mmaxs=Mmaxs, phi_Ms=phi_Ms, phis=phis), iterable=splits), total=len(splits)):
            results += result
    
    results  = np.array(results) # [1/dM]

    # pool.close()

    #  dM * phis.   
    return  results

def lumfn_stepwise(vmax, Mcol='MALL_0P0', Mmin_col='DDPMALL_0P0_VISZ', survey='gama', tolerance=1.e-3):
    dM        = 0.2
    phi_Ms    = np.arange(-25.5, -15.5, dM)

    phi_init  = dM * 1. * np.ones(len(phi_Ms), dtype=float)

    diff      = 1.e99
    phis      = phi_init 

    iteration = 0
    
    while (diff > tolerance):
        print('Solving for iteration {:d} with diff. {:.6e}'.format(iteration, diff))
        
        new_phis = []
    
        for i, (phi_M, phi) in enumerate(zip(phi_Ms, phis)):
            Ms       = vmax[Mcol]
            num      = np.count_nonzero(lum_binner(phi_M - Ms))

            weights  = lumfn_stepwise_eval(vmax, phi_M, phi, phis, phi_Ms, dM, Mcol=Mcol, survey=survey)
            
            phi_hat  = num / np.sum(weights.data)
            new_phis.append(phi_hat)

            # print(num, len(weights), np.sort(weights), phi_hat)

        new_phis    = np.array(new_phis)
        
        diff        = np.sum((new_phis - phis)**2.)
    
        #  Update previous estimate. 
        phis        = new_phis

        iteration  += 1

    phis  = phis / dM 
    nn    = dM * np.sum(phis)
    phis /= nn

    # print('Final M={} recovers weights for all galaxies in vmax ({} weights for {} galaxies).'.format(phi_M, len(weights), len(vmax)))

    weights = lumfn_stepwise_eval(vmax, -30.0, phi, phis, phi_Ms, dM, Mcol=Mcol, survey=survey)
    
    return  phi_Ms, phis, weights
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Gold stepwise luminosity function.')
    parser.add_argument('--log',          help='Create a log file of stdout.', action='store_true')
    parser.add_argument('-s', '--survey', help='Select survey', default='gama')
    parser.add_argument('--dryrun',       help='Dryrun', action='store_true')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')
    parser.add_argument('--version',      help='Add version', default='GAMA4')
    
    start       = time.time() 

    args        = parser.parse_args()
    log         = args.log
    survey      = args.survey
    dryrun      = args.dryrun
    nooverwrite = args.nooverwrite
    version     = args.version

    if log:
        logfile = findfile(ftype='lumfn_step', dryrun=False, survey=survey, log=True)

        print(f'Logging to {logfile}')

        sys.stdout = open(logfile, 'w')

    fpath       = findfile('ddp', dryrun=dryrun, survey=survey, version=version)
    opath       = findfile('lumfn_step', dryrun=dryrun, survey=survey, version=version)

    if nooverwrite:
        overwrite_check(opath)

    ddp    = Table.read(fpath)
    ddp.pprint()

    phi_Ms, phis, weights  = lumfn_stepwise(ddp, survey=survey)
    result                 = Table(np.c_[phi_Ms, phis], names=['Ms', 'PHI_STEP'])

    runtime                = calc_runtime(start, 'Writing {}'.format(opath))    

    result.meta['IMMUTABLE'] = 'FALSE'
    
    write_desitable(opath, result)

    ddp                    = Table.read(fpath) 
    ddp['WEIGHT_STEPWISE'] = weights

    runtime = calc_runtime(start, 'Writing {}'.format(fpath))

    write_desitable(fpath, ddp)

    runtime = calc_runtime(start, 'Finished')

    if log:
        sys.stdout.close()
