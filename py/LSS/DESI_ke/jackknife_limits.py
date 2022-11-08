import numpy             as np
import matplotlib.pyplot as plt
from   collections       import OrderedDict


jk_limits = {'JK0': {'ra_min': 129.,  'ra_max': 133.,  'dec_min': -2., 'dec_max': 3.},
             'JK1': {'ra_min': 133.,  'ra_max': 137.,  'dec_min': -2., 'dec_max': 3.},
             'JK2': {'ra_min': 137.,  'ra_max': 141.,  'dec_min': -2., 'dec_max': 3.},
             'JK3': {'ra_min': 174.,  'ra_max': 178.,  'dec_min': -3., 'dec_max': 2.},
             'JK4': {'ra_min': 178.,  'ra_max': 182.,  'dec_min': -3., 'dec_max': 2.},
             'JK5': {'ra_min': 182.,  'ra_max': 186.,  'dec_min': -3., 'dec_max': 2.},
             'JK6': {'ra_min': 211.5, 'ra_max': 215.5, 'dec_min': -2., 'dec_max': 3.},
             'JK7': {'ra_min': 215.5, 'ra_max': 219.5, 'dec_min': -2., 'dec_max': 3.},
             'JK8': {'ra_min': 219.5, 'ra_max': 223.5, 'dec_min': -2., 'dec_max': 3.}}

def set_jackknife(ras, decs, limits=None, debug=True):
    result       = np.array(['None'] * len(ras), dtype=str)

    if limits == None:
        limits   = jk_limits

    for strip in limits.keys():
        ra_min   = limits[strip]['ra_min']
        ra_max   = limits[strip]['ra_max']

        dec_min  = limits[strip]['dec_min']
        dec_max  = limits[strip]['dec_max']

        in_strip = (ras >= ra_min) & (ras <= ra_max) & (decs >= dec_min) & (decs <= dec_max)

        result[in_strip] = strip

        if debug:
            print(strip, limits[strip])

    return result

def plot_jackknife(dat):
    fig, ax = plt.subplots(1,1, figsize=(10,10))

    for idx in np.unique(dat['JK']):
        sub = dat[dat['JK'] == idx]

        ax.scatter(sub['RA'], sub['DEC'], s=0.25, label=idx)

    ax.set_xlabel('RA [deg.]')
    ax.set_ylabel('DEC [deg.]')

    ax.legend(frameon=True, ncol=6)

    return fig, ax

def solve_jackknife(rand, ndiv=4):
    '''
    Splits up dat and rand into jackknife areas based on (ra, dec) in (ndiv x ndiv) chunks.
    '''           
    njack         = ndiv * ndiv
    jk_volfrac    = (njack - 1.) / njack 

    dpercentile   = 100. / ndiv
    percentiles   = np.arange(dpercentile, 100. + dpercentile, dpercentile)
    
    jk            = 0
    limits        = OrderedDict()
 
    for ra_per in percentiles:
        # Given a vector V of length N, the q-th percentile of V is the q-th ranked value in a sorted copy of V. 
        # https://docs.scipy.org/doc/numpy-1.9.2/reference/generated/numpy.percentile.html
        rahigh    = np.percentile(rand[f'RANDOM_RA'], ra_per)
        ralow     = np.percentile(rand[f'RANDOM_RA'], ra_per - dpercentile)

        print('{:.6f}\t{:.6f}'.format(ralow, rahigh))

        for dec_per in percentiles:
            isin    = (rand['RANDOM_RA'] >= ralow) & (rand['RANDOM_RA'] <= rahigh)

            dechigh = np.percentile(rand[f'RANDOM_DEC'][isin], dec_per)
            declow  = np.percentile(rand[f'RANDOM_DEC'][isin], dec_per - dpercentile)

            print('\t{:.6f}\t{:.6f}'.format(declow, dechigh))

            limits[f'JK{jk}'] = {'ra_min': float(ralow), 'ra_max': float(rahigh), 'dec_min': float(declow), 'dec_max': float(dechigh)}

            jk     += 1

    jks = set_jackknife(rand['RANDOM_RA'], rand['RANDOM_DEC'], limits=limits)

    return  njack, jk_volfrac, limits, jks

if __name__ == '__main__':
    import pylab as pl

    from   astropy.table import Table
    from   findfile      import findfile


    field   = 'G9'
    version = 'v2'

    dpath   = findfile(ftype='gold', dryrun=False, survey='gama', version=version)
    rpath   = findfile(ftype='randoms', dryrun=False, field=field, survey='gama', version=version)
    
    dat     = Table.read(dpath)
    dat     = dat[dat['FIELD'] == field] 

    rand    = Table.read(rpath) 
    
    dat, rand = set_jackknife(dat, rand, ndiv=4)

    print(np.unique(dat['JK'].data))

    fig, ax = plot_jackknife(dat) 
    
    fig.savefig('fig.pdf')
