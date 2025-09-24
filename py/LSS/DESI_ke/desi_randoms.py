import os
import numpy         as     np

from   bitmask       import lumfn_mask
from   astropy.table import Table, unique, vstack
from   ros_tools     import tile2rosette, calc_rosr, ros_limits


def desi_randoms(ros, nrealz=4, oversample=8, dryrun=False):
    assert  'NERSC_HOST' in os.environ.keys()

    # Randoms uniform on the sphere with density 2500 per sq. deg., available to an assigned fiber.      
    root             = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'

    rpaths           = ['{}/random{}/rancomb_brightwdup_Alltiles.fits'.format(root, idx) for idx in range(nrealz)]
    rpaths          *= oversample

    rand             = vstack([Table.read(xx) for xx in rpaths])
    rand.meta        = {}
  
    # TODO:  Check TARGETID is a unique identifier, or bug.  If not, use RA. 
    rand             = unique(rand, keys='TARGETID')
    rand['ROS']      = tile2rosette(rand['TILEID'])

    rand['ROS_DIST'] = 1.e99

    rand             = rand[rand['ROS'] == ros]
    
    for rosn in np.unique(rand['ROS']):
        isin = (rand['ROS'].data == rosn)

        new_dist = calc_rosr(rosn, rand['RA'][isin], rand['DEC'][isin])
    
        rand['ROS_DIST'][isin] = np.minimum(rand['ROS_DIST'][isin], new_dist)

    # rand.pprint()

    limits                    = ros_limits(dryrun)

    hi_comp             = (rand['ROS_DIST'].data > limits[0]) & (rand['ROS_DIST'].data < limits[1])
    rand['IN_D8LUMFN']  = ~hi_comp * lumfn_mask.DESI_HICOMP

    rand                = rand[rand['IN_D8LUMFN'].data == 0]
    
    rand.rename_column('RA',  'RANDOM_RA')
    rand.rename_column('DEC', 'RANDOM_DEC')

    # Must come after dryrun.
    rand.meta['AREA'] = len(rand) / 2500. / nrealz
    rand.meta['NRAND'] = len(rand)
    rand.meta['NREALZ'] = nrealz
    rand.meta['IMMUTABLE'] = 'TRUE'
        
    return  rand
