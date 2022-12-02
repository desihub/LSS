import fitsio
import numpy as np

from   astropy.table import Table


def fetch_lss(pprint=False, sort=False, survey='sv3', release='fuji'):
    # https://desi.lbl.gov/trac/wiki/ClusteringWG/LSScat/SV3/version2.1/
    #
    # /global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/LSScats/2.1/

    survey     = survey.upper()

    # latest sv3
    versions   = {'SV3-everest': 2.1, 'SV3-fuji': 3}
    version    = versions['{}-{}'.format(survey, release)]
        
    # ('RA', 'DEC', 'TARGETID', 'Z', 'NTILE', 'TILES', 'rosette_number', 'rosette_r', 'FRACZ_TILELOCID', 'BITWEIGHTS', 'PROB_OBS', 'WEIGHT_ZFAIL', 'WEIGHT', 'flux_r_dered', 'NZ')
    clustering = Table.read(f'/global/cfs/cdirs/desi/survey/catalogs/{survey}/LSS/{release}/LSScats/{version}/BGS_BRIGHT_clustering.dat.fits')

    # NTILE, TILES, TILELOCIDS, LOCATION_ASSIGNED, TILELOCID_ASSIGNED, sort, COMP_TILE, rosette_number, rosette_r, FRACZ_TILELOCID, BITWEIGHTS, PROB_OBS
    full       = Table.read(f'/global/cfs/cdirs/desi/survey/catalogs/{survey}/LSS/{release}/LSScats/{version}/BGS_BRIGHT_full.dat.fits')

    if pprint:        
        full_cols = np.array(full.dtype.names).astype(str)

        print(full_cols)
        print(clustering.dtype.names)
        
        clustering.pprint()
        full.pprint()

    if sort:
        clustering.sort('TARGETID')
        full.sort('TARGETID')
        
    return  clustering, full


if __name__ == '__main__':
    fetch_lss(pprint=True)
    
