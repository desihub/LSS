import numpy as np
from astropy.table import Table, join, vstack
import os
import sys
import fitsio
from LSS import common_tools as common

fsse = fitsio.read(
    '/dvs_ro/cfs/cdirs/desicollab/users/zhaoc/spec_catas_dr2_elg/output/ELG_merged_properties_v2.fits.gz')
selbad = fsse['SSE_WITH_MASK'] >= 98.96

bad_tids = fsse['TARGETID'][selbad]

input_dir = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2/nonKP/'
output_dir = '/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2/testsse/'
regl = ['NGC', 'SGC']
nran = 18
for reg in regl:
    print('Processing region: ', reg)
    data = fitsio.read(os.path.join(
        input_dir, f'ELG_LOPnotqso_{reg}_clustering.dat.fits'))
    sel = np.isin(data['TARGETID'], bad_tids)
    data = data[~sel]
    print('removed ', np.sum(sel), ' objects from data')
    common.write_LSS_scratchcp(data, os.path.join(
        output_dir, f'ELG_LOPnotqso_{reg}_clustering.dat.fits'))
    for nr in range(nran):
        print('Processing randoms: ', nr)
        randoms = fitsio.read(os.path.join(
            input_dir, f'ELG_LOPnotqso_{reg}_{nr}_clustering.ran.fits'))
        sel = np.isin(randoms['TARGETID_DATA'], bad_tids)
        randoms = randoms[~sel]
        common.write_LSS_scratchcp(randoms, os.path.join(
            output_dir,  f'ELG_LOPnotqso_{reg}_{nr}_clustering.ran.fits'))
        print('removed ', np.sum(sel), ' objects from randoms')

    print('Done with region: ', reg)
print('Script finished!')
