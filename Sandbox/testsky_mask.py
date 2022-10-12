from desitarget import randoms
from mockfactory.desi import get_brick_pixel_quantities
import fitsio

from mpi4py import MPI

mpicomm = MPI.COMM_WORLD

specrel = 'guadalupe'
prog = 'dark'
specf = fitsio.read('/global/cfs/cdirs/desi/spectro/redux/'+specrel+'/zcatalog/ztile-main-'+prog+'-cumulative.fits')
sel_sky = specf['OBJTYPE'] == 'SKY'
specf = specf[sel_sky]
if mpicomm.rank == 0:
    print(len(specf))
columns = {}
columns['maskbits'] = {'fn': '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{region}/coadd/{brickname:.3s}/{brickname}/legacysurvey-{brickname}-maskbits.fits.fz', 'dtype': 'i2', 'default': 1}

columns[randoms.quantities_at_positions_in_a_brick] = {'drdir': '/global/project/projectdirs/cosmo/data/legacysurvey/dr9/{region}/', 'aprad': 1e-9}  # skip apflux
columns['brickname'] = None
columns['photsys'] = None
catalog = get_brick_pixel_quantities(specf['TARGET_RA'], specf['TARGET_DEC'], columns, mpicomm=mpicomm)
if mpicomm.rank == 0:
    for bit in range(0,13):
        sel = catalog['maskbits'] & 2**bit > 0
        print(bit,len(catalog[sel])/len(catalog))

