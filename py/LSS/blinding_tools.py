#tools to help apply blinding routines to LSS catalogs
#requires cosmology calculations from cosmoprimo

import fitsio
from astropy.table import Table

from LSS.common_tools import write_LSS

from cosmoprimo.fiducial import DESI
from cosmoprimo.utils import DistanceToRedshift
from cosmoprimo import Cosmology

cosmo_fid = DESI()
dis_fid = cosmo_fid.comoving_radial_distance


def apply_zshift_DE(data,out_file,w0=-1,wa=0,zcol='Z'):
    #data is table of LSS catalog info
    #out_file is the full path for where to write the output
    #w0 is the equation of state today to use to shift the redshifts
    #wa is the change in w0 w.r.t. the scale factor
    #zcol is the column name
    #data = Table(fitsio.read(in_file))
    cosmo_d = Cosmology(h=cosmo_fid['h'],Omega_m=cosmo_fid['Omega_m'],w0_fld=w0,wa_fld=wa)
    dis_val = dis_fid(data[zcol])
    d2z = DistanceToRedshift(cosmo_d.comoving_radial_distance)
    z_shift = d2z(dis_val)
    data['Z'] = z_shift
    write_LSS(data,out_file,comments=None)

