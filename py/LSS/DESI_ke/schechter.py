import numpy as np
import scipy.integrate as integrate

from   data.schechters import schechters


def schechter(M, phistar, Mstar, alpha):
    expa         = 10. ** (0.4 * (Mstar - M) * (1. + alpha))
    expb         = np.exp(-10. ** (0.4 * (Mstar - M)))

    return  np.log(10.) * phistar * expa * expb / 2.5

def named_schechter(M, named_type='TMR', zz=None, evolve=False):
    params       = schechters[named_type]
    
    log10phistar = params['log10phistar']
    Mstar        = params['Mstar']
    alpha        = params['alpha'] 

    P            = params['P']
    Q            = params['Q']
    zref         = params['zref']

    if zz == None:
        zz = zref

    zz           = np.array([zz], copy=True)[0]
    
    phistar      = 10. ** log10phistar

    if evolve:
        Mstar       -= Q * (zz - zref)
        phistar     *= 10. ** (0.4 * P * (zz - zref))

    return schechter(M, phistar, Mstar, alpha)
