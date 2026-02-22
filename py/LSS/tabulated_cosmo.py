import os
import numpy as np


_dir_tabulated = os.path.join(os.path.dirname(__file__), 'data')


class CosmologyError(Exception):

    """Exception related to cosmology."""


class TabulatedDESI(object):
    """
    Class to load tabulated z->E(z) and z->comoving_radial_distance(z) relations within DESI fiducial cosmology
    (in LSS/data/desi_fiducial_cosmology.dat) and perform the (linear) interpolations at any z.

    >>> cosmo = TabulatedDESI()
    >>> distance = cosmo.comoving_radial_distance([0.1, 0.2])
    >>> efunc = cosmo.efunc(0.3)

    The cosmology is defined in https://github.com/abacusorg/AbacusSummit/blob/master/Cosmologies/abacus_cosm000/CLASS.ini
    and the tabulated file was obtained using https://github.com/adematti/cosmoprimo/blob/main/cosmoprimo/fiducial.py.

    Note
    ----
    Redshift interpolation range is [0, 100].
    """
    _filename = os.path.join(_dir_tabulated, 'desi_fiducial_cosmology.dat')

    def __init__(self):
        self._z, self._efunc, self._comoving_radial_distance = np.loadtxt(self._filename, comments='#', usecols=None, unpack=True)

    def efunc(self, z):
        r"""Return :math:`E(z)`, where the Hubble parameter is defined as :math:`H(z) = H_{0} E(z)`, unitless."""
        z = np.asarray(z)
        mask = (z < self._z[0]) | (z > self._z[-1])
        if mask.any(): raise CosmologyError('Input z outside of tabulated range.')
        return np.interp(z, self._z, self._efunc, left=None, right=None)

    def comoving_radial_distance(self, z):
        r"""Return comoving radial distance, in :math:`\mathrm{Mpc}/h`."""
        z = np.asarray(z)
        mask = (z < self._z[0]) | (z > self._z[-1])
        if mask.any(): raise CosmologyError('Input z outside of tabulated range.')
        return np.interp(z, self._z, self._comoving_radial_distance, left=None, right=None)


if __name__ == '__main__':

    cosmo = TabulatedDESI()
    distance = cosmo.comoving_radial_distance([0.1, 0.2])
    efunc = cosmo.efunc(0.3)
    assert np.allclose(distance, [292.58423977, 570.41981713], rtol=1e-9, atol=0)
    assert np.allclose(efunc, 1.1736407657972752, rtol=1e-9, atol=0)