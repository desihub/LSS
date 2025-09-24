'''
Copied from Julian Bautista
'''

import numpy as N

class Cosmo:

    def __init__(self, OmegaM=0.31, h=0.676):
        print 'Initializing cosmology with OmegaM = %.2f'%OmegaM
        c = 299792458. #m/s
        OmegaL = 1.-OmegaM
        ztab = N.linspace(0., 4., 10000)
        E_z = N.sqrt(OmegaL + OmegaM*(1+ztab)**3)
        rtab = N.zeros(ztab.size)
        for i in range(1, ztab.size):
            rtab[i] = rtab[i-1] + c*1e-3 * (1/E_z[i-1]+1/E_z[i])/2. * (ztab[i]-ztab[i-1]) / 100.

        self.h = h
        self.c = c
        self.OmegaM = OmegaM
        self.OmegaL = OmegaL
        self.ztab = ztab
        self.rtab = rtab 

    #-- comoving distance in Mpc/h
    def get_comoving_distance(self, z):
        return N.interp(z, self.ztab, self.rtab)

    def get_redshift(self, r):
        return N.interp(r, self.rtab, self.ztab)


    #-- comoving spherical volume between zmin and zman in (Mpc/h)**3
    def shell_vol(self, zmin, zmax):
        rmin = self.get_comoving_distance(zmin)
        rmax = self.get_comoving_distance(zmax)
        return 4*N.pi/3.*(rmax**3-rmin**3)

    def get_box_size(self, ra, dec, zmin=0.5, zmax=1.0):

        dmin = get_comoving_distance(zmin)
        dmax = get_comoving_distance(zmax)

        theta = (-dec+90)*N.pi/180.
        phi = (ra)*N.pi/180

        xmin = dmin * N.sin(theta)*N.cos(phi)
        ymin = dmin * N.sin(theta)*N.sin(phi)
        zmin = dmin * N.cos(theta)
        xmax = dmax * N.sin(theta)*N.cos(phi)
        ymax = dmax * N.sin(theta)*N.sin(phi)
        zmax = dmax * N.cos(theta)

        for pair in [[xmin, xmax], [ymin, ymax], [zmin, zmax]]:
            print abs(N.array(pair).min() - N.array(pair).max())
