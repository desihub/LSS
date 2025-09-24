# f.pyx: numpy arrays -> extern from "fc.h"
# 3 steps:
# cython calc.pyx  -> f.c
# link: python f-setup.py build_ext --inplace  -> f.so, a dynamic library
# py test-f.py: import f gets f.so, f.fpy below calls fc()

from __future__ import division
import numpy as np
cimport numpy as np
from matplotlib import pyplot as plt

cdef extern from "./corr2dc.h" nogil:
    void corr2d(double *xc_2d, double *p1, long np1, double *p2, long np2, double *rlim,
            int nbins0, int nbins1, int nhocells, double *blen,  
	    double *posmin, int samp, int njn, int pbc, int los, int interactive);


def corr2dpy(np.ndarray[np.double_t,ndim=2] p1,
        np.ndarray[np.double_t,ndim=2] p2, np.ndarray[np.double_t,ndim=1] rlim, 
	nbins0,nbins1, nhocells, 
	np.ndarray[np.double_t,ndim=1] blen, 
	np.ndarray[np.double_t,ndim=1] posmin, samp,njn,pbc,los,interactive):
    """
    Calculate the 2d correlation function.

    Parameters
    ----------
    p1: (np1,6/7) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 1.
    p2: (np1,6/7) array of double
        Particle x, y, z, vx, vy, vz, (weight) of Group 2.
    rlim[4]: double, two dimensional correlation function sampling
        minimum,Maximum of first and second co-ordinate range to look at.
    nbins0,nbins1: int, two dimensional correlation function sampling
        Number of bins in the estimator.
    nhocells: int
        Number of cells.
    blen[3]: double, 3d box size or the extent of data
        Box length.
    posmin[3]: double, minimum value of x,y,z co-ordinate
        Box begining.
    samp : int
        To decide the sampling of the correlation function
    njn : int 0 , >0 (number of jacknife regions)
        To decide whether the number of jacknife regions
    interactive:0 can be used to print information

    Returns
    -------
    xc: (nbins,nbins) array of long
        2d correlation function.

    """
    cdef np.ndarray[np.double_t,ndim=2] xc = np.zeros((nbins0, nbins1), 
	                                dtype=np.float64)

    cdef np.ndarray[np.double_t,ndim=3] xc_jn = np.zeros((nbins0, nbins1,njn+1),
	                                          dtype=np.float64)
	                               
    nattr = 5 if njn > 0 else 4
    assert p1.shape[1] == p2.shape[1] == nattr
    np1, np2 = p1.shape[0], p2.shape[0]

    if(njn==0):
       corr2d(<double*> xc.data, <double*> p1.data, np1, <double*> p2.data, np2,
            <double*> rlim.data, nbins0,nbins1, nhocells, <double*> blen.data, 
	    <double*> posmin.data,samp, njn, pbc, los,interactive) 

       return xc
    else:
       corr2d(<double*> xc_jn.data, <double*> p1.data, np1, <double*> p2.data, np2,
            <double*> rlim.data, nbins0,nbins1, nhocells, <double*> blen.data, 
	    <double*> posmin.data,samp, njn, pbc, los,interactive) 

       return xc_jn
