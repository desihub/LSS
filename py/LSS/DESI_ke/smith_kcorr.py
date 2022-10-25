import os
import numpy as np
import matplotlib.pyplot as plt
import os

from   scipy.interpolate import interp1d
from   pkg_resources     import resource_filename

raw_dir = os.environ['CODE_ROOT'] + '/data/'        

class GAMA_KCorrection(object):
    def __init__(self, band, kind="linear"):
        """
        Colour-dependent polynomial fit to the GAMA K-correction (Fig. 13 of Smith+17), 
        used to convert between SDSS r-band Petrosian apparent magnitudes, and rest 
        frame absolute manigutues at z_ref = 0.1
        
        Args:
            k_corr_file: file of polynomial coefficients for each colour bin
            z0: reference redshift. Default value is z0=0.1
            kind: type of interpolation between colour bins,
                  e.g. "linear", "cubic". Default is "linear"
        """
        k_corr_file = raw_dir + '/ajs_kcorr_{}band_z01.dat'.format(band.lower())
        
        # read file of parameters of polynomial fit to k-correction
        # polynomial k-correction is of the form
        # A*(z-z0)^4 + B*(z-z0)^3 + C*(z-z0)^2 + D*(z-z0) + E 
        col_min, col_max, A, B, C, D, E, col_med = \
            np.loadtxt(k_corr_file, unpack=True)
    
        self.z0 = 0.1             # reference redshift
        self.nbins = len(col_min) # number of colour bins in file
        
        self.colour_min = np.min(col_med)
        self.colour_max = np.max(col_med)
        self.colour_med = col_med

        # functions for interpolating polynomial coefficients in rest-frame color.
        self.__A_interpolator = self.__initialize_parameter_interpolator(A, col_med, kind=kind)
        self.__B_interpolator = self.__initialize_parameter_interpolator(B, col_med, kind=kind)
        self.__C_interpolator = self.__initialize_parameter_interpolator(C, col_med, kind=kind)
        self.__D_interpolator = self.__initialize_parameter_interpolator(D, col_med, kind=kind)
        self.__E = E[0]

        # Linear extrapolation for z > 0.5
        self.__X_interpolator = lambda x: None
        self.__Y_interpolator = lambda x: None
        self.__X_interpolator, self.__Y_interpolator = self.__initialize_line_interpolators() 
   
    def __initialize_parameter_interpolator(self, parameter, median_colour, kind="linear"):
        # returns function for interpolating polynomial coefficients, as a function of colour
        return interp1d(median_colour, parameter, kind=kind, fill_value="extrapolate")
    
    def __initialize_line_interpolators(self):
        # linear coefficients for z>0.5
        X = np.zeros(self.nbins)
        Y = np.zeros(self.nbins)
        
        # find X, Y at each colour
        redshift = np.array([0.48,0.5])
        arr_ones = np.ones(len(redshift))
        for i in range(self.nbins):
            k = self.k(redshift, arr_ones*self.colour_med[i])
            X[i] = (k[1]-k[0]) / (redshift[1]-redshift[0])
            Y[i] = k[0] - X[i]*redshift[0]
        
        X_interpolator = interp1d(self.colour_med, X, kind='linear', fill_value="extrapolate")
        Y_interpolator = interp1d(self.colour_med, Y, kind='linear', fill_value="extrapolate")
        
        return X_interpolator, Y_interpolator

    def __A(self, colour):
        # coefficient of the z**4 term
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__A_interpolator(colour_clipped)

    def __B(self, colour):
        # coefficient of the z**3 term
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__B_interpolator(colour_clipped)

    def __C(self, colour):
        # coefficient of the z**2 term
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__C_interpolator(colour_clipped)

    def __D(self, colour):
        # coefficient of the z**1 term
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__D_interpolator(colour_clipped)

    def __X(self, colour):
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__X_interpolator(colour_clipped)

    def __Y(self, colour):
        colour_clipped = np.clip(colour, self.colour_min, self.colour_max)
        return self.__Y_interpolator(colour_clipped)


    def k(self, redshift, restframe_colour, median=False):
        """
        Polynomial fit to the GAMA K-correction for z<0.5
        The K-correction is extrapolated linearly for z>0.5

        Args:
            redshift: array of redshifts
            colour:   array of ^0.1(g-r) colour
        Returns:
            array of K-corrections
        """
        K   = np.zeros(len(redshift))
        idx = redshift <= 0.5
        
        if median:
            restframe_colour = np.copy(restframe_colour)
            
            # Fig. 13 of https://arxiv.org/pdf/1701.06581.pdf
            restframe_colour = 0.603 * np.ones_like(restframe_colour)

        K[idx] = self.__A(restframe_colour[idx])*(redshift[idx]-self.z0)**4 + \
                 self.__B(restframe_colour[idx])*(redshift[idx]-self.z0)**3 + \
                 self.__C(restframe_colour[idx])*(redshift[idx]-self.z0)**2 + \
                 self.__D(restframe_colour[idx])*(redshift[idx]-self.z0) + self.__E 

        idx = redshift > 0.5
        
        K[idx] = self.__X(restframe_colour[idx])*redshift[idx] + self.__Y(restframe_colour[idx])
        
        return  K    

    def k_nonnative_zref(self, refz, redshift, restframe_colour, median=False):
        refzs = refz * np.ones_like(redshift)
        
        return  self.k(redshift, restframe_colour, median=median) - self.k(refzs, restframe_colour, median=median) - 2.5 * np.log10(1. + refz)

    def rest_gmr_index(self, rest_gmr, kcoeff=False):
        bins = np.array([-100., 0.18, 0.35, 0.52, 0.69, 0.86, 1.03, 100.])
        idx  = np.digitize(rest_gmr, bins)
        '''
        if kcoeff==True:
            for i in enumerate(rest_gmr):
                ddict = {i:{col_med, A[0], B[0], C[0], D[0]}}
        '''
        return idx

class GAMA_KCorrection_color():
    def __init__(self):
        self.kRcorr = GAMA_KCorrection(band='R')
        self.kGcorr = GAMA_KCorrection(band='G')

    def obs_gmr(self, rest_gmr):        
        return  rest_gmr + self.kRcorr.k(z, rest_gmr) - self.kGcorr.k(z, rest_gmr)

    def rest_gmr_nonnative(self, native_rest_gmr):
        refzs = np.zeros_like(native_rest_gmr)
        
        return  native_rest_gmr + self.kGcorr.k(refzs, native_rest_gmr) - self.kRcorr.k(refzs, native_rest_gmr) 

    
def test_plots(axes):
    kcorr_r = GAMA_KCorrection(band='R')
    kcorr_g = GAMA_KCorrection(band='G')

    z    = np.arange(-0.01,0.601,0.01)
    cols = 0.130634, 0.298124, 0.443336, 0.603434, 0.784644, 0.933226, 1.06731

    # make r-band k-correction plot                                                                                                                      
    for c in cols:
        col = np.ones(len(z)) * c
        k = kcorr_r.k(z, col)
        axes[0].plot(z, k, label=r"$^{0.1}(g-r)_\mathrm{med}=%.3f$"%c)

    axes[0].set_xlabel(r"$z$")
    axes[0].set_ylabel(r"$^{0.1}K_r(z)$")
    axes[0].set_xlim(0,0.6)
    axes[0].set_ylim(-0.6,1)
    axes[0].legend(loc="upper left").draw_frame(False)

    # make g-band k-correction plot                                                                                                                      
    for c in cols:
        col = np.ones(len(z)) * c
        k = kcorr_g.k(z, col)
        axes[1].plot(z, k, label=r"$^{0.1}(g-r)_\mathrm{med}=%.3f$"%c)

    axes[1].set_xlabel(r"$z$")
    axes[1].set_ylabel(r"$^{0.1}K_g(z)$")
    axes[1].set_xlim(-0.01,0.6)
    axes[1].set_ylim(-0.4,1.4)
    axes[1].legend(loc="upper left").draw_frame(False)
    axes[0].set_ylabel(r"$^{0.1}K_r(z)$")
    axes[0].set_xlim(0,0.6)
    axes[0].set_ylim(-0.6,1)
    axes[0].legend(loc="upper left").draw_frame(False)

    # make g-band k-correction plot
    for c in cols:
        col = np.ones(len(z)) * c
        k = kcorr_g.k(z, col)
        axes[1].plot(z, k, label=r"$^{0.1}(g-r)_\mathrm{med}=%.3f$"%c)

    axes[1].set_xlabel(r"$z$")
    axes[1].set_ylabel(r"$^{0.1}K_g(z)$")
    axes[1].set_xlim(-0.01,0.6)
    axes[1].set_ylim(-0.4,1.4)

def test_nonnative_plots(axes, zref):    
    kcorr_r = GAMA_KCorrection(band='R')
    kcorr_g = GAMA_KCorrection(band='G')

    z    = np.arange(-0.01,0.601,0.01)
    cols = 0.130634, 0.298124, 0.443336, 0.603434, 0.784644, 0.933226, 1.06731

    # make r-band k-correction plot                                                                                                                     
    for c in cols:
        col = np.ones(len(z)) * c
        k = kcorr_r.k_nonnative_zref(zref, z, col)
        axes[0].plot(z, k, label=r"$^{0.0}(g-r)_\mathrm{med}=%.3f$"%c)

    axes[0].set_xlabel(r"$z$")
    axes[0].set_ylabel(r"$^{0.0}K_r(z)$")
    axes[0].set_xlim(0,0.6)
    axes[0].set_ylim(-0.6,1)
    axes[0].legend(loc="upper left").draw_frame(False)

    # make g-band k-correction plot                                                                                                                     
    for c in cols:
        col = np.ones(len(z)) * c
        k = kcorr_g.k_nonnative_zref(zref, z, col)
        axes[1].plot(z, k, label=r"$^{0.0}(g-r)_\mathrm{med}=%.3f$"%c)

    axes[1].set_xlabel(r"$z$")
    axes[1].set_ylabel(r"$^{0.0}K_g(z)$")
    axes[1].set_xlim(-0.01,0.6)
    axes[1].set_ylim(-0.4,1.4)
    
