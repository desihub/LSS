'''
Copied from Julian Baustista
'''

import numpy as np
import pylab as P
import os
import json
from numba import jit
from scipy.ndimage.filters import gaussian_filter
from scipy.fftpack import fftn, ifftn, fftshift, fftfreq
from CosmoJB import Cosmo
import pyfftw

#def getras(xa,ya):
#	ra = np.zeros(xa.size)
#	for i in range(0,xa.size):
#		x = xa[i]
#		y = ya[i]
#		rai = np.arctan(x/y)*180./np.pi# + 180 
#		if x < 0 and y > 0:
#			rai += 360.
#		if x < 0 and y < 0:
#			rai += 180.
#		if x > 0 and y < 0:
#			rai += 180	
#		ra[i] = rai
#	return ra

class MiniCat:

   def __init__(self, ra, dec, z, we):
        self.ra = ra
        self.dec = dec
        self.Z = z
        self.we = we
        self.size = ra.size 

class Recon:

    def __init__(self, \
                 data_ra, data_dec, data_z, data_we, \
                 rand_ra, rand_dec, rand_z, rand_we, \
                 bias=2.3, f=0.817, smooth=15., nbins=256, \
                 padding=200., opt_box=1, nthreads=1):

        beta = f/bias

        #-- parameters of box
        cosmo = Cosmo(OmegaM=0.31)
        print 'Num bins:', nbins
        print 'Smoothing [Mpc/h]:', smooth

        #-- getting weights
        dat = MiniCat(data_ra, data_dec, data_z, data_we)
        ran = MiniCat(rand_ra, rand_dec, rand_z, rand_we)

        #-- computing cartesian positions
        dat.dist = cosmo.get_comoving_distance(dat.Z)
        ran.dist = cosmo.get_comoving_distance(ran.Z)
        dat.x = dat.dist * np.cos(dat.dec*np.pi/180)*np.cos(dat.ra*np.pi/180)
        dat.y = dat.dist * np.cos(dat.dec*np.pi/180)*np.sin(dat.ra*np.pi/180)
        dat.z = dat.dist * np.sin(dat.dec*np.pi/180) 
        dat.newx = dat.x*1.
        dat.newy = dat.y*1.
        dat.newz = dat.z*1.
        ran.x = ran.dist * np.cos(ran.dec*np.pi/180)*np.cos(ran.ra*np.pi/180)
        ran.y = ran.dist * np.cos(ran.dec*np.pi/180)*np.sin(ran.ra*np.pi/180)
        ran.z = ran.dist * np.sin(ran.dec*np.pi/180) 
        ran.newx = ran.x*1.
        ran.newy = ran.y*1.
        ran.newz = ran.z*1.

        #print 'Randoms min of x, y, z', min(ran.x), min(ran.y), min(ran.z)
        #print 'Randoms max of x, y, z', max(ran.x), max(ran.y), max(ran.z)

        sum_wgal = np.sum(dat.we)
        sum_wran = np.sum(ran.we)
        alpha = sum_wgal/sum_wran
        ran_min = 0.01*sum_wran/ran.size

        self.nbins=nbins
        self.bias=bias
        self.f=f
        self.beta=beta
        self.smooth=smooth
        self.dat = dat
        self.ran = ran
        self.ran_min = ran_min
        self.alpha=alpha
        self.cosmo = cosmo
        self.nthreads = nthreads

        self.compute_box(padding=padding, optimize_box=opt_box)

    def compute_box(self, padding=200., optimize_box=1):
    
        if optimize_box:
            minx = np.min(self.ran.x)
            maxx = np.max(self.ran.x)
            miny = np.min(self.ran.y)
            maxy = np.max(self.ran.y)
            minz = np.min(self.ran.z)
            maxz = np.max(self.ran.z)

            dx = maxx-minx
            dy = maxy-miny
            dz = maxz-minz
            x0 = 0.5*(maxx+minx)
            y0 = 0.5*(maxy+miny) 
            z0 = 0.5*(maxz+minz) 

            box = max([dx, dy, dz])+2*padding
            self.xmin = x0-box/2 
            self.ymin = y0-box/2 
            self.zmin = z0-box/2 
            self.box = box
            self.binsize = box/self.nbins
        else:
            box = self.cosmo.get_comoving_distance(1.05)
            self.xmin=-box
            self.ymin=-box
            self.zmin=-box
            self.box = box*2.
            self.binsize = 2.*box/self.nbins
        
        print 'Box size [Mpc/h]:', self.box
        print 'Bin size [Mpc/h]:', self.binsize

    #@profile
    def iterate(self, iloop, save_wisdom=1):
        dat = self.dat
        ran = self.ran
        smooth = self.smooth
        binsize = self.binsize
        beta = self.beta
        bias = self.bias
        f = self.f
        nbins = self.nbins

        #-- Creating arrays for FFTW
        if iloop==0:
            delta  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
            deltak = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
            psi_x  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
            psi_y  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
            psi_z  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
            #delta = np.zeros((nbins, nbins, nbins), dtype='complex128')
            #deltak= np.zeros((nbins, nbins, nbins), dtype='complex128')
            #psi_x = np.zeros((nbins, nbins, nbins), dtype='complex128')
            #psi_y = np.zeros((nbins, nbins, nbins), dtype='complex128')
            #psi_z = np.zeros((nbins, nbins, nbins), dtype='complex128')
            
            print 'Allocating randoms in cells...'
            deltar = self.allocate_gal_cic(ran)
            print 'Smoothing...'
            deltar = gaussian_filter(deltar, smooth/binsize)

            #-- Initialize FFT objects and load wisdom if available
            wisdomFile = "wisdom."+str(nbins)+"."+str(self.nthreads)               
            if os.path.isfile(wisdomFile) :
                print 'Reading wisdom from ',wisdomFile
                g = open(wisdomFile, 'r')
                wisd=json.load(g)
                pyfftw.import_wisdom(wisd)
                g.close()
            print 'Creating FFTW objects...'
            fft_obj  = pyfftw.FFTW(delta, delta, axes=[0, 1, 2], \
                                   threads=self.nthreads)
            ifft_obj = pyfftw.FFTW(deltak, psi_x, axes=[0, 1, 2], \
                                   threads=self.nthreads, \
                                   direction='FFTW_BACKWARD')
        else:
            delta = self.delta
            deltak = self.deltak
            deltar=self.deltar
            psi_x = self.psi_x
            psi_y = self.psi_y
            psi_z = self.psi_z
            fft_obj = self.fft_obj
            ifft_obj = self.ifft_obj

        #fft_obj = pyfftw.FFTW(delta, delta, threads=self.nthreads, axes=[0, 1, 2])
        #-- Allocate galaxies and randoms to grid with CIC method
        #-- using new positions
        print 'Allocating galaxies in cells...'
        deltag = self.allocate_gal_cic(dat)
        print 'Smoothing...'
        deltag = gaussian_filter(deltag, smooth/binsize)

        print 'Computing fluctuations...'
        delta[:]  = deltag - self.alpha*deltar
        w=np.where(deltar>self.ran_min)
        delta[w] = delta[w]/(self.alpha*deltar[w])
        w2=np.where((deltar<=self.ran_min)) 
        delta[w2] = 0.
        #-- removing this to not change statistics
        #w3=np.where(delta>np.percentile(delta[w].ravel(), 99))
        #delta[w3] = 0.
        del(w)
        del(w2)
        #del(w3)
        del(deltag)

        print 'Fourier transforming delta field...'
        norm_fft = 1.#binsize**3 
        fft_obj(input_array=delta, output_array=delta)
        #delta = pyfftw.builders.fftn(\
        #                delta, axes=[0, 1, 2], \
        #                threads=self.nthreads, overwrite_input=True)()

        #-- delta/k**2 
        k = fftfreq(self.nbins, d=binsize)*2*np.pi
        #-- adding 1e-100 in order to avoid division by zero
        delta /= (k[:, None, None]**2 + k[None, :, None]**2 + k[None, None, :]**2 + 1e-100)
        delta[0, 0, 0] = 0. 

        #-- Estimating the IFFT in Eq. 12 of Burden et al. 2015
        print 'Inverse Fourier transforming to get psi...'
        norm_ifft = 1.#(k[1]-k[0])**3/(2*np.pi)**3*nbins**3

        deltak[:] = delta*-1j*k[:, None, None]/bias
        ifft_obj(input_array=deltak, output_array=psi_x)
        deltak[:] = delta*-1j*k[None, :, None]/bias
        ifft_obj(input_array=deltak, output_array=psi_y)
        deltak[:] = delta*-1j*k[None, None, :]/bias
        ifft_obj(input_array=deltak, output_array=psi_z)

        #psi_x = pyfftw.builders.ifftn(\
        #                delta*-1j*k[:, None, None]/bias, \
        #                axes=[0, 1, 2], \
        #                threads=self.nthreads, overwrite_input=True)().real
        #psi_y = pyfftw.builders.ifftn(\
        #                delta*-1j*k[None, :, None]/bias, \
        #                axes=[0, 1, 2], \
        #                threads=self.nthreads, overwrite_input=True)().real
        #psi_z = pyfftw.builders.ifftn(\
        #                delta*-1j*k[None, None, :]/bias, \
        #                axes=[0, 1, 2], \
        #                threads=self.nthreads, overwrite_input=True)().real
        #psi_x = ifftn(-1j*delta*k[:, None, None]/bias).real*norm_ifft
        #psi_y = ifftn(-1j*delta*k[None, :, None]/bias).real*norm_ifft
        #psi_z = ifftn(-1j*delta*k[None, None, :]/bias).real*norm_ifft

        #-- removing RSD from galaxies
        shift_x, shift_y, shift_z =  \
                self.get_shift(dat, psi_x.real, psi_y.real, psi_z.real, \
                               use_newpos=True)
        print "Few displacement values: "
        for i in range(10):
            print shift_x[i], shift_y[i], shift_z[i], dat.x[i]

        #-- for first loop need to approximately remove RSD component 
        #-- from Psi to speed up calculation
        #-- first loop so want this on original positions (cp), 
        #-- not final ones (np) - doesn't actualy matter
        if iloop==0:
            psi_dot_rhat = (shift_x*dat.x + \
                            shift_y*dat.y + \
                            shift_z*dat.z ) /dat.dist
            shift_x-= beta/(1+beta)*psi_dot_rhat*dat.x/dat.dist
            shift_y-= beta/(1+beta)*psi_dot_rhat*dat.y/dat.dist
            shift_z-= beta/(1+beta)*psi_dot_rhat*dat.z/dat.dist
        #-- remove RSD from original positions (cp) of 
        #-- galaxies to give new positions (np)
        #-- these positions are then used in next determination of Psi, 
        #-- assumed to not have RSD.
        #-- the iterative procued then uses the new positions as 
        #-- if they'd been read in from the start
        psi_dot_rhat = (shift_x*dat.x+shift_y*dat.y+shift_z*dat.z)/dat.dist
        dat.newx = dat.x + f*psi_dot_rhat*dat.x/dat.dist 
        dat.newy = dat.y + f*psi_dot_rhat*dat.y/dat.dist 
        dat.newz = dat.z + f*psi_dot_rhat*dat.z/dat.dist 

        self.deltar = deltar
        self.delta = delta
        self.deltak = deltak
        self.psi_x = psi_x
        self.psi_y = psi_y
        self.psi_z = psi_z
        self.fft_obj = fft_obj
        self.ifft_obj = ifft_obj



        #-- save wisdom
        wisdomFile = "wisdom."+str(nbins)+"."+str(self.nthreads)               
        if iloop==0 and save_wisdom and not os.path.isfile(wisdomFile):
            wisd=pyfftw.export_wisdom()
            f = open(wisdomFile, 'w')
            json.dump(wisd,f)
            f.close()
            print 'Wisdom saved at', wisdomFile

    def apply_shifts(self, verbose=1):
        
        for c in [self.dat, self.ran]:
            shift_x, shift_y, shift_z =  \
                self.get_shift(c, \
                    self.psi_x.real, self.psi_y.real, self.psi_z.real, \
                    use_newpos=True) #-- should be True here, thanks to Seshadri Nadathur
            c.newx += shift_x 
            c.newy += shift_y 
            c.newz += shift_z

    def summary(self):

        dat = self.dat
        sx = dat.newx-dat.x
        sy = dat.newy-dat.y
        sz = dat.newz-dat.z
        print 'Shifts stats for each dimension: std  16th  84th  min  max'
        for s in [sx, sy, sz]:
            print np.std(s), np.percentile(s, 16), np.percentile(s, 84), \
                    np.min(s), np.max(s)


    def allocate_gal_cic(self, c):
        xmin=self.xmin
        ymin=self.ymin
        zmin=self.zmin
        binsize=self.binsize
        nbins=self.nbins

        xpos = (c.newx-xmin)/binsize
        ypos = (c.newy-ymin)/binsize
        zpos = (c.newz-zmin)/binsize

        i = xpos.astype(int)
        j = ypos.astype(int)
        k = zpos.astype(int)

        ddx = xpos-i
        ddy = ypos-j
        ddz = zpos-k

        delta = np.zeros((nbins, nbins, nbins))
        edges = [np.linspace(0, nbins, nbins+1), \
                 np.linspace(0, nbins, nbins+1), \
                 np.linspace(0, nbins, nbins+1)]

        for ii in range(2):
            for jj in range(2):
                for kk in range(2):
                    pos = np.array([i+ii, j+jj, k+kk]).transpose()
                    weight = ( ((1-ddx)+ii*(-1+2*ddx))*\
                               ((1-ddy)+jj*(-1+2*ddy))*\
                               ((1-ddz)+kk*(-1+2*ddz)) ) *c.we
                    delta_t, edges = np.histogramdd(pos, \
                                     bins=edges, weights=weight) 
                    delta+=delta_t

        return delta

    def get_shift(self, c, f_x, f_y, f_z, use_newpos=False):

        xmin = self.xmin
        ymin = self.ymin
        zmin = self.zmin
        binsize = self.binsize

        if use_newpos:
            xpos = (c.newx-xmin)/binsize
            ypos = (c.newy-ymin)/binsize
            zpos = (c.newz-zmin)/binsize
        else: 
            xpos = (c.x-xmin)/binsize
            ypos = (c.y-ymin)/binsize
            zpos = (c.z-zmin)/binsize

        i = xpos.astype(int)
        j = ypos.astype(int)
        k = zpos.astype(int)

        ddx = xpos-i
        ddy = ypos-j
        ddz = zpos-k

        shift_x = np.zeros(c.size)
        shift_y = np.zeros(c.size)
        shift_z = np.zeros(c.size)

        for ii in range(2):
            for jj in range(2):
                for kk in range(2):
                    weight = ( ((1-ddx)+ii*(-1+2*ddx))*\
                               ((1-ddy)+jj*(-1+2*ddy))*\
                               ((1-ddz)+kk*(-1+2*ddz))  )
                    pos = (i+ii, j+jj, k+kk)
                    shift_x += f_x[pos]*weight
                    shift_y += f_y[pos]*weight
                    shift_z += f_z[pos]*weight

        return shift_x, shift_y, shift_z

    def export_cart(self, root):
        
        lines = ['%f %f %f %f'%(xx, yy, zz, ww) for (xx, yy, zz, ww) \
               in zip(self.dat.x, self.dat.y, self.dat.z, self.dat.we)]
        fout = open(root+'.dat.txt', 'w') 
        fout.write('\n'.join(lines))
        fout.close()
        lines = ['%f %f %f %f'%(xx, yy, zz, ww) for (xx, yy, zz, ww) \
               in zip(self.ran.x, self.ran.y, self.ran.z, self.ran.we)]
        fout = open(root+'.ran.txt', 'w') 
        fout.write('\n'.join(lines))
        fout.close()


    def cart_to_radecz(self, x, y, z):

        dist = np.sqrt(x**2+y**2+z**2)
        dec = 90 - np.degrees(np.arccos(z / dist))
        ra = np.degrees(np.arctan2(y, x))
        ra[ra < 0] += 360
        redshift = self.cosmo.get_redshift(dist)
        return ra, dec, redshift


    #def cart_to_radecz(self, x, y, z):

        #dist = np.sqrt(x**2+y**2+z**2)
        #dec = np.arcsin(z/dist)*180./np.pi
        #ra = getras(x,y)
        #ra = np.arctan(x/y)*180./np.pi# + 180 
        #if x < 0 and y > 0:
        #	ra += 180.
        #if x < 0 and y < 0:
        #	ra += 180.
        #if x > 0 and y < 0:
        #	ra += 360. 	
        #redshift = self.cosmo.get_redshift(dist)
        #return ra, dec, redshift

    def get_new_radecz(self, c):

        return self.cart_to_radecz(c.newx, c.newy, c.newz) 




