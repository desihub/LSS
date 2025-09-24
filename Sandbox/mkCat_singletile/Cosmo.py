# Robert J. Brunner Aug 2005
# Added to by Adam D. Myers Sep 2005+
# Useful cosmological calculations
# Changed by Ashley Ross to allow for w, 2010
#from math import sqrt,cos,log10,modf
from math import *
#from scipy.integrate import romberg
from romberg import rom
#from scipy.special import gamma as gamfunc
#from gamma import gamfunc
#from scipy import interpolate
#from histo import histo
from time import sleep
#import biggles
#import Numeric

def Gamma(m,mb,h):
	k = .1
	ag = 1. -.328*log(431.*m*h**2.)*mb/m+.38*log(22.3*m*h**2.)*(mb/m)**2.
	s =44.5*log(9.83/m*h**2.)/sqrt(1.+10.*(mb*h**2.)**.75)
	#print s
	return m*h*(ag+(1.-ag)/(1.+(.43*k*s)**4.))
	#return m*h*ag


def testdom(om_min=.2,om_max=.4,om_step=.001,zmin=.1,zmax=1.,zstep=.001):
	fo = open('dvom'+str(zmax)+'.dat','w')
	om = om_min
	while om < om_max:
		d = distance(omega=om,lamda=1.-om)
		z = zmin
		while z < zmax:
			dis = d.dc(z)
			disnc = d.intevnc(z)
			fo.write(str(dis)+' '+str(om)+' '+str(disnc)+' '+str(z)+'\n')
			z += zstep
		om += om_step
	fo.close()
	return True

def chifunc4ameoba(list,a,b,em=1.,ap1=1.):
	f = open('dvom1.0.dat').readlines()
	ssq = 0
	for i in range(0,len(f)):
		ln = f[i].split()
		rd = float(ln[0])
		ssq += (rd-(float(ln[1])**a*float(ln[2])*float(ln[3])**b))**2.
	return ssq

def chifunc4ameobaA(list,a,b,em=1.,ap1=1.):
	f = open('dvom1.0.dat').readlines()
	ssq = 0
	for i in range(0,len(f)):
		ln = f[i].split()
		rd = float(ln[0])
		ssq += (rd-(float(ln[1])**a*float(ln[2])*b*float(ln[3])**.1+50))**2.
	return ssq


def alph2DA(z,alph,err,om,lam,h=.7,w=-1.):
	df = distance(omega=om,lamda=lam,h=h,w=w)
	DAf = df.da(z)
	rs = 153.19
	DA = alph*DAf
	errA = err*DA
	return DA,errA

def alph(z,omf,hf,obhhf,om,h,obhh):
	dfid = distance(omf,1.-omf,h=hf,obhh=obhhf)
	Dvrs_fid = dfid.dV(z)/dfid.rs
	dn = distance(om,1.-om,h=h,obhh=obhh)
	Dvrs = dn.dV(z)/dn.rs
	Hfrs = dfid.Hz(z)*dfid.rs
	Hrs = dn.Hz(z)*dn.rs
	DAfrs = dfid.da(z)/dfid.rs
	DArs = dn.da(z)/dn.rs
	return Dvrs/Dvrs_fid,Hrs/Hfrs,DArs/DAfrs,dn.dV(z)/dfid.dV(z),dn.dc(z)/dfid.dc(z),(dfid.rs*hf)/(dn.rs*h)

def alphtheta(z,omf,hf,obhhf,om,h,obhh):
	dfid = distance(omf,1.-omf,h=hf,obhh=obhhf)
	dn = distance(om,1.-om,h=h,obhh=obhh)
	thbaofid = dfid.rs/dfid.dc(z)
	thbaon = dn.rs/dn.dc(z)
	print(dfid.rs/dn.rs)
	return thbaon/thbaofid



class distance:
# I know I spelled lambda wrong but do YOU know WHY?
	def __init__(self,omega=0.3,lamda=0.7,h=1,w=-1.,gam=.557,obhh=.0224):
		self.gamma = gam
		self.w = w
		self.fa = -3.*(1.+w)
		self.obhh = obhh
		self.om = float(omega)
		self.ol = float(lamda)
		self.h0 = float(h)
		self.c = 2997.9 # speed of light in km/centi-second
		self.K = (h**2.)*(omega+lamda-1.) # Curvature parameter
		self.g0 = 2.5*self.om/(self.om**0.571428571-self.ol+((1+(0.5*self.om))*(1+(0.014285714*self.ol))))
   # above is approximate Local linear growth factor
	#below are factors to obtain physical bao scale		
		self.tfac = 2.725/2.7
		tp =2.725/2.73
		omhh = self.om*h*h
		
		#zeq = 25000.*omhh/self.tfac**4.		
		zeq = 23900*omhh/tp**4.-1.
		#keq = .0746*omhh/self.tfac/self.tfac
		keq = 1./self.c*sqrt(2.*omhh*zeq)
		b1 = 0.313*(omhh)**(-.419)*(1.+.607*omhh**.674)
		b2 = .238*omhh**.223
		zd = 1291.*omhh**.251/(1.+.659*omhh**.828)*(1.+b1*obhh**b2)
		Req = self.RR(zeq)
		Rd = self.RR(zd)
		self.rs = 2./(3.*keq)*sqrt(6./Req)*log((sqrt(1.+Rd)+sqrt(Rd+Req))/(1.+sqrt(Req)))
		#print self.rs
	
	def RR(self,z):
		return 31.5*self.obhh*self.tfac**-4./(.001*z)
		
	def omz(self,z): # Omega_m (as a function of z FLAT COSMOLOGY!)
		return self.om*(1+z)**3./(self.om*(1.+z)**3.+self.ol*(1.+z)**self.fa)

	def olz(self,z): # Omega_Lamda (as a function of z FLAT COSMO!)
		return self.ol/(self.om*(1.+z)**3.+self.ol*(1.+z)**self.fa)*(1.+z)**self.fa
	def Hz(self,z):  # Hubble parameter h(z) (FLAT COSMOLOGY!!!!).
		return self.h0*sqrt(self.om*(1.+z)**3.+self.ol*(1.+z)**self.fa)

	def cHz(self,z):  # speed of light divided by Hubble parameter h(z) (FLAT COSMOLOGY!!!!).
		return self.c/self.Hz(z)

	def evolution(self, z):
		return 1./(sqrt(self.om*(1.+z)**3.+self.ol*(1.+z)**self.fa))

	def dV(self,z):
		#spherically averaged distance quantity
		return (z*self.cHz(z)*self.dc(z)**2.)**(1./3.)
	def da(self, z):   # Angular diameter distance from now to z
		return self.dc(z)/(1.+z)
	def dl(self, z):   # Luminosity distance from now to z
		return self.dc(z)*(1.+z)
	def dc(self, z): # Comoving distance from now to z
		return ((self.c/self.h0)*rom(0,z,self.evolution))    
#        
	def dc2z(self,x):
		zup = 1.
		zdown = .001
		xup = self.dc(zup)
		xdown = self.dc(zdown)
		zt = (zup+zdown)/2.
		xt = self.dc(zt)
		while abs(xt - x) > .001:
			while xdown > x:
				zdown = zdown - .1
				xdown = self.dc(zdown)
				print( zdown,xdown,x)
			if abs(x - xdown) < abs(xup - x):
				zup = zt
				xup = xt
				zt = (zdown+zt)/2.
				xt = self.dc(zt)
			else:
				zdown = zt
				xdown = xt
				zt = (zup+zt)/2.
				xt = self.dc(zt)
			print( xt,zt,zdown,zup,xdown,xup)
		return zt

	def mkD(self,zb=.01):
		Dl = []
		oldD = 1.
		olda = 1.
		for i in range(0,2000000):
			z = .001*i
			a = 1./(1.+z)
			da = a - olda
			aev = 1./(1.+z+.0005)
			D = oldD + da*oldD*(self.om/(self.om+self.ol*aev**3.))**self.gamma/aev
			Dl.append(D)
			oldD = D
			olda = a
		return Dl


    
	def pvolfunc(self,z):
		ans = self.dc(z)**2.*(self.c/self.h0)*self.evolution(z)*(1.+z)**-3.
		#print ans
		return ans
    
	def pvol(self,z1,z2): #proper volume element (dsolidangle) in shell between z1 and z2
		return rom(z1,z2,self.pvolfunc)

	def covolfunc(self,z):
		#ans = 4.*pi*self.dc(z)**2.*(self.c/self.h0)*self.evolution(z)
		ans = 4.*pi*self.dc(z)**2.*self.c/self.Hz(z)
		#print ans
		return ans
    
	def covol(self,z1,z2): #full-sky comoving volume in shell between z1 and z2
		return rom(z1,z2,self.covolfunc)

		
		
    
	def ngobsfunc(self,z,dndzfile='nzDR5VL20.5_mrs'):
    	
		dndzcomfile=dndzfile+'_ref.dat'
		ind = int(z*1000)
		#print dndzfile,dndzcomfile
		f = open(dndzfile+'.dat').readlines()
		dz = .001
		fcom = open(dndzcomfile).readlines()
		nz = float(f[ind].split()[1])

		if nz == 0:
			return 0,0
		indcom = ind-100 #if using file that starts from z = .1
		#indcom = ind
		if indcom<0:
			return 'problem with your zrange'
		nzcom = float(fcom[indcom].split()[1])
		if nzcom == 0:
			#zweight = 0
			zfac = 0
			return 0,0
		else:
			#zweight = nzcom**2.
			zweight = nzcom
			#zweight = 1.
			#if z > .4:
			#	zweight = (nz/float(fcom[400].split()[1]))**2.
			#zfac = nz/nzcom
    		
		#volel = self.dc(z)**2.*(self.c/self.h0)*self.evolution(z)*(1.+z)**-3.*dz*zfac
		#volel = 4.*pi*self.dc(z)**2.*(self.c/self.h0)*self.evolution(z)*(1.+z)**-3.*dz
		volel = 4.*pi*self.dc(z)**2.*(self.c/self.h0)*self.evolution(z)*dz
		nel = nz*zweight/volel
		#nel = nz/volel
		#print z,nel,nz,volel,zfac
		return nel,zweight
		#return nel
    	
	def ngobs(self,z1,z2,dzfile='nzDR5VL20.5_mrs'):
		obdeg = 5406.5453944537448
		spheredeg = 360.*360./pi
		sphfac = obdeg/spheredeg
		zweighttot = 0
		ntot = 0
		print(dzfile)
		it = 1000.*(z2-z1)
		for i in range(0,it):
			z = z1 + .001*i
			ev = self.ngobsfunc(z,dndzfile=dzfile)
			#print ev
			#print z,ev
			ntot += ev[0]

			#ntot += ev[0]
			zweighttot += ev[1]
		ng = ntot/zweighttot/sphfac
		print( ng,ntot,zweighttot)
		#ng = ntot*spheredeg/obdeg/(4.*pi)/float(it)
		return ng

	def nzvolref(self,zref=.336,dndzfile='nzDR5VL.3.4_mrs'):
		indref = int(zref*1000)
		#print dndzfile
		f = open(dndzfile+'.dat').readlines()
		dz = .001

		nzref = float(f[indref].split()[1])    	
	
		volref = self.dc(zref)**2.*(self.c/self.h0)*self.evolution(zref)*(1.+zref)**-3.*dz
		nref = nzref/volref
		fileout = dndzfile+'_ref.dat'
		fo = open(fileout,'w')
		print( dndzfile,fileout)
		for i in range(100,len(f)): #use if you want to start from z=.1
		#for i in range(1,len(f)):
			z = i/1000.+.001
			vol = self.dc(z)**2.*(self.c/self.h0)*self.evolution(z)*(1.+z)**-3.*dz
			nVL = nzref*vol/volref
			dndzf = float(f[i].split()[1])/nVL
			fo.write(str(z)+' '+str(dndzf)+'\n')
    
	def sepc(self,z1,z2,theta): # Comoving distance between two objects
# with redshifts z1, z2 and angular separation of theta degrees in
# any cosmology
		d2r = 0.017453292519943295   # factor to convert degree to radian
		s1 = self.dc(z1)
		s2 = self.dc(z2)
		s1sq = s1**2.
		s2sq = s2**2.
		costh = cos(d2rtheta)
		return sqrt(s1sq+s2sq-(self.K*s1sq*s2sq*(1+(costh**2.)))-(2.*s1*s2*sqrt(1-(self.K*s1sq))*sqrt(1-(self.K*s2sq))*costh))

	def sepcflat(self,z1,z2,theta): # Comoving distance between objects
# with redshifts z1, z2 and angular separation of theta in a flat
# cosmology
		s1 = self.dc(z1)
		s2 = self.dc(z2)
		s1sq = s1**2.
		s2sq = s2**2.
		costh = cos(theta)
		return sqrt(s1sq+s2sq-(2.*s1*s2*costh))

	def scflat(self,z1,z2): # Comoving distance between objects
# with redshifts z1, z2 and no angular separation in a flat
# cosmology
		s1 = self.dc(z1)
		s2 = self.dc(z2)
		s1sq = s1**2.
		s2sq = s2**2.
		return sqrt(s1sq+s2sq-(2.*s1*s2))

	def saflat(self,z1,z2): # Angular diameter distance between objects
# with redshifts z1, z2 and no angular separation in a flat
# cosmology
		return self.scflat(z1,z2)/(1.+z2)

	def sep(self,s1,s2,costh): # Comoving distance between objects
# with comoving distances s1, s2 and a cosine of their angular
# separation of costh, in a flat cosmology. (This is just the cosine
# rule!)
		s1sq = s1**2.
		s2sq = s2**2.
		return sqrt(s1sq+s2sq-(2.*s1*s2*costh))

	def dm(self,z):
# Distance Modulus by redshift in a given cosmology
		return 5.*log10(self.dl(z)) + 25.

	def Kcorr(self,z,alph=0):
# Simple k-correction in given cosmology for given spectral slope
		return 2.5*(1+alph)*log10(1.+z)
	def Kcorr2(self,z,alph=-0.45,K0=-0.42):
# Two parameter k-correction a la Wisotzki 2000A&A...353..861
		return K0+self.Kcorr(z,alph)
	def KcorrLRG(self,z):
		return 4.42*z-2.459
	def AbsMag(self,mag,z):
# Given an apparent magnitude and redshift return absolute magnitude
# with 'standard' K-correction in a given cosmology 
		return mag - self.Kcorr2(z) - self.dm(z)
    
	def AbsMag_nk(self,mag,z):
		return mag - self.dm(z)
    
	def AbsMagLRG(self,mag,z):
# Given an apparent magnitude and redshift return absolute magnitude
# with 'standard' K-correction in a given cosmology 
		return mag - self.KcorrLRG(z) - self.dm(z)
	def AppMag(self,mag,z):
# Given an absolute magnitude and redshift return apparent magnitude
# with 'standard' K-correction in a given cosmology 
		return mag + self.Kcorr2(z) + self.dm(z)

	def gam(self,z):
#The exponent for f = omz**gamma RSD models
		return 3./(5.-self.w/(1.-self.w)) + 3./125.*(1.-self.w)*(1.-1.5*self.w)/(1.-1.2*self.w)**3.*(1.-self.omz(z))

	def Dg(self,z):
#The approximate linear growth rate for the w/gamma
		bA = -.28/(self.w+.08)-.3
		return 2.5*self.omz(z)/(1.+z)*(self.omz(z)**self.gam(z)-self.olz(z)+(1.+.5*self.omz(z))*(1.+bA*self.olz(z)))**-1.

	def Dgn(self,z):
		return self.Dg(z)/self.Dg(0)

	def D(self,z):
# The linear growth factor in the class cosmology as a function of
# redshift, as approximated by Carrol, Press and Turner (1992)
# in Equation 29	
		return self.g(z)/(1+z)


	def g(self,z):
# The suppresion of the linear growth factor in the class
# cosmology as a function of redshift, as given by
# Carrol, Press and Turner (1992)
		return (2.5*self.omz(z)/(self.omz(z)**0.571428571-self.olz(z)+((1+(0.5*self.omz(z)))*(1+(0.014285714*self.olz(z))))))/self.g0
	def Daccurate(self,z):
# Linear growth factor as a function of redshift in the class
# cosmology as given by Carrol, Press and Turner (1992)
# by integrating their equation 28 in full, using da/dTau as
# given by their equation 9. 1e-15 prevents division by zero.
		a = 1./(1.+z)
		return 2.5*self.om*self.dadt(a)*(1.+z)*rom(1e-15,a,self.dadtint)

	def DaccurateRenorm(self,z):
# Linear growth factor as a function of redshift in the class
# cosmology as given by Carrol, Press and Turner (1992)
# by integrating their equation 28 in full, using da/dTau as
# given by their equation 9. 1e-15 prevents division by zero.
# Renormalized to give a=1, which is self-evidently true
		return self.Daccurate(z)/self.Daccurate(0)
    
	def Dsimp(self,z):
		a = 1./(1.+z)
		return exp(log(a)*self.omz(z)**self.gamma)
    
	def dadtint(self,a):
# equation 9 of Carrol, Press and Turner (1992) in the class
# cosmology, with an extra integration power of -3.  This can
# be integrated to give the linear growth factor
		return (1.+(self.om*((1./a)-1.))+(self.ol*(a**2.-1)))**-1.5

	def dadt(self,a):
# equation 9 of Carrol, Press and Turner (1992) in the class
# cosmology
		return (1.+(self.om*((1./a)-1.))+(self.ol*(a**2.-1)))**0.5

	def b0(self,b,z1,z2):
		#given the bias at one redshift, return the passively evolved bias at the other
		return 1.+ (b-1.)*self.DaccurateRenorm(z2)/self.DaccurateRenorm(z1)

	def epfac(self,b0,z1,z2):
    	#return the fraction by which the clustering strength should change for passively evolving objects
		return (b0/(b0-1.+self.DaccurateRenorm(z2)/self.DaccurateRenorm(z1)))**2.

	def evolution_nocos(self,z):
		return 1./(sqrt((1.+z)**3.+(1.+z)**self.fa))
	
	def intevnc(self,z):
		return ((self.c/self.h0)*rom(0,z,self.evolution_nocos))

class Limber:
# Limber's Equation given a cosmology, a dN/dz and a value of the
# real space correlation function slope.  dN/dz should be passed
# as a Numeric array interpolated spline, which can be generated
# from data using interpolate.splev/.splrep. I (Adam Myers) have
# a program called histo.histo that'll do this. It requires a file
# that contains a single column of redshifts, which, in the nominal
# case, should be called dNdZ.dat. This file will be binned at a
# resolution of 'bins' and interpolated. If sigma > 0.0 this dNdZ
# file is broadened by a Gaussian at either end of the redshift
# range a la Brunner, Szalay & Connolly (2000).
    def __init__(self,omega=0.3,lamda=0.7,h=1,gamma=1.8,dNdZfile="dNdZ.dat",display=False,tolerance=1.48e-8,epsilon=-3,bins=25,sigma=0.0):
# Note that if display is true, the interpolation of dNdZ will
# be plotted to an x11 window, if biggles is installed
# Tolerance is the accuracy of Romberg integrations
# Epsilon is the Epsilon form of the clustering evolution, which is
# used in EpsLimber
        self.om = float(omega)
        self.ol = float(lamda)
        self.h0 = float(h)
        self.c = 2997.9 # speed of light in km/centi-second
        self.dNdZ = histo(file=dNdZfile,label='z',nox=not display,bins=bins,sigma=sigma)
        self.gam = float(gamma)
        up = float(self.dNdZ[0][len(self.dNdZ[0])-1])
        down = float(self.dNdZ[0][0])
        self.max = up+(0.5*(up-down)/(bins-1.))
        self.min = down-(0.5*(up-down)/(bins-1.))
        if self.min <1e-4:
            self.min = 1e-4
            print( "WARNING:********************************************")
            print("Numerical Integration of Limber's Equation diverges")
            print("at z = 0 - you *might* choose a larger value of minz")
            print("or write a better integration program!  For now:\n")
            print("resetting minimum z to 0.0001")
            print( "WARNING:********************************************")
        self.Hgamma = gamfunc(0.5)*gamfunc(0.5*(float(gamma)-1))/gamfunc(0.5*float(gamma))
        self.tol = float(tolerance)
        self.eps = float(epsilon)
# Yields the gamma function prefix on Limber's Equation.  Don't
# confuse the gamma power-law slope with the Gamma Function itself!
    def LinearLimber(self):
# Variant of the integration of Limber's Equation in the linear
# regime. There's no epsilon factor - we assume the linear regime
# value and evolve it with z, using Carrol, Press and Turner (1992)
        check = self.CheckNorm()
        print( "%s %s%s" %('Normalization accurate to',check,'%'))
        check = 1.-(0.01*check)
        print( 'Integrating........')
        sleep(5)
        return self.Hgamma*romberg(self.TopFunc,(self.min,self.max),self.tol)/self.c/check
    def TopFunc(self,z): # numerator of Limbers Eqn*self.c
        N=interpolate.splev(z,self.dNdZ)
        Hub = distance(self.om,self.ol,self.h0).Hz(z)
        lgro = distance(self.om,self.ol,self.h0).D(z)
        com = distance(self.om,self.ol,self.h0).dc(z)
        return N**2.*Hub*lgro**2.*(com**(1.-self.gam))
    def BotFunc(self,z): # sqrt of denominator of Limber's Eqn
        N=interpolate.splev(z,self.dNdZ)
        return N
    def CheckNorm(self):
# Check that the denominator is normalized to unity (it should
# be because we call histo.py asking for a normalized interpolation
# This is particularly useful in determining a value of zmax!
        print( '***********************************************************')
        print( 'Returning % integrated error on normalization of dN/dZ.....')
        return (1.-romberg(self.BotFunc,(self.min,self.max),self.tol))*100.
    def EpsLimber(self):
# Variant of the integration of Limber's Equation using the
# Peebles epsilon formalism
        check = self.CheckNorm()
        print( "%s %s%s" %('Normalization accurate to',check,'%'))
        if check > 0.1 :
            print( 'Normalization fairly inaccurate......')
            print( 'You could try running with a different value of zmax')
        print( 'Integrating........')
        sleep(5)
        return self.Hgamma*romberg(self.TopFuncEps,(self.min,self.max),self.tol)/self.c
    def TopFuncEps(self,z): # numerator of Limbers Eqn*self.c
        N=interpolate.splev(z,self.dNdZ)
        Hub = distance(self.om,self.ol,self.h0).Hz(z)
        lgro = (1.+z)**(-3.-self.eps+self.gam)
        com = distance(self.om,self.ol,self.h0).dc(z)
        return N**2.*Hub*lgro*(com**(1.-self.gam))
    def r0fromAEps(self,A):
# Given a value of A, the amplitude of the projected correlation
# function in arcmin^(gamma-1), return the real-space scale-length
        arcmin2rad=0.00029088820866572158
# First convert the amplitude to radians
        return 10.**(log10((A*arcmin2rad**(self.gam-1.))/self.EpsLimber())/self.gam)
    def r0fromALinear(self,A):
# Given a value of A, the amplitude of the projected correlation
# function in arcmin^(gamma-1), return the real-space scale-length
        arcmin2rad=0.00029088820866572158
# First convert the amplitude to radians
        return 10.**(log10((A*arcmin2rad**(self.gam-1.))/self.LinearLimber())/self.gam)
    def biasinfofromALinear(self,A,err):
# Given a value of A, the amplitude of the projected correlation
# function in arcmin^(gamma-1), and associated error return
# and print to screen in TeX format the input amplitude&error
# real-space scale-length&error, and the bias parameter&error
# Assumes CDM scale-length of 5Mpc/h at z=0 and that self.gam
# is a good approximation to the underlying matter spectrum
        arcmin2rad=0.00029088820866572158
# First convert the amplitude to radians
        Apcerr = (float(err))/A
        r_0 = 10.**(log10((A*arcmin2rad**(self.gam-1.))/self.LinearLimber())/self.gam)
        r_0err = Apcerr*r_0/self.gam
        b = (r_0/5.)**(0.5*self.gam)
        berr = 0.5*Apcerr*b
# How much to round to (assume 1 sf more than error input)
# PROBLEM IF r_o IS BELOW ZERO WHICH IT never SHOULD BE, IN THIS
# CASE!!!
        stringerr = str(err)
        flag = False
        while stringerr[0] == '0' or stringerr[0] == '.':
            stringerr = stringerr.lstrip('0')
            stringerr = stringerr.lstrip('.')
            flag = True
        sf = len(stringerr)-1  # sig figs (implicit +1 'cos of the
# decimal point)
        if flag: sf += 1

        print( sf, 'significant figures')

        rounder = self.round2sf(err,sf)
        string1 = "%s %s %s %s"%(str(round(A,rounder)),'$\pm$',str(round(err,rounder)),'&')
        rounder = self.round2sf(r_0err,sf)
        string2 = "%s %s %s %s"%(str(round(r_0,rounder)),'$\pm$',str(round(r_0err,rounder)),'&')
        rounder = self.round2sf(berr,sf)
        string3 = "%s %s %s %s"%(str(round(b,rounder)),'$\pm$',str(round(berr,rounder)),'&')

        print ('A=   ', string1)
        print ('r_0=     ', string2)
        print ('b=   ' , string3)

        return string1, string2, string3

    def r0fromtheta0Eps(self,th0):
# Given a value of th0, the scale-length of the projected
# correlation function in arcmin, give the real-space scale-length
        arcmin2rad=0.00029088820866572158
# First convert the amplitude to radians
        return 10.**(log10(((th0*arcmin2rad)**(self.gam-1.))/self.EpsLimber())/self.gam)
    def r0fromtheta0Linear(self,th0):
# Given a value of th0, the scale-length of the projected
# correlation function in arcmin, give the real-space scale-length
        arcmin2rad=0.00029088820866572158
# First convert the amplitude to radians
        return 10.**(log10(((th0*arcmin2rad)**(self.gam-1.))/self.LinearLimber())/self.gam)
    def PlotFunc(self,delay=5): # Plot the cosmological functions that
# contribute to Limber's Equation.  Delay is how many secs to wait
# between plots
        step = (self.max-self.min)/10000.
        zval = []
        Hubval = []
        lgroval = []
        comval = []
        HubEdS = []
        lgroEdS = []
        comEdS = []

        for i in range(10000):
            zed = i*step
            zval.append(zed)
            Hubval.append(distance(self.om,self.ol,self.h0).Hz(zed))
            lgroval.append(distance(self.om,self.ol,self.h0).D(zed))
            comval.append(distance(self.om,self.ol,self.h0).dc(zed))
            HubEdS.append(distance(1,0,1).Hz(zed))
            lgroEdS.append(distance(1,0,1).D(zed))
            comEdS.append(distance(1,0,1).dc(zed))

        titstring = '($\Omega_m$='+str(self.om)+' $\Omega_\Lambda$='+str(self.ol)+' h='+str(self.h0)+')'
        message = 'Dotted line is EdS cosmology'
        author = "Adam D. Myers, UIUC"

        plot = bigles.FramedPlot()
        plot.xlabel = 'z'
        plot.ylabel = 'Hubble Parameter/$100kms^{-1}Mpc^{-1}, h(z)'
        plot.title_style["color"] = "red"
        plot.title = titstring 
        plot.add(biggles.Curve(Numeric.array(zval),Numeric.array(Hubval),color="red"))
        plot.add(biggles.Curve(Numeric.array(zval),Numeric.array(HubEdS),type="dotted"))
        plot.add( biggles.PlotLabel(.23, .95, author,size=2.3) )
        plot.add( biggles.PlotLabel(.72, .05, message,size=2.3) )
        plot.show()
        plot.write_eps("HubPar.eps")
        sleep(delay)

        plot = biggles.FramedPlot()
        plot.xlabel = 'z'
        plot.ylabel = 'Linear Growth Factor, D(z)'
        plot.title_style["color"] = "green"
        plot.title = titstring 
        plot.add(biggles.Curve(Numeric.array(zval),Numeric.array(lgroval),color="green"))
        plot.add(biggles.Curve(Numeric.array(zval),Numeric.array(lgroEdS),type="dotted"))
        plot.add( biggles.PlotLabel(.8, .95, author,size=2.3) )
        plot.add( biggles.PlotLabel(.25, .05, message,size=2.3) )
        plot.show()
        plot.write_eps("LinearGrowth.eps")
        sleep(delay)

        plot = biggles.FramedPlot()
        plot.xlabel = 'z'
        plot.ylabel = 'Comoving Distance, $\Chi$(z), (Mpc)'
        plot.title_style["color"] = "blue"
        plot.title = titstring 
        plot.add(biggles.Curve(Numeric.array(zval),Numeric.array(comval),color="blue"))
        plot.add(biggles.Curve(Numeric.array(zval),Numeric.array(comEdS),type="dotted"))
        plot.add( biggles.PlotLabel(.23, .95, author,size=2.3) )
        plot.add( biggles.PlotLabel(.72, .05, message,size=2.3) )
        plot.show()
        plot.write_eps("ComovDist.eps")

    def round2sf(self,err,sf):
#
# Given an parameter, and a number of significant figures, work out
# What the round command integer is to get to that number of sf
#
        stringerr = str(err)
        a = stringerr.split(".")
        if a[0] != '0': rounder = len(a[0])
        else:
            stringerr = a[1]
            rounder = 0
            while stringerr[0] == '0':
                rounder -= 1
                stringerr = stringerr.lstrip('0')
        rounder-=sf
        rounder*=-1

        return rounder

class QSO:
# models for QSOs based on links between QSO bias, Black Hole Mass
# Dark Matter Halo Mass, Accretion Efficiency etc.
    def __init__(self,omega=0.3,lamda=0.7,h=1,MDMHtimesh=3.0e12):
# Here omega,lamda and h are the usual cosmological parameters
# MDMHtimesh is the Mass of a dark matter halo in units of
# Msolar/h.  At some point we'll want to write a function to work
# this out from the QSO bias using Sheth,Mo,Tormen(2001)
#
        self.om = float(omega)
        self.ol = float(lamda)
        self.h0 = float(h)
        self.MDMH = MDMHtimesh/float(h)
        self.c = 2997.9 # speed of light in km/centi-second
    def MBHWyLo(self,z,instance=1):
# Determine the Black Hole Mass from the Dark Matter Halo mass
# Units are same as those of the DMH mass, so typically M_solar
# Assuming the relationship of Wyithe and Loeb (2004)
# This assumes that bulge velocity dispersion and MBH is constant
# With redshift (sigma-MBH)
# instance is either:
#       1 means we think the DMH profile is SIS
#       2 means we think the DMH profile is NFW
#       3 means we think the DMH profile is S02 (Seljak 2002)
# z is the redshift

       eps = 10.**-5.1
       if instance == 2: eps *= 3.7
       if instance == 3: eps *= 25

# This is the ol' concentration parameter of NFW fame the preamble
# constant is, interestingly, 18*pi^2
 
       a = 177.65287921960845
       omzed = distance(self.om,self.ol).omz(z)
       DELC = a  +  82.*(omzed-1)  -  39.*(omzed-1)**2.

# Now piece it all together to work out MBH

       BHM = (self.MDMH*1e-12)**(2./3.)
       BHM *= (DELC*self.om/(a*omzed))**(5./6.)
       BHM *= (1+z)**(5./2.)
       BHM *= eps*self.MDMH
       ex = modf(log10(BHM))[1]
       print( "Mass is",BHM*10.**-ex,'x 10 **',ex)
       return BHM
    def MBHFerr(self,z,instance=1):
# Determine the Black Hole Mass from the Dark Matter Halo mass
# Units are same as those of the DMH mass, so typically M_solar
# Assuming the relationship of Ferrarese (2002) ApJ 578, 90
# This assumes that bulge velocity dispersion and MBH is NOT
# constant with redshift
# instance is either:
#       1 means we think the DMH profile is SIS
#       2 means we think the DMH profile is NFW
#       3 means we think the DMH profile is S02 (Seljak 2002)
# z is the redshift

       if instance == 1:
           BHM=1e8*0.027*(self.MDMH/1e12)**1.82
           ex = modf(log10(BHM))[1]
           print( "Mass is",BHM*10.**-ex,'x 10 **',ex)
           return BHM
       if instance == 2:
           BHM=1e8*0.1*(self.MDMH/1e12)**1.65
           ex = modf(log10(BHM))[1]
           print( "Mass is",BHM*10.**-ex,'x 10 **',ex)
           return BHM
       if instance == 3:
           BHM=1e8*0.67*(self.MDMH/1e12)**1.82
           ex = modf(log10(BHM))[1]
           print( "Mass is",BHM*10.**-ex,'x 10 **',ex)
           return BHM
    def LEdd(self,BHM):
# Eddington Limiting Luminosity in Watts given the Black Hole Mass
# in solar masses
        L = (BHM/1e8)*10.**39.1
        ex = modf(log10(L))[1]
        print( "Eddington Luminosity is",L*10.**-ex,'x 10 **',ex)
        return L
    def BolfromM(self,M,alph=-0.5,reflam=2500,measlam=4669,Lbolmult=[6.3,3.2,9.4]):
# Given a bolometric luminosity, return the absolute magnitude
# in the band centred at a wavelength of measlam, following the
# schema of Richards et al. 2005 and Elvis 1994
#
# The conversions from M to L_bol are:
#
# first convert from K-corrected M to L:
# log(L) = -0.4*(M-5+48.6)+log(4pi)+2log(3.086*10^18)
# L is in erg/s/Hz/cm/cm see Oke&Gunn (1983)
#
# The rest of the example if for g at a reference wavelength of 2500
# Then convert from L_g (at 4669) to L_2500, using a slope of alpha
# L_2500 = log(L_g)-(alpha*log(4669/2500))
#
# Then convert to L_Bol using (Elvis 1994)
# L_bol = 6.3*L_2500
#
# The luminosity is then in erg/s/Hz so convert to Watts by 
# multiplying by 10^7

        measlam=float(measlam)
        reflam=float(reflam)
 
        logL = (-0.4*(M+43.6))+log10(12.566370614359172)+(2.*log10(3.086e18))
        L = 10.**(logL+(alph*log10(measlam/reflam)))
        L *= 1.e7

        return Lbolmult[0]*L, Lbolmult[1]*L, Lbolmult[2]*L

    def MfromBol(self,L,alph=-0.5,reflam=2500,measlam=4669,Lbolmult=6.3):
# Essentially, reverse of BolFromL

        measlam=float(measlam)
        reflam=float(reflam)

        L /= (Lbolmult*1.e7)
        logL = log10(L)-(alph*log10(measlam/reflam))
        M = 2.5*log10(12.566370614359172)+5.*log10(3.086e18)-2.5*logL-43.6
        return M
    
