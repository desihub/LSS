import math
import numpy as np
from scipy import integrate
from matplotlib import pyplot
from scipy import interpolate


def H_z(z, H_0, omega_m, omega_l):
  Hz = H_0*np.sqrt(omega_m*(1+z)**3 + omega_l)

  return Hz

def r_comoving(z, H_0, omega_m, omega_l):
  c= 299792.45 #km/s
  try:
     n=H_0.size
     rcom=np.zeros(n)
     for ii in range(0,n):
        hfact = lambda x: 1.0/H_z(x, H_0[ii], omega_m[ii], omega_l[ii])
        rtemp = integrate.quad(hfact,0,z)
        rcom[ii] = c*rtemp[0]
  except:
     n=1
     hfact = lambda x: 1.0/H_z(x, H_0, omega_m, omega_l)
     rtemp = integrate.quad(hfact,0,z)
     rcom = c*rtemp[0]

  return rcom

def V_comoving(z1, z2, H_0, omega_m, omega_l, angle):
  c= 299792.45 #km/s
  hfact = lambda x: (r_comoving(x,H_0,omega_m,omega_l))**2./H_z(x, H_0, omega_m, omega_l)
  vtemp = integrate.quad(hfact,z1,z2)

  vcom = c*vtemp[0]*angle
  return vcom

def D_A(z, H_0, omega_m, omega_l):
  try:
     n=H_0.size
     DA=np.zeros(n)
     for ii in range(0,n):
        r_c = r_comoving(z, H_0[ii], omega_m[ii], omega_l[ii])
        DA[ii] = r_c/(1.+z)
  except:
     n=1
     r_c = r_comoving(z, H_0, omega_m, omega_l)
     DA  = r_c/(1.+z)

  return DA

def r_s(omega_m,omega_b, H_0):
   h=H_0/100.0
   OM_mh2=omega_m*h*h
   OM_bh2=omega_b*h*h
   s=44.5*np.log(9.83/OM_mh2)/np.sqrt(1+10*np.power(OM_bh2,0.75))
   #print 'Sound Horizon at drag Epoch from Eisenstein and HU 1998=',s
   return s

def HZ(z,Ho,omM,omL):
   Hz=Ho*np.sqrt(omM*np.power(1+z,3)+(1-omM-omL)*np.power(1.0+z,2)+omL )
   return Hz

def D_V(z,Ho,omM,omL):
   c= 299792.45 #km/s
   DA=D_A(z, Ho, omM, omL)
   Hz=HZ(z,Ho,omM,omL)
   DV=np.power(c*z*np.power(1.0+z,2)*np.power(DA,2)/Hz,1.0/3.0)

   return DV

def omLz(z,Ho,omM,omL):
   Hz=HZ(z,Ho,omM,omL)
   Lz=omL*np.power(Ho/Hz,2)
   return Lz

def omMz(z,Ho,omM,omL):
   Hz=HZ(z,Ho,omM,omL)
   Mz=omM*np.power(1.0+z,3)*np.power(Ho/Hz,2)
   return Mz

def gz(z,Ho,omM,omL):
   Mz=omMz(z,Ho,omM,omL)
   Lz=omLz(z,Ho,omM,omL)
   gz=(5.0*Mz/2.0)*(1/(np.power(Mz,4.0/7.0)-Lz+((1+Mz/2.0)*(1+Lz/70.0))))
   return gz

def sigma8z(z,Ho,omM,omL):
   Dz=gz(z,Ho,omM,omL)/(1+z)
   Do=gz(0,Ho,omM,omL)
   sigma8z=Dz/Do
   return sigma8z


def zintp(redshift,H0=67,omM=0.27,interp='',zmax=1.0):
   if(interp==''): #interpolate for redshift to comoving distance
      omL=1.0-omM
      Nintp=1000
      zsamp_intp=np.linspace(0,zmax,Nintp)
      rcomovz=np.zeros(Nintp)
      h=H0/100.0
      for ii in range(0,Nintp): #comoving distance in Mpc/h
         rcomovz[ii]=r_comoving(zsamp_intp[ii],H0,omM,omL)*h
      #create interpolation object
      interp = interpolate.splrep(zsamp_intp,rcomovz, s=0,k=1)

   rco=interpolate.splev(redshift,interp, der=0)

   return rco, interp

def RDZ2XYZ(RA,DEC,redshift,H0,omM,omL,interp=''):
   ngal=RA.size
   XYZ=np.zeros(ngal*3).reshape(ngal,3)

   #interpolate for redshift to comoving distance
   rco, interp=zintp(redshift,H0=H0,omM=omM,interp=interp,zmax=5.0)
   #rco=interpolate.splev(redshift,interp, der=0)
   theta=np.pi*(90-DEC)/180
   phi  =np.pi*RA/180

   XYZ[:,0]=rco*np.sin(theta)*np.cos(phi)
   XYZ[:,1]=rco*np.sin(theta)*np.sin(phi)
   XYZ[:,2]=rco*np.cos(theta)

   return XYZ,interp


if __name__ == "__main__":
 #Fiducial cosmology
 if(0):
   H_0 = 70.
   omega_m = 0.274
   omega_l = 1.0-omega_m
   omega_b = 0.045714

 #best fit without sigma8
 if(0):
   H_0 = 67.92
   omega_m = 0.307
   omega_l = 1.0-omega_m
   omega_b = 0.045714

 #PTHALO mocks:
 if(0):
   H_0=70.0
   omega_m=0.274
   omega_l=1-omega_m
   omega_b=0.0458571
   s8z0=0.80

 #QPM fiducial
 if(0):
   H_0 = 70.0
   omega_m = 0.29
   omega_l = 1.0-omega_m
   omega_b=0.0458571

 #Mock Challenge fiducial
 if(0):
   H_0 = 67.6
   omega_m = 0.31
   omega_l = 1.0-omega_m
   omega_b=0.048142

 #Patchy fiducial
 if(0):
   H_0 = 67.77
   omega_m = 0.307115
   omega_l = 1.0-omega_m
   omega_b=0.048206
   s8z0=0.8288

 #Acacia fiducial
 if(1):
   #OmM=0.31 h=0.676 Ol=0.69 Obh2=0.022
   H_0 = 67.6
   omega_m = 0.31
   omega_l = 1.0-omega_m
   omega_b=0.022/np.power(H_0/100,2)
   s8z0=0.80

 #WMAP cosmology
 if(0):
   H_0 = 69.7
   omega_m = 0.2815
   omega_l = 1.0-omega_m
   omega_b=0.022/np.power(H_0/100,2)
   s8z0=0.82

 #sound horizon
 #r_s(omega_m,omega_b,H_0)

 #z = np.arange(0,0.61,0.01)
 z=np.array([0,0.27,0.4,0.43,0.55,0.7,0.61,1.1])
 nz = np.size(z)
 r = np.zeros(nz)
 if(1):
   print 'z[i],r[i],Hz,DAz,rs,OM_z, np.power(OM_z,0.55),s8z*s8z0, np.power(OM_z,0.55)*s8z*s8z0'
   for i in range(0,nz):
      r[i] = r_comoving(z[i], H_0, omega_m, omega_l)*H_0/100.0
      Hz=H_z(z[i], H_0, omega_m, omega_l),
      DAz=D_A(z[i], H_0, omega_m, omega_l)
      OM_z=omMz(z[i],H_0,omega_m, omega_l)
      s8z=sigma8z(z[i],H_0,omega_m,omega_l)
      rs= r_s(omega_m,omega_b, H_0)
      print z[i],r[i],Hz,DAz,rs,OM_z, np.power(OM_z,0.55),s8z*s8z0, np.power(OM_z,0.55)*s8z*s8z0

 #pyplot.plot(z, r)
 #pyplot.show()
 #pyplot.savefig("com.ps")
 #pyplot.clf()
