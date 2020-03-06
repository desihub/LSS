minisvdir = '/global/cscratch1/sd/ajross/miniSV/'
datadir = minisvdir+'LSScats/'
dirpcadw = '/global/cscratch1/sd/ajross/pcadw/'
dirpc = '/global/cscratch1/sd/ajross/paircounts/'
import fitsio
import numpy as np

from matplotlib import pyplot as plt

def P2(mu):
	return .5*(3.*mu**2.-1.)
	
def P4(mu):
	return .125*(35.*mu**4.-30.*mu**2.+3.)

def P6(mu):
	return 1./16.*(231.*mu**6.-315.*mu**4.+105.*mu**2.-5.)

def P8(mu):
	return 1./128.*(6435.*mu**8.-12012.*mu**6.+6930.*mu**4.-1260.*mu**2.+35.)


om = 0.31

def createSourcesrd_ad(sample,tile,date,zmin=.5,zmax=1.1):
	'''
	prepare minisv files for paircounts
	'''
	#from healpix import healpix,radec2thphi
	from random import random
	from Cosmo import distance
	d = distance(om,1.-om) #cosmology assumed in final BOSS analyses, make sure this conforms to current
	#h = healpix()
	file = sample+tile+'_'+date

	fd = fitsio.read(datadir+file+'_clustering.dat.fits')
	wz = (fd['Z'] > zmin) & (fd['Z'] < zmax)
	zw = '_zm'+str(zmin)+'zx'+str(zmax)
	zl = fd[wz]['Z']

	cdl = np.zeros(len(zl))
	for i in range(0,len(zl)):
		cdl[i] = d.dc(zl[i])	
	sral = np.sin(np.radians(fd[wz]['RA']))
	cral = np.cos(np.radians(fd[wz]['RA']))
	sdecl = np.sin(np.radians(fd[wz]['DEC']))
	cdecl = np.cos(np.radians(fd[wz]['DEC']))
	wl = fd[wz]['WEIGHT']
	print(str(len(cdl))+' data objects going out for paircounts')
	fdo = open(dirpcadw+ 'g'+file+zw+'pcadw.dat','w')
	for i in range(0,len(cdl)):
		fdo.write(str(sral[i])+' '+str(cral[i])+' '+str(sdecl[i])+' '+str(cdecl[i])+' '+str(cdl[i])+' '+str(wl[i])+'\n')

	fdo.close()

	fr = fitsio.read(datadir+file+'_clustering.ran.fits')
	wz = (fr['Z'] > zmin) & (fr['Z'] < zmax)

	zl = fr[wz]['Z']

	cdl = np.zeros(len(zl))
	for i in range(0,len(zl)):
		cdl[i] = d.dc(zl[i])	
	sral = np.sin(np.radians(fr[wz]['RA']))
	cral = np.cos(np.radians(fr[wz]['RA']))
	sdecl = np.sin(np.radians(fr[wz]['DEC']))
	cdecl = np.cos(np.radians(fr[wz]['DEC']))
	wl = np.ones(len(cdl))
	print(str(len(cdl))+' randdom objects going out for paircounts')
	fdo = open(dirpcadw+ 'r'+file+zw+'pcadw.dat','w')
	for i in range(0,len(cdl)):
		fdo.write(str(sral[i])+' '+str(cral[i])+' '+str(sdecl[i])+' '+str(cdecl[i])+' '+str(cdl[i])+' '+str(wl[i])+'\n')

	fdo.close()

	
	return True

def ppxilcalc_LSDfjack_bs(sample,tile,date,zmin=.5,zmax=1.1,bs=5,start=0,rmaxf=250,rmax=50,mumin=0,mumax=1.,wmu='counts',mom=0):
	fl = sample+tile+'_'+date+'_zm'+str(zmin)+'zx'+str(zmax)
	DDnl = []	
	DDnorml = 0
	DDnormt = 0
	DRnl = []
	DRnorml = 0
	DRnormt = 0
	RRnl = []
	RRnl0 = []
	nmubin = 100
	nbin = rmax/bs
	for i in range(0,rmaxf*nmubin):
		DDnl.append(0)
		DRnl.append(0)
		RRnl.append(0)
		RRnl0.append(0)
	RRnorml = 0
	RRnormt = 0
	pl = []
	nmut = 0
	for i in range(0,nmubin):
		mu = i/float(nmubin)+.5/float(nmubin)
		mub = int(mu*nmubin)
		##print mu
		if mu > mumin and mu < mumax:
			pl.append((1.,P2(mu),P4(mu),P6(mu),P8(mu),mub))
			nmut += 1.
		else:
			pl.append((0,0,0,0,0,0))	
	fdp = open(dirpc+'g'+fl+'g'+fl+'2ptdmu.dat').readlines()
	DDnormt += float(fdp[0])
	fdnp = open(dirpc+'g'+fl+'r'+fl+'2ptdmu.dat').readlines()
	fr = open(dirpc+'r'+fl+'r'+fl+'2ptdmu.dat').readlines()
	DRnormt += float(fdnp[0])
	RRnormt += float(fr[0])
	for k in range(1,len(fdp)):
		dp = float(fdp[k])
		dr = float(fdnp[k])
		rp = float(fr[k])
		DDnl[k-1] += dp
		DRnl[k-1] += dr
		RRnl[k-1] += rp
					
	xil = np.zeros(int(nbin),'f')
	for i in range(start,rmax,bs):
		xi = 0
		dd = 0
		dr = 0
		rr = 0
	
		ddt = 0
		drt = 0
		rrt = 0
		for j in range(0,nmubin):
			if wmu != 'counts':
				dd = 0
				dr = 0
				rr = 0
			for k in range(0,bs):
				bin = nmubin*(i+k)+j			
				if bin < len(RRnl):
					#if RRnl[bin] == 0:
					#	pass
			
					#else:
					dd += DDnl[bin]
					rr += RRnl[bin]
					dr += DRnl[bin]
					ddt +=dd
					rrt += rr
					drt += dr
			
			#if rr != 0 and wm == 'muw':			
			if wmu != 'counts':
				xi += pl[j][mom]/float(nmut)*(dd/DDnormt-2*dr/DRnormt+rr/RRnormt)*RRnormt/rr
		if wmu == 'counts':
			xi = (dd/DDnormt-2*dr/DRnormt+rr/RRnormt)*RRnormt/rr		
		if i/bs < nbin:
			xil[i//bs] = xi
		#print ddt/DDnormt,drt/DRnormt,rrt/RRnormt
	rl = []
	for i in range(0,len(xil)):
		rl.append(start+bs/2.+bs*i)
	rl = np.array(rl)
	#plt.plot(rl,xil)
	#plt.show()
	bsst = str(bs)+'st'+str(start)
	fo = open('xi'+fl+bsst+'.dat','w')
	for i in range(0,len(rl)):
		fo.write(str(rl[i])+' '+str(xil[i])+'\n')
	fo.close()		
	return xil

def plotxi():
	dl = np.loadtxt('/Users/ashleyross/Dropbox/DESI/xiLRG70003_20200219_zm0.5zx1.1bsc.dat').transpose() 
	de = np.loadtxt('/Users/ashleyross/Dropbox/DESI/xiELG70004_20200219_zm0.8zx1.6bsc.dat').transpose()
	ml = (dl[0]/7.78)**-1.98
	plt.loglog(dl[0],dl[1],'r-',label='MINI_SV_LRG, Tile 70003')
	plt.loglog(de[0],de[1],'b-',label='MINI_SV_ELG, Tile 70004')
	plt.loglog(dl[0],ml,'r:',label=r'$(r/7.78)^{-1.98}$ (Kitanidis et al.)')
	plt.legend()
	plt.xlabel(r'$r$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$\xi$')
	plt.title('From 202200219')
	plt.savefig('/Users/ashleyross/Dropbox/DESI/miniSVxi.png')
	plt.show()


if __name__ == '__main__':
	createSourcesrd_ad('LRG','70003','20200219')
	createSourcesrd_ad('ELG','70004','20200219',zmin=.8,zmax=1.6)