import fitsio
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import os

matplotlib.rcParams['font.size'] = 12

njack = '60'
trs = ['ELG_LOPnotqso','QSO','LRG','BGS_BRIGHT']
bsl = [10,10,5,5]
dirxi = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/2.1/xi/smu/'
xit = 'poles'
baotempdir = '/global/homes/a/ajross/code/BAOfit_xs/BAOtemplates/'

outdir = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/2.1/xi/plots/blinded/'
xlim = 15,185
xlab = 'scale (s; arbitrary units)'
dirc = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/2.1/'
cfn = fitsio.read(dirc+'LRG_N_clustering.dat.fits')
cfs = fitsio.read(dirc+'LRG_S_clustering.dat.fits')
prelim = 'Raw Data \n Not For Scientific Analysis'
labeld = 'data with jack-knife errors'

tr = 'LRG'
zw = '0.4_1.1'
bs = '5'
pt = 's'
cl = 'firebrick'
mfc = cl
modpar = 'DESI0.43.04.08.015.00'
ylim0 = -50,110
bv = 2

def plotxis():
	fn_txt = dirxi+'xi'+xit+'_'+tr+'_NScomb_'+zw+'_default_FKP_lin'+str(bs)+'_njack'+njack+'.txt'
	xi = np.loadtxt(fn_txt).transpose()
	plt.errorbar(xi[0],xi[0]**2.*xi[2],xi[0]**2.*xi[5],fmt=pt,color=cl,label=labeld,mfc=mfc)
	ax = plt.gca()
	ax.axes.xaxis.set_ticks([])
	ax.axes.yaxis.set_ticks([])
	plt.xlabel(xlab)
	plt.ylabel(r'$s^2\times$ configuration-space monopole')
	(zmin,zmax) = zw.strip('lowz').split('_')
	cfn = fitsio.read(dirc+tr+'_N_clustering.dat.fits')
	cfs = fitsio.read(dirc+tr+'_S_clustering.dat.fits')

	seln = cfn['Z'] > float(zmin)
	seln &= cfn['Z'] < float(zmax)
	sels = cfs['Z'] > float(zmin)
	sels &= cfs['Z'] < float(zmax)
	ngal = len(cfn[seln])+len(cfs[sels])
	plt.title('1st two months of DESI '+tr[:3]+'s; '+str(ngal)+ ' with '+zmin+'<z<'+zmax)
	xilin = np.loadtxt(baotempdir+'xi0'+modpar+'.dat').transpose()

	plt.plot(xilin[0],xilin[0]**2.*xilin[1]*bv,'k-.',label='a simplified BAO model (not a fit)')

	plt.xlim(xlim)
	plt.ylim(ylim0)
	plt.legend()
	plt.text(xlim[0]+5,ylim0[0]+20,prelim,rotation=35,size=24,alpha=0.2)
	plt.savefig(outdir+'xi0'+tr+zw+str(bs)+'_noaxisnum.png')
	plt.show()

	ax = plt.gca()
	ax.axes.xaxis.set_ticks([])
	ax.axes.yaxis.set_ticks([])

	plt.errorbar(xi[0],xi[0]*xi[3],xi[0]*xi[6],fmt=pt,color=cl,label=labeld,mfc=mfc)
	plt.xlabel(xlab)
	plt.ylabel(r'$s\times$ configuration-space quadrupole')
	#plt.minorticks_on()
	plt.title('1st two months of DESI '+tr[:3]+'s; '+str(ngal)+ ' with '+zmin+'<z<'+zmax)
	xilin = np.loadtxt(baotempdir+'xi2'+modpar+'.dat').transpose()
	plt.plot(xilin[0],xilin[0]*xilin[1]*bv,'k-.',label='a simplified BAO model (not a fit)')
	ylim = plt.gca().get_ylim()
	plt.text(xlim[0]+5,ylim[0]+0.1*(ylim[1]-ylim[0]),prelim,rotation=35,size=24,alpha=0.2)
	plt.legend()
	plt.xlim(xlim)
	#plt.ylim(-50,100)
	plt.savefig(outdir+'xi2'+tr+zw+str(bs)+'_noaxisnum.png')
	plt.show()

plotxis()



tr = 'QSO'
zw = '0.8_2.1lowz'
bs = '10'
pt = '^'
cl = 'forestgreen'
mfc = cl
modpar = 'DESI0.3954715.00'
ylim0 = -50,75
bv = 1.2
plotxis()


tr = 'BGS_BRIGHT'
zw = '0.1_0.5'
bs = '5'
pt = 'o'
cl = 'k'
mfc = 'w'
modpar = 'DESI0.43.06.010.015.00'
ylim0 = -50,85
bv = 1.6
plotxis()

tr = 'ELG_LOPnotqso'
zw = '0.8_1.6'
bs = '10'
pt = '*'
cl = 'steelblue'
mfc = cl
modpar = 'Challenge_matterpower0.704815.00'
ylim0 = -25,50
bv = 0.6
baotempdir = os.environ['HOME']+'/code/BAOfit/BAOtemplates/'
plotxis()
