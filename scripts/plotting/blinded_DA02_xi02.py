import fitsio
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import os

matplotlib.rcParams['font.size'] = 12

njack = '60'
trs = ['ELG_LOPnotqso','QSO','LRG','BGS_BRIGHT']
bsl = [10,10,5,5]
dirxi = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/test/xi/smu/'
xit = 'poles'
outdir = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/test/xi/plots/blinded/'
xlim = 15,185
xlab = 'scale (s; arbitrary units)'
dirc = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/test/'
cfn = fitsio.read(dirc+'LRG_N_clustering.dat.fits')
cfs = fitsio.read(dirc+'LRG_S_clustering.dat.fits')
tr = 'LRG'
zw = '0.4_1.1'
bs = '5'
prelim = 'Raw Data. \n Not For Scientific Analysis.'

fn_txt = dirxi+'xi'+xit+'_'+tr+'_NScomb_'+zw+'_default_FKP_lin'+str(bs)+'_njack'+njack+'.txt'
xi = np.loadtxt(fn_txt).transpose()
plt.errorbar(xi[0],xi[0]**2.*xi[2],xi[0]**2.*xi[5],fmt='ro',label='data with jack-knife errors')
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
plt.xlabel(xlab)
plt.ylabel(r'$s^2\times$ configuration-space monopole')
(zmin,zmax) = zw.split('_')
seln = cfn['Z'] > float(zmin)
seln &= cfn['Z'] < float(zmax)
sels = cfs['Z'] > float(zmin)
sels &= cfs['Z'] < float(zmax)
ngal = len(cfn[seln])+len(cfs[sels])
plt.title('1st two months of DESI LRGs; '+str(ngal)+ ' with '+zmin+'<z<'+zmax)
xilin = np.loadtxt(os.environ['HOME']+'/code/BAOfit/BAOtemplates/xi0Challenge_matterpower0.43.04.08.015.00.dat').transpose()
plt.plot(xilin[0],xilin[0]**2.*xilin[1]*2,'k-.',label='a simplified BAO model (not a fit)')
ylim = -50,110
plt.xlim(xlim)
plt.ylim(ylim)
plt.legend()
plt.savefig(outdir+'xi0'+tr+zw+str(bs)+'_noaxisnum.png')
plt.text(xlim[0]+5,ylim[0]+20,prelim,rotation=35,size=24,alpha=0.2)
plt.show()
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

plt.errorbar(xi[0],xi[0]*xi[3],xi[0]*xi[6],fmt='ro',label='data with jack-knife errors')
plt.xlabel(xlab)
plt.ylabel(r'$s\times$ configuration-space quadrupole')
#plt.minorticks_on()
plt.title('1st two months of DESI LRGs; '+str(ngal)+ ' with '+zmin+'<z<'+zmax)
xilin = np.loadtxt(os.environ['HOME']+'/code/BAOfit/BAOtemplates/xi2Challenge_matterpower0.43.04.08.015.00.dat').transpose()
plt.plot(xilin[0],xilin[0]*xilin[1]*2,'k-.',label='a simplified BAO model (not a fit)')
ylim = plt.gca().get_ylim()
plt.text(xlim[0]+5,ylim[0]+0.1*(ylim[1]-ylim[0]),prelim,rotation=35,size=24,alpha=0.2)
plt.legend()
plt.xlim(xlim)
#plt.ylim(-50,100)
plt.savefig(outdir+'xi2'+tr+zw+str(bs)+'_noaxisnum.png')
plt.show()


