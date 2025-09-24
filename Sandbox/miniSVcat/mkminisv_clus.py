'''
script to take full files and make clustering for type, tile, night to be edited below
'''

from astropy.table import Table, join,unique,vstack
import numpy as np
import fitsio
from matplotlib import pyplot as plt
import sys

minisvdir = '/project/projectdirs/desi/users/ajross/catalogs/minisv2/'
dirout = minisvdir+'LSScats/'

try:
	type = str(sys.argv[1])
	print(type)
	tile = int(sys.argv[2])
	print(tile)
	night = str(sys.argv[3])
	print(night)
except:
	print('requires three arguments: type=LRG/QSO/ELG, tile, night')	

elgandlrgbits = [1,5,6,7,8,9,11,12,13]


def cutphotmask(aa):
	keep = (aa['NOBS_G']>0) & (aa['NOBS_R']>0) & (aa['NOBS_Z']>0)
	for biti in elgandlrgbits:
		keep &= ((aa['MASKBITS'] & 2**biti)==0)
	aa = aa[keep]
	print(str(len(aa)) +' after imaging veto' )
	return aa
	
def countloc(aa):
	locs = aa['LOCATION']
	la = np.max(locs)+1
	nl = np.zeros(la)
	for i in range(0,len(aa)):
		nl[locs[i]] += 1
	return nl

def assignweights(aa,nl):
	wts = np.ones(len(aa))
	for i in range(0,len(aa)):
		loc = aa[i]['LOCATION']
		wts[i] = nl[loc]
	return wts	

df = dirout+type +str(tile)+'_'+night+'_full.dat.fits'
dfout = dirout+type +str(tile)+'_'+night+'_clustering.dat.fits'
rf = dirout+type +str(tile)+'_'+night+'_full.ran.fits'
rfout = dirout+type +str(tile)+'_'+night+'_clustering.ran.fits'	
dd = fitsio.read(df)	
print(np.unique(dd['ZWARN']))
maxp = np.max(dd['PRIORITY'])

dchi2 = 0

if type != 'LRG':
	wfail = (dd['ZWARN'] != 999999) & (dd['ZWARN'] > 0)	
else:
	wfail = (dd['ZWARN'] != 999999) & ((dd['DELTACHI2'] <= dchi2) | (dd['ZWARN'] > 0)	)

loc_fail = dd[wfail]['LOCATION']	
print(len(loc_fail))

ddm = cutphotmask(dd)
nl = countloc(ddm)

wg = (ddm['ZWARN'] == 0) 
if type == 'LRG':
	wg &= ddm['DELTACHI2'] > dchi2

ddzg = ddm[wg]

print('clustering catalog will have '+str(len(ddzg))+ ' objects in it')

ddclus = Table()
ddclus['RA'] = ddzg['RA']
ddclus['DEC'] = ddzg['DEC']
ddclus['Z'] = ddzg['Z']
ddclus['WEIGHT'] = assignweights(ddzg,nl)

print('minimum,maximum weight')
print(np.min(ddclus['WEIGHT']),np.max(ddclus['WEIGHT']))	

ddclus.write(dfout,format='fits',overwrite=True)

#plt.hist(ddclus['Z'],normed=True,bins=20,range=(0.5,1.1),histtype='step')
#plt.hist(ddclus['Z'],weights=ddclus['WEIGHT'],normed=True,bins=20,range=(0.5,1.1),histtype='step')
#plt.show()

dr = fitsio.read(rf)
drm = cutphotmask(dr)

wpr = drm['PRIORITY'] <= maxp
wzf = np.isin(drm['LOCATION'],loc_fail)
wzt = wpr & ~wzf

drmz = drm[wzt]
print(str(len(drmz))+' after cutting based on failures and priority')
plt.plot(drmz['RA'],drmz['DEC'],'k,')
plt.plot(drm[~wpr]['RA'],drm[~wpr]['DEC'],'b,')
plt.plot(drm[wzf]['RA'],drm[wzf]['DEC'],'g,')
plt.plot(ddclus['RA'],ddclus['DEC'],'r.')
plt.show()
rclus = Table()
rclus['RA'] = drmz['RA']
rclus['DEC'] = drmz['DEC']
rclus['Z'] = drmz['Z']
rclus['WEIGHT'] = np.ones(len(drmz))

rclus.write(rfout,format='fits',overwrite=True)

