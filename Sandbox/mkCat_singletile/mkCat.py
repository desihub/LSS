'''
one executable to catalogs from data from a single tile
'''



#standard python
import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import glob
import argparse
from astropy.table import Table,join,unique,vstack

#from this package
import cattools as ct


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--tile", help="observed tile to use")
parser.add_argument("--night", help="date of observation")
args = parser.parse_args()

type = args.type
tile = args.tile
night = args.night

if type == 'LRG':
	tarbit = 10 #targeting bit
	pr = 4000 #priority; anything with higher priority vetos fiber in randoms
if type == 'QSO':
	tarbit = 12
	pr = 10000
if type == 'ELG':
	tarbit = 11
	pr = 10000


print(type,tile,night)
tp = 'CMX_TARGET'
print('targeting bit, priority, target type; CHECK THEY ARE CORRECT!')
print(tarbit,pr,tp)

sys.path.append("../")

minisvdir = '/project/projectdirs/desi/users/ajross/catalogs/minisv2/'
dirout = minisvdir+'LSScats/test/'
randir = minisvdir+'random/'
tardir = minisvdir+'targets/'
fadir = tardir
coaddir = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/'

if type != 'ELG':
	mtlf = tardir+'MTL_Tile_'+str(tile)+'_0.37.0_all.fits'
else:
	mtlf = tardir+'MTL_TILE_ELG_'+str(tile)+'_0.37.0.fits'	

print('using '+mtlf +' as the mtl file; IS THAT CORRECT?')

elgandlrgbits = [1,5,6,7,8,9,11,12,13] #these get used to veto imaging area



mkfulld = False
mkfullr = True

'''
Will need to add in lines for running fiber assign on randoms for future observations
Not eager to add them in now since old data was observed with bugged version of fiberassign
'''

if mkfulld:
    tspec = ct.combspecdata(tile,night,coaddir)
    pdict,goodloc = ct.goodlocdict(tspec)
    tfa = ct.gettarinfo_type(fadir,tile,goodloc,mtlf,tarbit,tp=tp)
    
    tout = join(tfa,tspec,keys=['TARGETID'],join_type='left')
    wz = tout['ZWARN']*0 == 0
    print('there are '+str(len(tout[wz]))+' rows with spec obs redshifts')

    fout = dirout+type+str(tile)+'_'+night+'_full.dat.fits'
    tout.write(fout,format='fits', overwrite=True) 
    print('wrote matched targets/redshifts to '+fout)
    
if mkfullr:
    tspec = ct.combspecdata(tile,night,coaddir)
    pdict,goodloc = ct.goodlocdict(tspec)
    ranall = ct.mkfullran(tile,goodloc,pdict,randir)
    fout = dirout+type+str(tile)+'_'+night+'_full.ran.fits'
    ranall.write(fout,format='fits', overwrite=True)
    
       

# dfout = dirout+type +str(tile)+'_'+night+'_clustering.dat.fits'
# rf = dirout+type +str(tile)+'_'+night+'_full.ran.fits'
# rfout = dirout+type +str(tile)+'_'+night+'_clustering.ran.fits'	
# dd = fitsio.read(df)	
# print(np.unique(dd['ZWARN']))
# maxp = np.max(dd['PRIORITY'])
# 
# dchi2 = 0
# 
# if type != 'LRG':
# 	wfail = (dd['ZWARN'] != 999999) & (dd['ZWARN'] > 0)	
# else:
# 	wfail = (dd['ZWARN'] != 999999) & ((dd['DELTACHI2'] <= dchi2) | (dd['ZWARN'] > 0)	)
# 
# loc_fail = dd[wfail]['LOCATION']	
# print(len(loc_fail))
# 
# ddm = cutphotmask(dd)
# nl = countloc(ddm)
# 
# wg = (ddm['ZWARN'] == 0) 
# if type == 'LRG':
# 	wg &= ddm['DELTACHI2'] > dchi2
# 
# ddzg = ddm[wg]
# 
# print('clustering catalog will have '+str(len(ddzg))+ ' objects in it')
# 
# ddclus = Table()
# ddclus['RA'] = ddzg['RA']
# ddclus['DEC'] = ddzg['DEC']
# ddclus['Z'] = ddzg['Z']
# ddclus['WEIGHT'] = assignweights(ddzg,nl)
# 
# print('minimum,maximum weight')
# print(np.min(ddclus['WEIGHT']),np.max(ddclus['WEIGHT']))	
# 
# ddclus.write(dfout,format='fits',overwrite=True)
# 
# #plt.hist(ddclus['Z'],normed=True,bins=20,range=(0.5,1.1),histtype='step')
# #plt.hist(ddclus['Z'],weights=ddclus['WEIGHT'],normed=True,bins=20,range=(0.5,1.1),histtype='step')
# #plt.show()
# 
# dr = fitsio.read(rf)
# drm = cutphotmask(dr)
# 
# wpr = drm['PRIORITY'] <= maxp
# wzf = np.isin(drm['LOCATION'],loc_fail)
# wzt = wpr & ~wzf
# 
# drmz = drm[wzt]
# print(str(len(drmz))+' after cutting based on failures and priority')
# plt.plot(drmz['RA'],drmz['DEC'],'k,')
# plt.plot(drm[~wpr]['RA'],drm[~wpr]['DEC'],'b,')
# plt.plot(drm[wzf]['RA'],drm[wzf]['DEC'],'g,')
# plt.plot(ddclus['RA'],ddclus['DEC'],'r.')
# plt.show()
# rclus = Table()
# rclus['RA'] = drmz['RA']
# rclus['DEC'] = drmz['DEC']
# rclus['Z'] = drmz['Z']
# rclus['WEIGHT'] = np.ones(len(drmz))
# 
# rclus.write(rfout,format='fits',overwrite=True)

