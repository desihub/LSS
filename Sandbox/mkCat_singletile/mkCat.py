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

ffd = dirout+type+str(tile)+'_'+night+'_full.dat.fits'
ffr = dirout+type+str(tile)+'_'+night+'_full.ran.fits'
fcd = dirout+type+str(tile)+'_'+night+'_clustering.dat.fits'
fcr = dirout+type+str(tile)+'_'+night+'_clustering.ran.fits'


if type != 'ELG':
	mtlf = tardir+'MTL_Tile_'+str(tile)+'_0.37.0_all.fits'
else:
	mtlf = tardir+'MTL_TILE_ELG_'+str(tile)+'_0.37.0.fits'	

print('using '+mtlf +' as the mtl file; IS THAT CORRECT?')

elgandlrgbits = [1,5,6,7,8,9,11,12,13] #these get used to veto imaging area

zfailmd = 'zwarn' #only option so far, but can easily add things based on delta_chi2 or whatever
weightmd = 'wloc' #only option so far, weight observed redshifts by number of targets that wanted fiber


mkfulld = True
mkfullr = False
mkclus = True

'''
Will need to add in lines for running fiber assign on randoms for future observations
Not eager to add them in now since old data was observed with bugged version of fiberassign
'''

if mkfulld:
    tspec = ct.combspecdata(tile,night,coaddir)
    pdict,goodloc = ct.goodlocdict(tspec)
    tfa = ct.gettarinfo_type(fadir,tile,goodloc,mtlf,tarbit,tp=tp)
    
    tout = join(tfa,tspec,keys=['TARGETID','LOCATION','PRIORITY'],join_type='left') #targetid should be enough, but all three are in both and should be the same
    wz = tout['ZWARN']*0 == 0
    print('there are '+str(len(tout[wz]))+' rows with spec obs redshifts')
    print(tout.dtype.names)

    
    tout.write(ffd,format='fits', overwrite=True) 
    print('wrote matched targets/redshifts to '+ffd)
    
if mkfullr:
    tspec = ct.combspecdata(tile,night,coaddir)
    pdict,goodloc = ct.goodlocdict(tspec)
    ranall = ct.mkfullran(tile,goodloc,pdict,randir)
    #fout = dirout+type+str(tile)+'_'+night+'_full.ran.fits'
    ranall.write(ffr,format='fits', overwrite=True)

if mkclus:
    maxp,loc_fail = ct.mkclusdat(ffd,fcd,zfailmd,weightmd,maskbits=elgandlrgbits)    
       
    ct.mkclusran(ffr,fcr,fcd,maxp,loc_fail,maskbits=elgandlrgbits)
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

