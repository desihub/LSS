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
from matplotlib import pyplot as plt

#from this package
import cattools as ct
import fa4lsscat as fa
import xitools as xt


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--tile", help="observed tile to use")
parser.add_argument("--night", help="date of observation")
parser.add_argument("--fadate", help="date for fiberassign run")
args = parser.parse_args()

type = args.type
tile = args.tile
night = args.night
fadate = args.fadate

if type == 'LRG':
    tarbit = 0 #targeting bit
    pr = 4000 #priority; anything with higher priority vetos fiber in randoms
if type == 'QSO':
    tarbit = 2
    pr = 10000
if type == 'ELG':
    tarbit = 1
    pr = 10000


print(type,tile,night)
tp = 'SV1_DESI_TARGET'
print('targeting bit, priority, target type; CHECK THEY ARE CORRECT!')
print(tarbit,pr,tp)

sys.path.append("../")
#print(sys.path)

svdir = '/project/projectdirs/desi/users/ajross/catalogs/SV/'
dirout = svdir+'LSScats/test/'
randir = svdir+'random'
rm = 0
rx = 10
for i in range(rm,rx):
    if not os.path.exists(svdir+'random'+str(i)):
        os.mkdir(svdir+'random'+str(i))
        print('made '+str(i)+' random directory')
fadir = '/global/cfs/cdirs/desi/users/raichoor/fiberassign-sv1/'+fadate+'/'
tardir = fadir
coaddir = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/'

ffd = dirout+type+str(tile)+'_'+night+'_full.dat.fits'

fcd = dirout+type+str(tile)+'_'+night+'_clustering.dat.fits'

mtlf = fadir+'/0'+tile+'-targ.fits'

print('using '+mtlf +' as the mtl file; IS THAT CORRECT?')

elgandlrgbits = [1,5,6,7,8,9,11,12,13] #these get used to veto imaging area

zfailmd = 'zwarn' #only option so far, but can easily add things based on delta_chi2 or whatever
weightmd = 'wloc' #only option so far, weight observed redshifts by number of targets that wanted fiber

mkranmtl = False #make a mtl file of randoms
runrfa = False #run randoms through fiberassign
mkfulld = True
mkfullr = True
mkclus = True
docatplots = True
doclus = True

tilef = fadir+'0'+tile+'-tiles.fits' #the tile file
fbaf = fadir+'fba-0'+tile+'.fits' #the tile file

if mkranmtl: #this cuts the random file to the tile and adds columns necessary for fiberassign, done here it is very inefficient (would be better to do all tiles at once)
    for i in range(rm,rx):
        ct.randomtilesi(tilef ,randir,i)

if runrfa:
    fbah = fitsio.read_header(fbaf)
    dt = fbah['FA_RUN']
    for i in range(rm,rx):
        fa.getfatiles(randir+str(i)+'/tilenofa-'+str(tile)+'.fits',tilef,dirout=randir+str(i)+'/',dt = dt)

if mkfulld:
    tspec = ct.combspecdata(tile,night,coaddir)
    pdict,goodloc = ct.goodlocdict(tspec)
    tfa = ct.gettarinfo_type(fadir,tile,goodloc,mtlf,tarbit,tp=tp)
    print(tspec.dtype.names)
    tout = join(tfa,tspec,keys=['TARGETID','LOCATION','PRIORITY'],join_type='left') #targetid should be enough, but all three are in both and should be the same
    print(tout.dtype.names)
    wz = tout['ZWARN']*0 == 0
    wzg = tout['ZWARN'] == 0
    print('there are '+str(len(tout[wz]))+' rows with spec obs redshifts and '+str(len(tout[wzg]))+' with zwarn=0')
        
    tout.write(ffd,format='fits', overwrite=True) 
    print('wrote matched targets/redshifts to '+ffd)
    
if mkfullr:
    tspec = ct.combspecdata(tile,night,coaddir)
    pdict,goodloc = ct.goodlocdict(tspec)
    for i in range(rm,rx):
        ranall = ct.mkfullran(tile,goodloc,pdict,randir+str(i)+'/')
        #fout = dirout+type+str(tile)+'_'+night+'_full.ran.fits'
        ffr = dirout+type+str(tile)+'_'+night+'_'+str(i)+'_full.ran.fits'
        ranall.write(ffr,format='fits', overwrite=True)

if mkclus:
    maxp,loc_fail = ct.mkclusdat(ffd,fcd,zfailmd,weightmd,maskbits=elgandlrgbits)    
    for i in range(rm,rx):
        ffr = dirout+type+str(tile)+'_'+night+'_'+str(i)+'_full.ran.fits'
        fcr = dirout+type+str(tile)+'_'+night+'_'+str(i)+'_clustering.ran.fits'      
        ct.mkclusran(ffr,fcr,fcd,maxp,loc_fail,maskbits=elgandlrgbits)

if docatplots:
    ii = 0
    fd = fitsio.read(ffd)
    plt.plot(fd['RA'],fd['DEC'],'r.',label='potential targets')
    fc = fitsio.read(fcd)
    plt.plot(fc['RA'],fc['DEC'],'bo',label='good redshifts')
    ffr = fitsio.read(dirout+type+str(tile)+'_'+night+'_'+str(ii)+'_clustering.ran.fits')
    plt.plot(ffr['RA'],ffr['DEC'],'k,',label='randoms')

    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.title(type+' on tile '+tile+' observed '+night)
    plt.legend()
    plt.show()
    if type == 'ELG':
        zr = (.3,2)
    if type == 'LRG':
        zr = (.4,1.7)
    if type == 'QSO':
        zr = (.1,4.5)    
    plt.hist(fc['Z'],bins=100,range=zr,histtype='step')
    plt.xlabel('redshift')
    plt.ylabel('# with zwarn == 0')
    plt.title(type+' on tile '+tile+' observed '+night)
    plt.show()

if doclus:
	import subprocess
	if type == 'ELG':
		zmin = .8
		zmax = 1.6
	if type == 'LRG':
		zmin = .5
		zmax = 1.1
	if type == 'QSO':
		zmin = 1.
		zmax = 2.

	rmax = 10
	gf = xt.createSourcesrd_ad(type,tile,night,zmin=zmin,zmax=zmax,datadir=dirout)
	subprocess.run(['chmod','+x','dopc'+gf+'.sh'])
	subprocess.run('./dopc'+gf+'.sh')
	for i in range(rm+1,rmax):
		gf = xt.createSourcesrd_ari(type,tile,night,i,zmin=zmin,zmax=zmax,datadir=dirout)
		subprocess.run(['chmod','+x','dopc'+gf+'.sh'])
		subprocess.run('./dopc'+gf+'.sh')
	xt.ppxilcalc_LSDfjack_bs(type,tile,night,zmin=zmin,zmax=zmax,nran=rmax)
	xt.ppxilcalc_LSDfjack_bs(type,tile,night,zmin=zmin,zmax=zmax,bs=5,nran=rmax)
        
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

