import fitsio
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import glob
import argparse

from astropy.table import Table,vstack,join,unique
from desitarget.targetmask import zwarn_mask


parser = argparse.ArgumentParser()
parser.add_argument("--prog", help="look for bright or dark tiles",default='dark')
parser.add_argument("--newdir", help="directory with new reductions, must include final /",default='/global/cfs/cdirs/desi/users/sjbailey/spectro/redux/skycte2/tiles/pernight/')
parser.add_argument("--fiddir", help="directory with fiducial reductions",default='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/')
args = parser.parse_args()

#get tile list
fls = glob.glob(args.newdir+'*')
tls = []
for fl in fls:
    tl = fl.strip(args.newdir)
    tl = int(tl)
    tls.append(tl)

print('found '+str(len(tls))+' tiles in new reductions')

mt = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
wd = mt['SURVEY'] == 'main'
wd &= mt['OBSSTATUS'] == 'obsend'
wd &= mt['FAPRGRM'] == args.prog
wd &= np.isin(mt['TILEID'],tls)
mtd = mt[wd]
print('found '+str(len(mtd))+' '+args.prog+' tiles with obsend status in new reductions')

tiles4comb = Table()
tiles4comb['TILEID'] = mtd['TILEID']
tiles4comb['ZDATE'] = mtd['ARCHIVEDATE']
tiles4comb['THRUDATE'] = mtd['LASTNIGHT']


def combspecdata_simp(tile,tdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/',thru='thru' ):
    #put data from different spectrographs together, one table for fibermap, other for z
    tdate = str(tdate)
    specs = []
    #find out which spectrograph have data
    zfn = 'zbest'
    zhdu = 'ZBEST'
    shdu = 'SCORES'
    if int(tdate) >  20210730:
        zfn = 'redrock'
        zhdu = 'REDSHIFTS'
        #shdu = 'TSNR2' 
        

    for si in range(0,10):
        ff = coaddir+str(tile)+'/'+tdate+'/'+zfn+'-'+str(si)+'-'+str(tile)+'-'+thru+tdate+'.fits'
        if os.path.isfile(ff):
            fq = coaddir+str(tile)+'/'+tdate+'/zmtl-'+str(si)+'-'+str(tile)+'-'+thru+tdate+'.fits'
            if os.path.isfile(fq):

                specs.append(si)
            else:
                print('did not find '+fq)    
        elif zfn == 'zbest':
            zfnt = 'redrock'
            ff = coaddir+str(tile)+'/'+tdate+'/'+zfnt+'-'+str(si)+'-'+str(tile)+'-'+thru+tdate+'.fits'
            if os.path.isfile(ff):
                fq = coaddir+str(tile)+'/'+tdate+'/zmtl-'+str(si)+'-'+str(tile)+'-'+thru+tdate+'.fits'
                zfn = zfnt
                zhdu = 'REDSHIFTS'
                if os.path.isfile(fq):

                    specs.append(si)
                else:
                    print('did not find '+fq)    
            else:
                print('did not find '+ff)            
        else:
            print('did not find '+ff)        
    print('spectrographs with data on tile '+str(tile)+':')
    print(specs)            
    if len(specs) == 0:
        return None
    for i in range(0,len(specs)):
        tn = Table.read(coaddir+str(tile)+'/'+tdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-'+trhu+tdate+'.fits',hdu=zhdu)
        tnq = Table.read(coaddir+str(tile)+'/'+tdate+'/zmtl-'+str(specs[i])+'-'+str(tile)+'-'+thru+tdate+'.fits')
        tnf = Table.read(coaddir+str(tile)+'/'+tdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-'+thru+tdate+'.fits',hdu='FIBERMAP')
        tns = Table.read(coaddir+str(tile)+'/'+tdate+'/coadd-'+str(specs[i])+'-'+str(tile)+'-'+thru+tdate+'.fits',hdu=shdu)
    
        if i == 0:
           tspec = tn
           tq = tnq
           tf = tnf
           ts = tns
        else:    
            ts = vstack([ts,tns],metadata_conflicts='silent')
            tq = vstack([tq,tnq],metadata_conflicts='silent')
            tspec = vstack([tspec,tn],metadata_conflicts='silent')
            tf = vstack([tf,tnf],metadata_conflicts='silent')
        
    
    tf = unique(tf,keys=['TARGETID'])
    tq.keep_columns(['TARGETID','Z_QN','Z_QN_CONF','IS_QSO_QN','ZWARN'])
    tq['ZWARN'].name = 'ZWARN_MTL'
    tspec = join(tspec,tf,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')
    tspec = join(tspec,ts,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')
    tspec = join(tspec,tq,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')

    print(len(tspec),len(tf))
    #tspec['LOCATION'] = tf['LOCATION']
    #tspec['FIBERSTATUS'] = tf['FIBERSTATUS']
    #tspec['PRIORITY'] = tf['PRIORITY']
    return tspec

def combtiles(tiles,coadddir,thru):
    s = 0
    n = 0
    nfail = 0

    for tile,tdate in zip(tiles['TILEID'],tiles['THRUDATE']):
        tdate = str(tdate)
        tspec = combspecdata_simp(tile,tdate,coadddir,thru=thru)
        if tspec:
            tspec['TILEID'] = tile
            if s == 0:
                specd = tspec
                s = 1
            else:
                specd = vstack([specd,tspec],metadata_conflicts='silent')
            specd.sort('TARGETID')
            kp = (specd['TARGETID'] > 0)
            specd = specd[kp]

            n += 1
            print(tile,n,len(tiles),len(specd)) 
        else:
                print(str(tile)+' failed')
                nfail += 1  
    return specd
#get fid and new data

specdnew = combtiles(tiles4comb,args.newdir,thru='')
specdfid = combtiles(tiles4comb,args.fiddir)

print('comparing lengths of combined data; old,new:')
print(len(specdnew),len(specdold))
