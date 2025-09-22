import fitsio
import datetime
import os
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--add", help="whether or not to add to the tracking file",default='n')
args = parser.parse_args()

today = datetime.date.today()

fn = 'temp.txt'
if args.add == 'y':
    fn = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/stats.txt'
   

if os.path.isfile(fn):
    fo = open(fn,'a')
else:
    fo = open(fn,'w')

fo.write('#####################\n')
fo.write('as of '+today.strftime("%B %d, %Y")+'\n')
fo.write('#####################\n')

tps = ['LGE','QSO','LRG','ELG','ELG_LOP','ELG_LOPnotqso','BGS_ANY','BGS_BRIGHT']

zcol = 'Z_not4clus'
for tp in tps:
    print(tp+':')
    prog='dark'
    if 'BGS' in tp:
        prog = 'bright'
    if 'LGE' in tp:
        prog = 'dark1b'
    rtnv = fitsio.read_header('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+prog+'_0_full_noveto.ran.fits',ext=1)
    areanv = rtnv['NAXIS2']/2500
    dtnv = fitsio.read_header('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'_full_noveto.dat.fits',ext=1)
    ngnv = dtnv['NAXIS2']
    print('number density pre veto '+str(ngnv/areanv))
    rt = fitsio.read_header('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'_0_full.ran.fits',ext=1)
    area = rt['NAXIS2']/2500
    cols = ['Z_not4clus','ZWARN','LOCATION_ASSIGNED','DELTACHI2']
    if tp[:3] == 'ELG':
        cols.append('o2c')
    dt = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'_full.dat.fits',columns=cols)
    print('number density after veto '+str(len(dt)/area))
    #wz = dt[zcol]*0 == 0
    #wz &= dt[zcol] != 999999
    wz = dt['ZWARN']*0 == 0
    wz &= dt['ZWARN'] != 1.e20
    wz &= dt['ZWARN'] != 999999
    wz &= dt['LOCATION_ASSIGNED'] == 1

    if tp == 'QSO':
        #good redshifts are currently just the ones that should have been defined in the QSO file when merged in full
        wg = dt[zcol]*0 == 0
        wg &= dt[zcol] != 999999
        wg &= dt[zcol] != 1.e20
    
    if tp[:3] == 'ELG':
        wg = dt['o2c'] > 0.9

    if tp == 'LGE':
        wg = dt['ZWARN'] == 0
        wg &= dt['DELTACHI2']>10
        print(np.sum(wz),np.sum(wg))
    
    if tp == 'LRG':
        # Custom DELTACHI2 vs z cut from Rongpu
        wg = dt['ZWARN'] == 0
        drz = (10**(3 - 3.5*dt[zcol]))
        mask_bad = (drz>30) & (dt['DELTACHI2']<30)
        mask_bad |= (drz<30) & (dt['DELTACHI2']<drz)
        mask_bad |= (dt['DELTACHI2']<10)
        wg &= dt[zcol]<1.4
        wg &= (~mask_bad)


    if tp[:3] == 'BGS':
        wg = dt['DELTACHI2'] > 40
    comp = len(dt[wz])/len(dt)#np.sum(1/dt['FRACZ_TILELOCID'][wz])
    print('area: '+str(area))
    print('# of good z: '+str(len(dt[wz&wg])))
    #print('completeness: '+str(round(len(dt[wz])/len(dt),3)))
    print('completeness: '+str(round(comp,3)))
    fo.write(tp+':\n')
    fo.write('area: '+str(area)+'\n')
    fo.write('# of good z: '+str(len(dt[wz&wg]))+'\n')
    #fo.write('completeness: '+str(round(len(dt[wz])/len(dt),3))+'\n')
    fo.write('completeness: '+str(round(comp,3))+'\n')
