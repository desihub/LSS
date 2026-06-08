import fitsio
import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table
import healpy as hp

def radec2thphi(ra,dec):
    return (-dec+90.)*np.pi/180.,ra*np.pi/180.

def get_stats(tp,veto='_noveto',prog='dark'):
    zcol = 'Z'
    if veto == '':
        zcol = 'Z_not4clus'
    print(zcol)    
    cols = ['RA','DEC',zcol,'NTILE','ZWARN','COMP_TILE','LOCATION_ASSIGNED','TILEID','GOODHARDLOC','DELTACHI2']
    if tp[:3] == 'ELG':
        cols.append('o2c')
    tilesall = Table.read("/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv")
    #pt = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/processed_tiles_'+prog+'.fits')
    dd = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'zdone_full'+veto+'.dat.fits',columns=cols)
    tids = np.unique(dd['TILEID'])
    asd = dd['LOCATION_ASSIGNED'] == 1
    asd &= dd['ZWARN'] != 999999
    asd &= dd['GOODHARDLOC'] == 1
    dd = dd[asd]
    wz = dd['ZWARN'] == 0
    if tp[:3] == 'BGS':
        print('applying extra cut for BGS')
        wz &= dd['DELTACHI2'] > 40

    if tp[:3] == 'ELG':
        wz = dd['o2c'] > 0.9
        wz &= dd['ZWARN']*0 == 0
        wz &= dd['ZWARN'] != 999999

    if tp == 'LRG':
        wz = dd['ZWARN'] == 0
        wz &= dd['ZWARN']*0 == 0
        wz &= dd['ZWARN'] != 999999
        wz &= dd[zcol]<1.5
        wz &= dd['DELTACHI2']>15
        zmin = 0.4
        zmax = 1.1
    if tp == 'QSO':
        #good redshifts are currently just the ones that should have been defined in the QSO file when merged in full
        wz = dd['Z']*0 == 0
        wz &= dd['Z'] != 999999
        wz &= dd['Z'] != 1.e20
        wz &= dd['ZWARN'] != 999999
        zmin = 0.8
        zmax = 3.5

    print(len(dd),len(dd[wz]))
    d = Table.read('/global/cfs/cdirs/desi/users/raichoor/main-status/skymaps/main-skymap-'+prog+'-goal.fits')
    tilesprog = tilesall[np.in1d(tilesall["TILEID"], tids)]               
    ntileobs = np.zeros(len(d), dtype=int) # store the nb of obs. tiles per pixel
    for p in range(d["TILEIDS"].shape[1]):
        sel = np.in1d(d["TILEIDS"][:, p], tilesprog["TILEID"][tilesprog["PASS"] == p])
        ntileobs[sel] += 1
    donefrac = np.zeros(1024*1024*12)
    sel = d["NPASS"] > 0
    donefrac[sel] = ntileobs[sel]/d[sel]["NPASS"]
    th,phi = radec2thphi(dd['RA'],dd['DEC'])
    pix = hp.ang2pix(1024,th,phi,nest=True)
    apx = np.arange(12*1024*1024)
    nol = np.unique(ntileobs)
    dl = []
    nt = 0
    ntz = 0
    for no in nol:
        if no > 0:
            spx = ntileobs == no
            px = apx[spx]
            seld = np.in1d(pix,px)
            print('number of '+tp+' in regions observed '+str(no)+' times is:'+str(len(dd[seld])))
            print('number of '+tp+' w. good z in regions observed '+str(no)+' times is:'+str(len(dd[seld&wz])))
            nppx = len(dd[seld])/len(apx[spx])
            nppxz = len(dd[seld&wz])/len(apx[spx])
            
            dl.append(nppx)
            sn = d['NPASS'] == no
            nf = len(d[sn])*nppx
            nfz = len(d[sn])*nppxz
            print('final expected number of '+tp+' in regions observed '+str(no)+' times is:'+str(nf))
            print('final expected number of '+tp+' w. good z in regions observed '+str(no)+' times is:'+str(nfz))
            nt += nf
            ntz += nfz
    print('total expected number is '+str(nt))
    print('total expected number w. good z is '+str(ntz))
              
tpl = ['LRG','QSO','ELG','ELG_LOPnotqso','BGS_ANY']
prog = 'dark'
for tp in tpl:
    if tp[:3] == 'BGS':
        prog = 'bright'
    get_stats(tp,prog=prog)                                
    