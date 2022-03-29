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
parser.add_argument("--newdir", help="directory with new reductions, must include final /",default='/global/cfs/cdirs/desi/spectro/redux/f5/tiles/pernight/')
parser.add_argument("--fiddir", help="directory with fiducial reductions",default='/global/cfs/cdirs/desi/spectro/redux/everest/tiles/pernight/')
parser.add_argument("--output",help='full path for output file',default='temp.fits')
parser.add_argument("--mktable",help='whether or not to make output (starts from previous output if n)',default='y')
parser.add_argument("--survey",help='sv1, sv3, or main',default='main')
args = parser.parse_args()

if args.survey == 'sv1':
    tarcol = 'SV1_DESI_TARGET'
    bgscol = 'SV1_BGS_TARGET'
if args.survey == 'sv3':
    tarcol = 'SV3_DESI_TARGET'
    bgscol = 'SV3_BGS_TARGET'
if args.survey == 'main':
    tarcol = 'DESI_TARGET'
    bgscol = 'BGS_TARGET'        


if args.mktable == 'y':
    #get tile list
    fls = glob.glob(args.newdir+'*')
    tls = []
    for fl in fls:
        #tl = fl.strip(args.newdir)
        tl = fl.replace(args.newdir, '') 
        if len(tl) > 0:
            tl = int(tl)
            tls.append(tl)

    print('found '+str(len(tls))+' tiles in new reductions')

    mt = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
    wd = mt['SURVEY'] == args.survey
    #wd &= mt['OBSSTATUS'] == 'obsend'
    wd &= mt['FAPRGRM'] == args.prog
    wd &= np.isin(mt['TILEID'],tls)
    mtd = mt[wd]
    print('found '+str(len(mtd))+' '+args.prog+' tiles with obsend status in new reductions')

    tiles4comb = Table()
    tiles4comb['TILEID'] = mtd['TILEID']
    tiles4comb['ZDATE'] = mtd['ARCHIVEDATE']
    tiles4comb['THRUDATE'] = mtd['LASTNIGHT']


    def combspecdata_simp(tile,tdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/',thru='' ):
        #put data from different spectrographs together, one table for fibermap, other for z
        tdate = str(tdate)
        specs = []
        #find out which spectrograph have data
        shdu = 'SCORES'
        zfn = 'redrock'
        zhdu = 'REDSHIFTS'
            #shdu = 'TSNR2' 
        

        for si in range(0,10):
            ff = coaddir+str(tile)+'/'+tdate+'/'+zfn+'-'+str(si)+'-'+str(tile)+'-'+thru+tdate+'.fits'
            if os.path.isfile(ff):
                #fq = coaddir+str(tile)+'/'+tdate+'/zmtl-'+str(si)+'-'+str(tile)+'-'+thru+tdate+'.fits'
                #if os.path.isfile(fq):

                specs.append(si)
                #else:
                #    print('did not find '+fq)    
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
            tn = Table.read(coaddir+str(tile)+'/'+tdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-'+thru+tdate+'.fits',hdu=zhdu)
            #tnq = Table.read(coaddir+str(tile)+'/'+tdate+'/zmtl-'+str(specs[i])+'-'+str(tile)+'-'+thru+tdate+'.fits')
            tnf = Table.read(coaddir+str(tile)+'/'+tdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-'+thru+tdate+'.fits',hdu='FIBERMAP')
            tns = Table.read(coaddir+str(tile)+'/'+tdate+'/coadd-'+str(specs[i])+'-'+str(tile)+'-'+thru+tdate+'.fits',hdu=shdu)
    
            if i == 0:
               tspec = tn
               #tq = tnq
               tf = tnf
               ts = tns
            else:    
                ts = vstack([ts,tns],metadata_conflicts='silent')
                #tq = vstack([tq,tnq],metadata_conflicts='silent')
                tspec = vstack([tspec,tn],metadata_conflicts='silent')
                tf = vstack([tf,tnf],metadata_conflicts='silent')
        
    
        tf = unique(tf,keys=['TARGETID'])
        #tq.keep_columns(['TARGETID','Z_QN','Z_QN_CONF','IS_QSO_QN','ZWARN'])
        #tq['ZWARN'].name = 'ZWARN_MTL'
        tspec = join(tspec,tf,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')
        tspec = join(tspec,ts,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')
        #tspec = join(tspec,tq,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')

        print(len(tspec),len(tf))
        #tspec['LOCATION'] = tf['LOCATION']
        #tspec['FIBERSTATUS'] = tf['FIBERSTATUS']
        #tspec['PRIORITY'] = tf['PRIORITY']
        return tspec

    def combtiles(tiles,coadddir,thru=''):
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

    specdnew = combtiles(tiles4comb,args.newdir,'')
    gtls = np.isin(tiles4comb['TILEID'],specdnew['TILEID'])
    print('number of tiles where night matches thru night '+str(np.sum(gtls)))
    specdfid = combtiles(tiles4comb[gtls],args.fiddir)

    print('comparing lengths of combined data; old,new:')
    print(len(specdnew),len(specdfid))

    combt = join(specdfid,specdnew,keys=['TARGETID'],table_names=['fid', 'new'])

    print('length after join is '+str(len(combt)))

    combt.write(args.output,format='fits',overwrite=True)
    print('wrote joined table to '+args.output)
else:
    combt = Table.read(args.output)
    
sel = combt['COADD_NUMEXP_fid'] == combt['COADD_NUMEXP_new']
print('number of rows where COADD_NUMEXP matches:' +str(len(combt[sel])))

def checkQA(dat,ver):            
    nodata = dat["ZWARN_MTL_"+ver] & zwarn_mask["NODATA"] != 0
    num_nod = np.sum(nodata)
    print('number with no data '+str(num_nod))
    badqa = dat["ZWARN_MTL_"+ver] & zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA") != 0
    num_badqa = np.sum(badqa)
    print('number with bad qa '+str(num_badqa))
    nomtl = nodata | badqa
    wfqa = ~nomtl
    return wfqa

#sqafid = checkQA(combt,'fid')
#sqanew = checkQA(combt,'new')
sqafid = combt['COADD_FIBERSTATUS_fid'] == 0
sqanew = combt['COADD_FIBERSTATUS_new'] == 0
selt = sel & sqafid & sqanew

combpass = combt[selt]

print('number of rows after keeping only those with COADD_NUMEXP match and good QA for both '+str(len(combpass)))

if args.prog == 'bright':
    wbgs = combpass[bgscol+'_fid'] > 0
    bgstar = combpass[wbgs]

    sbn = []
    sbf = []
    pts = []
    sdchi2fid = bgstar['DELTACHI2_fid'] > 40
    #print(len(combpass),len(combpass[sdchi2fid]))
    sdchi2new = bgstar['DELTACHI2_new'] > 40

    for pt in range(0,10): 
        pts.append(pt)
        sp = bgstar['PETAL_LOC_fid'] == pt
        spbn = len(bgstar[sp&sdchi2new])/len(bgstar[sp])
        #print(len(combpass[sp&welg&sdchi2new]),len(combpass[sp&welg]))
        sbn.append(spbn)
        spbf = len(bgstar[sp&sdchi2fid])/len(bgstar[sp])
        sbf.append(spbf)

    plt.plot(pts,sbn,'^-',color='purple',label='f5')  
    plt.plot(pts,sbf,'^--',color='purple',label='everest')
    plt.title('BGS_ANY targets from '+args.survey)
    plt.xlabel('PETAL_LOC')
    plt.ylabel(r'fraction with $\Delta\chi^2>40$')
    plt.legend()
    plt.show()   


if args.prog == 'dark':
    #now specifically look at LRGs

    wlrg = (combpass[tarcol+'_fid'] & 1) > 0
    lrgtar = combpass[wlrg]

    def glrg(dat,ver):
        drz = (10**(3 - 3.5*dat['Z_'+ver]))
        mask_bad = (drz>30) & (dat['DELTACHI2_'+ver]<30)
        mask_bad |= (drz<30) & (dat['DELTACHI2_'+ver]<drz)
        mask_bad |= (dat['DELTACHI2_'+ver]<10)
        wz = dat['ZWARN_'+ver] == 0
        wz &= dat['Z_'+ver]<1.4
        wz &= (~mask_bad)

        wzwarn = wz#zmtlf['ZWARN'] == 0
        gzlrg = dat[wzwarn]
        return gzlrg

    glfid = glrg(lrgtar,'fid')
    glnew = glrg(lrgtar,'new')

    print('old LRG redshift success was '+str(len(glfid)/len(lrgtar)))
    print('new LRG redshift success was '+str(len(glnew)/len(lrgtar)))

    #do per petal
    sln = []
    slf = []
    for pt in range(0,10):
        sp = lrgtar['PETAL_LOC_fid'] == pt
        spf = glfid['PETAL_LOC_fid'] == pt
        spn = glnew['PETAL_LOC_fid'] == pt
        print('for petal '+str(pt))
        print('old LRG success rate was '+str(len(glfid[spf])/len(lrgtar[sp])))
        print('new LRG success rate is '+str(len(glnew[spn])/len(lrgtar[sp])))
        sln.append(len(glnew[spn])/len(lrgtar[sp]))
        slf.append(len(glfid[spf])/len(lrgtar[sp]))

    #plot QSO and ELG using deltachi2 > 25 threshold

    sdchi2fid = combpass['DELTACHI2_fid'] > 25
    #print(len(combpass),len(combpass[sdchi2fid]))
    sdchi2new = combpass['DELTACHI2_new'] > 25
    #print(len(combpass),len(combpass[sdchi2fid]),len(combpass[sdchi2new]))
    welg = (combpass[tarcol+'_fid'] & 2) > 0
    wqso = (combpass[tarcol+'_fid'] & 4) > 0
    sen = []
    sqn = []
    sef = []
    sqf = []
    pts = []

    for pt in range(0,10): 
        pts.append(pt)
        sp = combpass['PETAL_LOC_fid'] == pt
        spen = len(combpass[sp&welg&sdchi2new])/len(combpass[sp&welg])
        #print(len(combpass[sp&welg&sdchi2new]),len(combpass[sp&welg]))
        sen.append(spen)
        spef = len(combpass[sp&welg&sdchi2fid])/len(combpass[sp&welg])
        sef.append(spef)
        spqn = len(combpass[sp&wqso&sdchi2new])/len(combpass[sp&wqso])
        sqn.append(spqn)
        spqf = len(combpass[sp&wqso&sdchi2fid])/len(combpass[sp&wqso])
        sqf.append(spqf)

    plt.plot(pts,sln,'o-r')  
    plt.plot(pts,slf,'o--r')
    plt.plot(pts,sen,'o-b')  
    plt.plot(pts,sef,'o--b') 
    plt.plot(pts,sqn,'o-',color='orange')  
    plt.plot(pts,sqf,'o--',color='orange')
    plt.show()   



