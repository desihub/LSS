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
import LSS.main.cattools as ct
import LSS.common_tools as common

from LSS.globals import main

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--min_real", help="minimum number for mock realization",default=0,type=int)
parser.add_argument("--max_real", help="maximum (+1) for mock realization",default=25,type=int) 

parser.add_argument("--base_dir", help="base directory for input/output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data_version",help="version for redshifts",default='v0.6')
parser.add_argument("--mockcatver", default=None, help = "if not None, gets added to the output path")
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18) 
parser.add_argument("--imsys_zbin",help="if yes, do imaging systematic regressions in z bins",default='y')
parser.add_argument("--add_imsys_ran",help="add sysnet weights to randoms",default='n')
parser.add_argument("--par",help="whether to run in parallel",default='n')

parser.add_argument("--imsys_nside",help="healpix nside used for imaging systematic regressions",default=256,type=int)
parser.add_argument("--imsys_colname",help="column name for fiducial imaging systematics weight, if there is one (array of ones by default)",default=None)


args = parser.parse_args()
print(args)

tp = args.tracer
rm = int(args.minr)
rx = int(args.maxr)


datadir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'+args.verspec+'/LSScats/'+args.data_version+'/'
    

if tp[:3] == 'BGS' or tp == 'bright' or tp == 'MWS_ANY':
    prog = 'BRIGHT'

else:
    prog = 'DARK'

progl = prog.lower()

mainp = main(args.tracer,args.verspec,survey=args.survey)

zmin = mainp.zmin
zmax = mainp.zmax


fit_maps = mainp.fit_maps


tpstr = tp
if tp == 'BGS_BRIGHT-21.5':
    tpstr = 'BGS_BRIGHT'
nside = 256


if tp[:3] == 'ELG':
    if args.imsys_zbin == 'y':
        zrl = [(0.8,1.1),(1.1,1.6)]
    else:
        zrl = [(0.8,1.6)]
if tp[:3] == 'QSO':
    if args.imsys_zbin == 'y':
        zrl = [(0.8,1.3),(1.3,2.1)]#,(2.1,3.5)] 
    else:
        zrl = [(0.8,3.5)]   
if tp[:3] == 'LRG':
    if args.imsys_zbin == 'y':
        zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)] 
    else:
        zrl = [(0.4,1.1)]  
if tp == 'BGS_BRIGHT-21.5':
    zrl = [(0.1,0.4)]
elif tp[:3] == 'BGS':
    zrl = [(0.01,0.5)]
    zmin = 0.01
    zmax = 0.5    

def splitGC_wo(flroot,datran='.dat',rann=0):
    import LSS.common_tools as common
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    app = 'clustering'+datran+'.fits'
    if datran == '.ran':
        app = str(rann)+'_clustering'+datran+'.fits'

    fn = fitsio.read(flroot+app)
    #c = SkyCoord(fn['RA']* u.deg,fn['DEC']* u.deg,frame='icrs')
    #gc = c.transform_to('galactic')
    sel_ngc = common.splitGC(fn)#gc.b > 0
    outf_ngc = flroot+'NGC_'+app
    common.write_LSS(fn[sel_ngc],outf_ngc)
    outf_sgc = flroot+'SGC_'+app
    common.write_LSS(fn[~sel_ngc],outf_sgc)



debv = common.get_debv()
    


lssmapdirout = datadir+'/hpmaps/'

def get_imlin(realization):
    mockdir = args.base_dir+'mock'+str(realization)+'/'
    if args.mockcatver is not None:
        mockdir += args.mockcatver+'/'
    dirout = mockdir
    from LSS.imaging import densvar
    use_maps = fit_maps
       
    dat = Table(fitsio.read(os.path.join(dirout.replace('global','dvs_ro') , tp+'_clustering.dat.fits')))
    ranl = []
    for i in range(0,1):
        #rann = fitsio.read(os.path.join(dirout.replace('global','dvs_ro') , tp+'_NGC'+'_'+str(i)+'_clustering.ran.fits'), columns=['RA', 'DEC','WEIGHT','WEIGHT_FKP']) 
        #rans = fitsio.read(os.path.join(dirout.replace('global','dvs_ro') , tp+'_SGC'+'_'+str(i)+'_clustering.ran.fits'), columns=['RA', 'DEC','WEIGHT','WEIGHT_FKP']) 
        #ran = np.concatenate((rann,rans))
        #ran = common.addNS(Table(ran))
        ran = fitsio.read(os.path.join(dirout.replace('global','dvs_ro') , tp+'_'+str(i)+'_clustering.ran.fits'), columns=['RA', 'DEC','WEIGHT','WEIGHT_FKP','PHOTSYS'])
        ranl.append(ran)
    rands = np.concatenate(ranl)
    regl = ['N','S']


    syscol = 'WEIGHT_IMLIN'
    
    dat[syscol] = np.ones(len(dat))
    for reg in regl:
        pwf = lssmapdirout+'QSO_mapprops_healpix_nested_nside'+str(nside)+'_'+reg+'.fits'
        sys_tab = Table.read(pwf)
        cols = list(sys_tab.dtype.names)
        for col in cols:
            if 'DEPTH' in col:
                bnd = col.split('_')[-1]
                sys_tab[col] *= 10**(-0.4*common.ext_coeff[bnd]*sys_tab['EBV'])
        for ec in ['GR','RZ']:
            if 'EBV_DIFF_'+ec in fit_maps: 
                sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
        #seld = dat['PHOTSYS'] == reg
        selr = rands['PHOTSYS'] == reg
        print(zrl)
        for zr in zrl:
            zmin = zr[0]
            zmax = zr[1]
            
            print('getting weights for region '+reg+' and '+str(zmin)+'<z<'+str(zmax))
            wsysl = densvar.get_imweight(dat,rands[selr],zmin,zmax,reg,fit_maps,use_maps,sys_tab=sys_tab,zcol='Z',wtmd='wt_comp',figname=dirout+args.tracer+'_'+reg+'_'+str(zmin)+str(zmax)+'_linimsysfit.png')
            sel = wsysl != 1
            dat[syscol][sel] = wsysl[sel]
    fname = os.path.join(dirout, tp+'_clustering.dat.fits')
    common.write_LSS(dat,fname)
    flroot = os.path.join(dirout, tp+'_')
    splitGC_wo(flroot)







    if args.add_imsys_ran == 'y':
        regl = ['NGC','SGC']
        
        wtcol = 'WEIGHT_IMLIN'
        fb = dirout+tp
        fcdn = fitsio.read(fb.replace('global','dvs_ro')+'_NGC_clustering.dat.fits',columns=['TARGETID',wtcol])
        fcds = fitsio.read(fb.replace('global','dvs_ro')+'_SGC_clustering.dat.fits',columns=['TARGETID',wtcol])
        indata = Table(np.concatenate((fcdn,fcds)))
        indata.rename_column('TARGETID', 'TARGETID_DATA')

        def addrancol(rn):
            for reg in regl:
                fname = dirout+tp+'_'+reg+'_'+str(rn)+'_clustering.ran.fits'
                cd = Table(fitsio.read(fname.replace('global','dvs_ro')))
                cols2rem = [wtcol,wtcol+'_1',wtcol+'_2']
                for col in cols2rem:
                    if col in list(cd.dtype.names):
                        cd.remove_column(col)
                cd = join(cd,indata,keys=['TARGETID_DATA'],join_type='left')
                common.write_LSS(cd,fname)
    
        if args.par == 'n':
            for rn in range(rm,rx):
                addrancol(rn)
        if args.par == 'y':
            nproc = 9
            nran = rx-rm
            inds = np.arange(nran)
            from multiprocessing import Pool
            with Pool(processes=nproc) as pool:
                res = pool.map(addrancol, inds)

if __name__ == '__main__':
    from multiprocessing import Pool
    from desitarget.internal import sharedmem
    import sys
    inds = []
    for i in range(args.min_real,args.max_real):
        inds.append(i)
    pool = sharedmem.MapReduce()
    with pool:
        pool.map(get_imlin,inds)#,reduce=reduce)

