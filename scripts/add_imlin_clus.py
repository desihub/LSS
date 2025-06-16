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

import logging
logname = 'mkCat'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


import LSS.main.cattools as ct
import LSS.common_tools as common

from LSS.globals import main
#except:

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='DA2')
parser.add_argument("--verspec",help="version for redshifts",default='loa-v1')
parser.add_argument("--usemaps", help="the list of maps to use; defaults to what is set by globals", type=str, nargs='*',default=None)
parser.add_argument("--extra_clus_dir", help="an optional extra layer of directory structure for clustering catalog",default='fNL/')

parser.add_argument("--relax_zbounds", help="whether or not to use less restricted redshift bounds",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18,type=int) 

parser.add_argument("--imsys_zbin",help="if yes, do imaging systematic regressions in z bins",default='y')
parser.add_argument("--imsys_finezbin",help="if yes, do imaging systematic regressions in dz=0.1 bins",default='n')
parser.add_argument("--imsys_1zbin",help="if yes, do imaging systematic regressions in just 1 z bin",default='n')
#parser.add_argument("--imsys_clus_fb",help="perform linear weight fits in fine redshift bins",default='n')
parser.add_argument("--imsys_clus",help="add weights for imaging systematics using eboss method, applied to clustering catalogs?",default='n')

parser.add_argument("--syscol",help="name for new systematic column (automatically determined based on other choices if None)",default=None)

parser.add_argument("--imsys_clus_ran",help="add weights for imaging systematics using eboss method, applied to clustering catalogs, to randoms?",default='n')

#parser.add_argument("--imsys_clus_fb_ran",help="add linear weight fits in fine redshift bins to randoms",default='n')

parser.add_argument("--replace_syscol",help="whether to replace any existing weight_sys with new",default='n')
parser.add_argument("--add_syscol2blind",help="whether to add the new weight column to the blinded catalogs",default='n')

parser.add_argument("--nran4imsys",help="number of random files to using for linear regression",default=10,type=int)

parser.add_argument("--par", help="run different random number in parallel?",default='y')


args = parser.parse_args()
common.printlog(str(args),logger)

type = args.type
tp = type
basedir = args.basedir
version = args.version
specrel = args.verspec
rm = int(args.minr)
rx = int(args.maxr)



    

if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    prog = 'BRIGHT'

else:
    prog = 'DARK'

progl = prog.lower()

mainp = main(args.type,args.verspec,survey=args.survey,relax_zbounds=args.relax_zbounds)
mdir = mainp.mdir+progl+'/' #location of ledgers
tdir = mainp.tdir+progl+'/' #location of targets
#mtld = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied

tsnrcut = mainp.tsnrcut
dchi2 = mainp.dchi2
tsnrcol = mainp.tsnrcol        
zmin = mainp.zmin
zmax = mainp.zmax

maindir = basedir +'/'+args.survey+'/LSS/'

ldirspec = maindir+specrel+'/'
    
#dirout = ldirspec+'LSScats/'+version+'/'
#dirout = '/pscratch/sd/x/xychen/imsys_tests/unblinded/'
dirout = '/pscratch/sd/x/xychen/imsys_tests/altmtl_mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl0/kibo-v1/mock0/LSScats/'





dirin = dirout
#lssmapdirout = dirout+'/hpmaps/'
lssmapdirout = '/pscratch/sd/x/xychen/imsys_tests/unblinded/'+'hpmaps/'

if args.syscol is None:
    if args.imsys_zbin == 'y':
        syscol = 'WEIGHT_IMLIN'
        #if args.usemaps[0] == 'all':
        #    syscol += '_ALL'
        #if args.usemaps[0] == 'allebv':
        #    syscol += '_ALLEBV'
    if args.imsys_1zbin == 'y':
        syscol = 'WEIGHT_IMLIN_1ZBIN'
        #if args.usemaps[0] == 'all':
        #    syscol += '_ALL'
        #if args.usemaps[0] == 'allebv':
        #    syscol += '_ALLEBV'

    if args.imsys_finezbin == 'y':
        syscol = 'WEIGHT_IMLIN_FINEZBIN'
        #if args.usemaps[0] == 'allebv':
        #    syscol += '_ALLEBV'

else:
    syscol = args.syscol


if args.usemaps == None:
    fit_maps = mainp.fit_maps
    if args.imsys_finezbin == 'y':
        mainp.fit_maps_all
elif args.usemaps[0] == 'all': 
    fit_maps = mainp.fit_maps_all
    syscol += '_ALL'
elif args.usemaps[0] == 'allebv':
    fit_maps = mainp.fit_maps_allebv
    syscol += '_ALLEBV'
elif args.usemaps[0] == 'allebvcmb':
    fit_maps = mainp.fit_maps_allebvcmb
    syscol += '_ALLEBVCMB'

else:
    fit_maps = [mapn for mapn in args.usemaps]

common.printlog('using '+str(fit_maps),logger)

zl = (zmin,zmax)

tracer_clus = args.type
tpstr = args.type
if 'BGS_BRIGHT' in tracer_clus:
    tpstr = 'BGS_BRIGHT'
if 'LRG' in tracer_clus:
    tpstr = 'LRG'
nside = 256
inds = np.arange(rm,rx)

if type[:3] == 'ELG':
    if args.imsys_zbin == 'y':
        zrl = [(0.8,1.1),(1.1,1.6)]
    elif args.imsys_1zbin == 'y':
        zrl = [(0.8,1.6)]
    elif args.imsys_finezbin == 'y':
        imsys_clus_fb = 'y'
    zsysmin = 0.8
    zsysmax = 1.6


if type[:3] == 'QSO':
    if args.imsys_zbin == 'y':
        zrl = [(0.8,1.3),(1.3,2.1),(2.1,3.5)] 
    elif args.imsys_1zbin == 'y':
        zrl = [(0.8,3.5)]   
    elif args.imsys_finezbin == 'y':
        imsys_clus_fb = 'y'
    zsysmin = 0.8
    zsysmax = 3.5
if type[:3] == 'LRG':
    if args.imsys_zbin == 'y':
        zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)] 
    elif args.imsys_1zbin == 'y':
        zrl = [(0.4,1.1)]
    elif args.imsys_finezbin == 'y':
        imsys_clus_fb = 'y'
    zsysmin = 0.4
    zsysmax = 1.1
    if args.relax_zbounds == 'y':
        zsysmax = 1.2
        zsysmin = 0.3      

if 'BGS_BRIGHT-' in type:
    zrl = [(0.1,0.4)]
elif type[:3] == 'BGS':
    zrl = [(0.01,0.5)]
    zmin = 0.01
    zmax = 0.5    



common.printlog('the added weight column will be '+syscol,logger)

debv = common.get_debv()
zcmb = common.mk_zcmbmap()
sky_g,sky_r,sky_z = common.get_skyres()


if args.imsys_clus == 'y':
    from LSS.imaging import densvar
    
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_clustering.dat.fits')
    dat_ngc = Table(fitsio.read(fname))
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_clustering.dat.fits')
    dat_sgc = Table(fitsio.read(fname))
    dat = vstack([dat_sgc,dat_ngc])
    foutname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_clustering.dat.fits')
    ranl = []
    for i in range(0,args.nran4imsys):#int(args.maxr)):
        ran = fitsio.read(os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_'+str(i)+'_clustering.ran.fits')) 
        ranl.append(ran)
        ran = fitsio.read(os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_'+str(i)+'_clustering.ran.fits')) 
        ranl.append(ran)

    rands = np.concatenate(ranl)
    #syscol = 'WEIGHT_IMLIN_CLUS'
    regl = ['S','N']
    if args.type == 'QSO':
        regl = ['DES','SnotDES','N']
    dat[syscol] = np.ones(len(dat))
    for reg in regl:
        regu = reg
        if reg == 'DES' or reg == 'SnotDES':
            regu = 'S'
       # pwf = lssmapdirout+tpstr+'_mapprops_healpix_nested_nside'+str(nside)+'_'+regu+'.fits'
        pwf = lssmapdirout+'QSO_mapprops_healpix_nested_nside'+str(nside)+'_'+regu+'.fits'
        sys_tab = Table.read(pwf)
        cols = list(sys_tab.dtype.names)
        for col in cols:
            if 'DEPTH' in col:
                bnd = col.split('_')[-1]
                sys_tab[col] *= 10**(-0.4*common.ext_coeff[bnd]*sys_tab['EBV'])
        for ec in ['GR','RZ']:
            if 'EBV_DIFF_'+ec in fit_maps: 
                sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
        if 'ZCMB' in fit_maps:
            sys_tab['ZCMB'] = zcmb
        #seld = dat['PHOTSYS'] == reg
        selr = rands['PHOTSYS'] == reg

        def _add_sysweight(zm,zx):
            common.printlog('getting weights for region '+reg+' and '+str(zm)+'<z<'+str(zx),logger)
            wsysl = densvar.get_imweight(dat,rands,zm,zx,reg,fitmapsbin,use_maps,sys_tab=sys_tab,zcol='Z',modoutname = dirout+args.extra_clus_dir+tracer_clus+'_'+reg+'_'+str(zm)+str(zx)+'_linfitparam.txt',figname=dirout+args.extra_clus_dir+tracer_clus+'_'+reg+'_'+str(zm)+str(zx)+'_linclusimsysfit.png',wtmd='clus',logger=logger)
            sel = wsysl != 1
            dat[syscol][sel] = wsysl[sel]
       
        if args.imsys_finezbin == 'y':
            dz = 0.1
            zm = zsysmin
            zx = zm + dz
            fitmapsbin = fit_maps
            while zm < zsysmax:
                zx = zm + dz
                zx = round(zx,1)
                #this is now controlled above
                #if type == 'LRG':
                #    fitmapsbin = mainp.fit_maps_all
                #else:
                #    fitmapsbin = fit_maps
                use_maps = fitmapsbin
                _add_sysweight(zm,zx)
                zm = zx
        elif args.imsys_1zbin == 'y':
            zm = zmin
            zx = zmax
            if type == 'LRG':
                fitmapsbin = mainp.fit_maps_all
            else:
                fitmapsbin = fit_maps
            use_maps = fitmapsbin
            _add_sysweight(zm,zx)
            
        elif args.imsys_zbin == 'y':
       
            for zr in zrl:
                zm = zr[0]
                zx = zr[1]
                if type == 'LRG':
                    if reg == 'N':
                        fitmapsbin = fit_maps
                    else:
                        if zmax == 0.6:
                            fitmapsbin = mainp.fit_maps46s
                        if zmax == 0.8:
                            fitmapsbin = mainp.fit_maps68s
                        if zmax == 1.1:
                            fitmapsbin = mainp.fit_maps81s
                else:
                    fitmapsbin = fit_maps
                use_maps = fitmapsbin
                _add_sysweight(zm,zx)
        else:
            sys.exit('no valid z binning choice in arguments, exiting...')
    #attach data to NGC/SGC catalogs, write those out
    dat.keep_columns(['TARGETID',syscol])
    if syscol in list(dat_ngc.colnames):
        dat_ngc.remove_column(syscol)
    dat_ngc = join(dat_ngc,dat,keys=['TARGETID'])    
    if args.replace_syscol == 'y':
        dat_ngc['WEIGHT'] /= dat_ngc['WEIGHT_SYS']
        dat_ngc['WEIGHT_SYS'] = dat_ngc[syscol]
        dat_ngc['WEIGHT'] *= dat_ngc['WEIGHT_SYS']
    common.write_LSS_scratchcp(dat_ngc,os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_clustering.dat.fits'),logger=logger)
    if syscol in list(dat_sgc.colnames):
        dat_sgc.remove_column(syscol)
    dat_sgc = join(dat_sgc,dat,keys=['TARGETID'])
    if args.replace_syscol == 'y':
        dat_sgc['WEIGHT'] /= dat_sgc['WEIGHT_SYS']
        dat_sgc['WEIGHT_SYS'] = dat_sgc[syscol]
        dat_sgc['WEIGHT'] *= dat_sgc['WEIGHT_SYS']
    common.write_LSS_scratchcp(dat_sgc,os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_clustering.dat.fits'),logger=logger)


if args.imsys_clus_ran == 'y':
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_NGC_clustering.dat.fits')
    dat_ngc = Table(fitsio.read(fname,columns=['TARGETID',syscol]))
    fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_SGC_clustering.dat.fits')
    dat_sgc = Table(fitsio.read(fname,columns=['TARGETID',syscol]))
    dat = vstack([dat_sgc,dat_ngc])
    dat.rename_column('TARGETID','TARGETID_DATA')
    regl = ['NGC','SGC']
    syscolr = syscol
    #if args.replace_syscol == 'y':
    #    syscolr = 'WEIGHT_SYS'
    def _add2ran(rann):
        for reg in regl:
            ran_fn = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_'+reg+'_'+str(rann)+'_clustering.ran.fits')
            ran = Table(fitsio.read(ran_fn))
            if syscolr in ran.colnames:
                ran.remove_column(syscolr)
            ran = join(ran,dat,keys=['TARGETID_DATA'])
            if args.replace_syscol == 'y':
                ran['WEIGHT'] /= ran['WEIGHT_SYS']
                ran['WEIGHT_SYS'] = ran[syscolr]
                ran['WEIGHT'] *= ran['WEIGHT_SYS']
            common.write_LSS_scratchcp(ran,ran_fn,logger=logger)

    if args.par == 'y':
        from multiprocessing import Pool
        with Pool() as pool:
            res = pool.map(_add2ran, inds)
    else:
        for rn in inds:#range(rm,rx):
             _add2ran(rn)
            
if args.add_syscol2blind == 'y':
    syscolr = syscol
    #if args.replace_syscol == 'y':
    #    syscolr = 'WEIGHT_SYS'

    dats = []
    for reg in ['NGC','SGC']:
        fname = os.path.join(dirout+args.extra_clus_dir, tracer_clus+'_'+reg+'_clustering.dat.fits')
        dati = Table(fitsio.read(fname,columns=['TARGETID',syscol]))
        dats.append(dati)
        fname_blind = os.path.join(dirout+args.extra_clus_dir+'/blinded/', tracer_clus+'_'+reg+'_clustering.dat.fits')
        dat_blind = Table(fitsio.read(fname_blind))
        if syscol in list(dat_blind.colnames):
            dat_blind.remove_column(syscol)

        dat_blind = join(dat_blind,dati,keys=['TARGETID'])
        if args.replace_syscol == 'y':
            dat_blind['WEIGHT'] /= dat_blind['WEIGHT_SYS']
            dat_blind['WEIGHT_SYS'] = dat_blind[syscol]
            dat_blind['WEIGHT'] *= dat_blind['WEIGHT_SYS']

        common.write_LSS_scratchcp(dat_blind,fname_blind,logger=logger)
    dat = vstack(dats)
    dat.rename_column('TARGETID','TARGETID_DATA')
    regl = ['NGC','SGC']
    def _add2ranblind(rann):
        for reg in regl:
            ran_fn = os.path.join(dirout+args.extra_clus_dir+'/blinded/', tracer_clus+'_'+reg+'_'+str(rann)+'_clustering.ran.fits')
            ran = Table(fitsio.read(ran_fn))
            if syscolr in ran.colnames:
                ran.remove_column(syscolr)
            ran = join(ran,dat,keys=['TARGETID_DATA'])
            if args.replace_syscol == 'y':
                ran['WEIGHT'] /= ran['WEIGHT_SYS']
                ran['WEIGHT_SYS'] = ran[syscolr]
                ran['WEIGHT'] *= ran['WEIGHT_SYS']

            common.write_LSS_scratchcp(ran,ran_fn,logger=logger)

    if args.par == 'y':
        from multiprocessing import Pool
        with Pool() as pool:
            res = pool.map(_add2ranblind, inds)
    else:
        for rn in inds:#range(rm,rx):
             _add2ranblind(rn)
    


