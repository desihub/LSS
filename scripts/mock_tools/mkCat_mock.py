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
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desitarget import targetmask
from desitarget.internal import sharedmem
from desimodel.footprint import is_point_in_desi
from desitarget import targetmask

import LSS.main.cattools as ct
import LSS.common_tools as common
import LSS.mocktools as mocktools
#import LSS.mkCat_singletile.fa4lsscat as fa
#from LSS.globals import main

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--mockver", help="type of mock to use",default='ab_firstgen')
parser.add_argument("--mockmin", help="number for the realization",default=1,type=int)
parser.add_argument("--mockmax", help="number for the realization",default=2,type=int)
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/main/mocks/')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='DA02')
#parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--combd", help="combine the data tiles together",default='n')
parser.add_argument("--combr", help="combine the random tiles together",default='n')
parser.add_argument("--combdr", help="combine the random tiles info together with the assignment info",default='n')
parser.add_argument("--fulld", help="make the 'full' data files ",default='n')
parser.add_argument("--fullr", help="make the random files associated with the full data files",default='n')
parser.add_argument("--add_gtl", help="whether to get the list of good tileloc from observed data",default='y')

parser.add_argument("--add_veto", help="add veto column to the full files",default='n')
parser.add_argument("--apply_veto", help="apply vetos to the full files",default='n')
parser.add_argument("--mkclusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--mkclusdat", help="make the data clustering files; these are cut to a small subset of columns",default='n')

parser.add_argument("--mkclusran_allpot", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--mkclusdat_allpot", help="make the data clustering files; these are cut to a small subset of columns",default='n')

parser.add_argument("--mkclusran_tiles", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--mkclusdat_tiles", help="make the data clustering files; these are cut to a small subset of columns",default='n')


parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=1,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 40 are available (use parallel script for all)",default=2,type=int) 
parser.add_argument("--par", help="run different random number in parallel?",default='n')

parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
parser.add_argument("--newspec",help="if y, merge in redshift info even if no new tiles",default='n')

args = parser.parse_args()
print(args)

rm = int(args.minr)
rx = int(args.maxr)
par = True
if args.par == 'n':
    par = False

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'

tracer = args.tracer
survey = args.survey

if tracer[:3] == 'BGS' or tracer == 'bright' or tracer == 'MWS_ANY':
    pr = 'BRIGHT'
    pdir = 'bright'
else:
    pr = 'DARK'
    pdir = 'dark'

pd = pdir

if args.mockver == 'ab_firstgen':
    mockdir = 'FirstGenMocks/AbacusSummit/'
    randens = 10460.

maindir = args.base_output +mockdir+args.survey+'/'

tiles = fitsio.read( '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/tiles-'+pr+'.fits')



gtl = None
if args.add_gtl == 'y':
    
    if args.survey == 'DA02':
        datarel = 'guadalupe'
    datadir = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+datarel+'/'    
    specdat = ct.get_specdat(datadir,pdir,datarel)
    gtl = np.unique(specdat['TILELOCID'])

def docat(mocknum,rannum):

    lssdir = maindir+'mock'+str(mocknum)+'/'
    if not os.path.exists(lssdir):
        os.mkdir(lssdir)
        print('made '+lssdir)

    dirout = lssdir+'LSScats/'
    if not os.path.exists(dirout):
        os.mkdir(dirout)
        print('made '+dirout)


    if args.tracer != 'dark' and args.tracer != 'bright':
        if args.tracer == 'BGS_BRIGHT':
            bit = targetmask.bgs_mask[args.tracer]
            desitarg='BGS_TARGET'
        else:
            bit = targetmask.desi_mask[args.tracer]
            desitarg='DESI_TARGET'


    if args.combr == 'y' and mocknum == 1:
        fbadir = maindir+'random_fba'+str(rannum)
        outdir = fbadir
        tarf = fbadir+'/targs.fits'
        common.combtiles_pa_wdup(tiles,fbadir,outdir,tarf,addcols=['TARGETID','RA','DEC'],fba=True,tp=pdir)

    if args.combd == 'y' and rannum == 1:
        fbadir = maindir+'fba'+str(mocknum)
        outdir = fbadir
        tarf = fbadir+'/targs.fits'
        asn = common.combtiles_assign_wdup(tiles,fbadir,outdir,tarf,tp=pdir)
        #if using alt MTL that should have ZWARN_MTL, put that in here
        asn['ZWARN_MTL'] = np.copy(asn['ZWARN'])
        pa = common.combtiles_pa_wdup(tiles,fbadir,outdir,tarf,addcols=['TARGETID','RA','DEC'],fba=True,tp=pdir,ran='dat')

        pa['TILELOCID'] = 10000*pa['TILEID'] +pa['LOCATION']
        tj = join(pa,asn,keys=['TARGETID','LOCATION','TILEID'],join_type='left')
        outfs = lssdir+'datcomb_'+pdir+'_tarspecwdup_zdone.fits'
        tj.write(outfs,format='fits', overwrite=True)
        tc = ct.count_tiles_better('dat',pdir,specrel='',survey=args.survey,indir=lssdir,gtl=gtl) 
        outtc =  lssdir+'Alltiles_'+pdir+'_tilelocs.dat.fits'
        tc.write(outtc,format='fits', overwrite=True)
    
    if args.combdr == 'y':
        fbadir_data = maindir+'fba'+str(mocknum)
        fbadir_ran = maindir+'random_fba'+str(rannum)
        specf = Table(fitsio.read(fbadir_data+'/datcomb_'+pdir+'assignwdup.fits'))
        specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
        specf.remove_columns(['TARGETID'])
        fgu = Table(fitsio.read(fbadir_ran+'/rancomb_'+pdir+'wdup.fits'))
        print(len(fgu))
        fgu = join(fgu,specf,keys=['LOCATION','TILEID'],join_type='left')
        print(len(fgu))
        print(fgu.dtype.names)
        fgu.sort('TARGETID')
        outf = lssdir+'/rancomb_'+str(rannum)+pdir+'wdupspec_zdone.fits'
        print(outf)
        fgu.write(outf,format='fits', overwrite=True)
        tc = ct.count_tiles_better('ran',pdir,rannum,specrel='',survey=args.survey,indir=lssdir,gtl=gtl)
        tc.write(lssdir+'/rancomb_'+str(rannum)+pdir+'_Alltilelocinfo.fits',format='fits', overwrite=True)

    specver = 'mock'    
    imbits = []    
    if args.fulld == 'y':
        
        ftar = None
        dz = lssdir+'datcomb_'+pdir+'_tarspecwdup_zdone.fits'
        tlf = lssdir+'Alltiles_'+pdir+'_tilelocs.dat.fits'
        ct.mkfulldat(dz,imbits,ftar,args.tracer,bit,dirout+args.tracer+notqso+'_full_noveto.dat.fits',tlf,desitarg=desitarg,specver=specver,notqso=notqso,gtl_all=gtl)

    maxp = 3400
    pthresh = 3000
    zmin = 0.8
    zmax = 3.5
    if args.tracer[:3] == 'LRG':# or notqso == 'notqso':
        maxp = 3200
        zmin = 0.4
        zmax = 1.1
    if args.tracer[:3] == 'ELG':
        maxp = 3000
        zmin = 0.8
        zmax = 1.6
    if args.tracer[:3] == 'BGS':
        maxp = 2100
        pthresh = 2000
        zmin = 0.1
        zmax = 0.5


        
    if args.fullr == 'y':
        zf = lssdir+'datcomb_'+pdir+'_tarspecwdup_zdone.fits'
        dz = Table.read(zf) 

        if gtl is None:
            specdat = common.cut_specdat(dz)
            gtlf = np.unique(specdat['TILELOCID'])  
            del specdat     
        else:
            gtlf = gtl
        wg = np.isin(dz['TILELOCID'],gtlf)
        dz = dz[wg]
        wtype = ((dz[desitarg] & bit) > 0)
        if notqso == 'notqso':
            wtype &= ((dz[desitarg] & 4) == 0)
        dz = dz[wtype]
        lznp,tlid_full = common.find_znotposs_tloc(dz,priority_thresh=pthresh)
        del dz
        

        outf = dirout+args.tracer+notqso+'_'+str(rannum)+'_full_noveto.ran.fits'
        
        ct.mkfullran(gtlf,lznp,lssdir,rannum,imbits,outf,args.tracer,pdir,notqso=notqso,maxp=maxp,tlid_full=tlid_full)
        


    if args.apply_veto == 'y':
        print('applying vetos to mock '+str(mocknum))
        fin = dirout+args.tracer+notqso+'_full_noveto.dat.fits'
        fout = dirout+args.tracer+notqso+'_full.dat.fits'
        common.apply_veto(fin,fout,ebits=None,zmask=False,maxp=maxp)
        print('applying vetos to random '+str(rannum))
        fin = dirout+args.tracer+notqso+'_'+str(rannum)+'_full_noveto.ran.fits'
        fout = dirout+args.tracer+notqso+'_'+str(rannum)+'_full.ran.fits'
        common.apply_veto(fin,fout,ebits=None,zmask=False,maxp=maxp)

        #print('random veto '+str(ii)+' done')

    regl = ['_N','_S']    
    
    #needs to happen before randoms so randoms can get z and weights
    
    if args.mkclusdat == 'y':
        ct.mkclusdat(dirout+args.tracer+notqso,tp=args.tracer,dchi2=None,tsnrcut=0,zmin=zmin,zmax=zmax)#,ntilecut=ntile)


    if args.mkclusran == 'y':
        rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL']
        tsnrcol = 'TSNR2_ELG'
        if args.tracer[:3] == 'BGS':
            tsnrcol = 'TSNR2_BGS'
        ct.mkclusran(dirout+args.tracer+notqso+'_',dirout+args.tracer+notqso+'_',rannum,rcols=rcols,tsnrcut=0,tsnrcol=tsnrcol)#,ntilecut=ntile,ccut=ccut)
        #for clustering, make rannum start from 0
        for reg in regl:
            ranf = dirout+args.tracer+notqso+reg+'_'+str(rannum)+'_clustering.ran.fits'
            ranfm = dirout+args.tracer+notqso+reg+'_'+str(rannum-1)+'_clustering.ran.fits'
            os.system('mv '+ranf+' '+ranfm)

    nztl = ['']
    if args.mkclusdat_allpot == 'y':
        fbadir = maindir+'fba'+str(mocknum)
        tarf = fbadir+'/targs.fits'

        ztab = Table(fitsio.read(tarf,columns=['TARGETID','RSDZ']))
        ztab.rename_column('RSDZ', 'Z')
        mocktools.mkclusdat_allpot(dirout+args.tracer+notqso,ztab,tp=args.tracer,dchi2=None,tsnrcut=0,zmin=zmin,zmax=zmax)#,ntilecut=ntile)

    if args.mkclusran_allpot == 'y':
        nztl.append('_complete')
        rcols=['Z','WEIGHT']
        tsnrcol = 'TSNR2_ELG'
        if args.tracer[:3] == 'BGS':
            tsnrcol = 'TSNR2_BGS'
        ct.mkclusran(dirout+args.tracer+notqso+'_',dirout+args.tracer+notqso+'_complete'+'_',rannum,rcols=rcols,tsnrcut=0,tsnrcol=tsnrcol)#,ntilecut=ntile,ccut=ccut)
        #for clustering, make rannum start from 0
        for reg in regl:
            ranf = dirout+args.tracer+notqso+'_complete'+reg+'_'+str(rannum)+'_clustering.ran.fits'
            ranfm = dirout+args.tracer+notqso+'_complete'+reg+'_'+str(rannum-1)+'_clustering.ran.fits'
            os.system('mv '+ranf+' '+ranfm)

    if args.mkclusdat_tiles == 'y':
        fbadir = maindir+'fba'+str(mocknum)
        tarf = fbadir+'/targs.fits'

        ztab = Table(fitsio.read(tarf,columns=['RA','DEC','RSDZ','DESI_TARGET']))
        ztab.rename_column('RSDZ', 'Z')
        mocktools.mkclusdat_tiles(dirout+args.tracer+notqso,ztab,bit,zmin=zmin,zmax=zmax)#,ntilecut=ntile)

    if args.mkclusran_tiles == 'y':
        nztl.append('_tiles')
        fbadir = maindir+'random_fba'+str(rannum)
        tarf = fbadir+'/targs.fits'
        ffc = Table(fitsio.read(tarf,columns=['RA','DEC']))
        rcols=['Z','WEIGHT']
        mocktools.mkclusran_tiles(ffc,dirout+args.tracer+notqso+'_tiles'+'_',rannum,rcols=rcols)#,ntilecut=ntile,ccut=ccut)
        #for clustering, make rannum start from 0
        for reg in regl:
            ranf = dirout+args.tracer+notqso+'_tiles'+reg+'_'+str(rannum)+'_clustering.ran.fits'
            ranfm = dirout+args.tracer+notqso+'_tiles'+reg+'_'+str(rannum-1)+'_clustering.ran.fits'
            os.system('mv '+ranf+' '+ranfm)


    if args.nz == 'y':
    
        if args.tracer == 'QSO':
            dz = 0.05
            P0 = 6000
        
        else:    
            dz = 0.02
    
        if args.tracer[:3] == 'LRG':
            P0 = 10000
        if args.tracer[:3] == 'ELG':
            P0 = 4000
        if args.tracer[:3] == 'BGS':
            P0 = 7000
    
        for zt in nztl:
            for reg in regl:
                fb = dirout+args.tracer+notqso+zt+reg
                fcr = fb+'_0_clustering.ran.fits'
                fcd = fb+'_clustering.dat.fits'
                fout = fb+'_nz.txt'
                common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax,randens=randens)
                nran = rannum
                ranmin = rannum-1
                common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0,ranmin=ranmin,nran=nran)


    #print('done with random '+str(ii))
    return True
        #ct.mkclusran(dirout+type+'Alltiles_',ii,zmask=zma)
    #logf.write('ran mkclusran\n')
    #print('ran mkclusran\n')
    
if __name__ == '__main__':
    rx = args.maxr
    rm = args.minr
    mockmin = args.mockmin
    mockmax = args.mockmax
    if args.par == 'y':
        from multiprocessing import Pool
        from desitarget.internal import sharedmem
        
        N = rx-rm+1
        inds = []
        for i in range(rm,rx):
            inds.append(i)
        pool = sharedmem.MapReduce(np=N)
        with pool:
        
            def reduce( r):
                print('chunk done')
                return r
            pool.map(prep,inds,reduce=reduce)

        #p.map(doran,inds)
    else:
        for mn in range(mockmin,mockmax):
            for i in range(rm,rx):
                print('processing mock '+str(mn)+' and random '+str(i))
                docat(mn,i)
