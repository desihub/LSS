#get unblinded clustering catalogs, in all of N, S, NGC, SGC regions
#python scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --survey Y1 --verspec iron --version v0.6 --clusd y --clusran y --NStoGC y --nz y --par n --resamp y

##/pscratch/sd/a/acarnero/codes/LSS/Sandbox/LSSpipe_Y1.txt

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
parser.add_argument("--mockver", help="type of mock to use",default='ab_firstgen')
parser.add_argument("--famd", help="whether to use the fiberassign split into passes",default='passes')

parser.add_argument("--mockmin", help="number for the realization",default=1,type=int)
parser.add_argument("--mockmax", help="number for the realization",default=2,type=int)
parser.add_argument("--base_output", help="base directory for output",default='/pscratch/sd/a/acarnero/SecondGen/')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='DA02')
parser.add_argument("--specdata", help="mountain range for spec prod",default='himalayas')
parser.add_argument("--combd", help="combine the data tiles together",default='n')
parser.add_argument("--combr", help="combine the random tiles together",default='n')
parser.add_argument("--combdr", help="combine the random tiles info together with the assignment info",default='n')
parser.add_argument("--countran", help="count instances of focal plane locations for randoms",default='n')
parser.add_argument("--fulld", help="make the 'full' data files ",default='n')
parser.add_argument("--fullr", help="make the random files associated with the full data files",default='n')
parser.add_argument("--add_gtl", help="whether to get the list of good tileloc from observed data",default='y')
parser.add_argument("--mkHPmaps", help="make healpix maps for imaging properties using sample randoms",default='n')
parser.add_argument("--add_veto", help="add veto column to the full files",default='n')
parser.add_argument("--apply_veto", help="apply vetos to the full files",default='n')
parser.add_argument("--apply_veto_ran", help="apply vetos to the full files",default='n')
parser.add_argument("--mkclusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--mkclusdat", help="make the data clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--apply_map_veto", help="apply vetos to data and randoms based on values in healpix maps",default='n')
parser.add_argument("--mkclusran_allpot", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--mkclusdat_allpot", help="make the data clustering files; these are cut to a small subset of columns",default='n')

parser.add_argument("--mkclusran_tiles", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--mkclusdat_tiles", help="make the data clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--FKPfull", help="add FKP weights to full catalogs",default='n')
parser.add_argument("--split_GC",help='whether to combine N/S and then split NGC/SGC',default='n')

parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=1,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 40 are available (use parallel script for all)",default=2,type=int) 
parser.add_argument("--par", help="run different random number in parallel?",default='n')

parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
parser.add_argument("--newspec",help="if y, merge in redshift info even if no new tiles",default='n')
parser.add_argument("--equal_data_dens", help="if y, make mock n(z) equal data n(z)", default = 'n')
parser.add_argument("--nran_clus_data", help="number of random catalogues to use for clustering data", default = 4)
parser.add_argument("--use_map_veto", help="Tag for extraveto added in name, for example, _HPmapcut", default = '')
parser.add_argument("--resamp",help="resample radial info for different selection function regions",default='n')
parser.add_argument("--getFKP", help="calculate n(z) and FKP weights on final clustering catalogs", default='n')

#--use_map_veto _HPmapcut

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

randens = 10460. #the number density of randoms in the 1st gen file getting used
if args.mockver == 'ab_firstgen':
    mockdir = 'FirstGenMocks/AbacusSummit/'
    mockz = 'RSDZ'
    maindir = args.base_output +mockdir+args.survey+'/'

if args.mockver == 'EZ_3gpc1year':
    mockdir = 'FA_EZ_1year/fiberassign_EZ_3gpc/'    
    mockz = 'TRUEZ'
    maindir = args.base_output +mockdir+args.survey+'/'

if args.mockver == 'ab_secondgen':
    maindir = args.base_output
    mockz = 'RSDZ'
    args.famd = 'ab_secondgen'

if args.survey == 'MVMY1':
    tile_fn = '/global/cfs/cdirs/desi/users/FA_EZ_1year/fiberassign_EZ_3gpc/fba001/inputs/tiles.fits'
else:
    tile_fn = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/tiles-'+pr+'.fits'


tiles = fitsio.read(tile_fn)



gtl = None
if args.add_gtl == 'y':
    print('adding gtl')
    datarel = args.specdata
    if args.survey == 'DA02':
        datarel = 'guadalupe'
    datadir = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+datarel+'/'   
    specdat = ct.get_specdat(datadir,pdir,datarel,badfib= main(args.tracer, args.specdata, survey=args.survey).badfib)
    tlocid = 10000*specdat['TILEID'] +specdat['LOCATION']
    gtl = np.unique(tlocid)#np.unique(specdat['TILELOCID'])
    del specdat

def docat(mocknum,rannum):

    lssdir = os.path.join(maindir, 'mock'+str(mocknum))
    if not os.path.exists(lssdir):
        os.mkdir(lssdir)
        print('made '+lssdir)

    dirout = os.path.join(lssdir, 'LSScats')
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
        
        tarf = fbadir+'/targs.fits'
        if args.famd == 'passes':
            fbadir = maindir+'/ran'+str(rannum)+'_'+pdir+'/faruns/'
            tarf = maindir+'/ran'+str(rannum)+'_'+pdir+'/inputs/targ.fits'

        outdir = fbadir
        common.combtiles_pa_wdup(tiles,fbadir,outdir,tarf,addcols=['TARGETID','RA','DEC'],fba=True,tp=pdir)

    if args.combd == 'y' and rannum == 1:
        fbadir = os.path.join(maindir,'fba'+str(mocknum))
        outdir = fbadir
        if not os.path.exists(outdir):
            os.mkdir(outdir)
            print('made '+outdir)

        if args.survey == 'MVMY1':
            tarf = '/global/cfs/cdirs/desi/users/FA_EZ_1year/fiberassign_EZ_3gpc/fba'+str(mocknum).zfill(3)+'/inputs/targ.fits'
            indir = '/global/cfs/cdirs/desi/users/FA_EZ_1year/fiberassign_EZ_3gpc/fba'+str(mocknum).zfill(3)+'/'
            asn = mocktools.combtiles_assign_wdup_7pass(indir,outdir,tarf,tp=pdir)
            asn['ZWARN_MTL'] = np.copy(asn['ZWARN'])
            pa = mocktools.combtiles_pa_wdup_7pass(indir,outdir,tarf,addcols=['TARGETID','RA','DEC'],fba=True,tp=pdir,ran='dat',dtar='SV3_')
        elif args.famd == 'passes':
            fbadir = maindir+'/multipass_mock'+str(mocknum)+'_'+pdir+'/faruns/'
            outdir = fbadir
            tarf = maindir+'/multipass_mock'+str(mocknum)+'_'+pdir+'/inputs/targ.fits'

            indir = maindir+'/multipass_mock'+str(mocknum)+'_'+pdir+'/'
            asn = mocktools.combtiles_assign_wdup_7pass(indir,outdir,tarf,tp=pdir)
            asn['ZWARN_MTL'] = np.copy(asn['ZWARN'])
            pa = mocktools.combtiles_pa_wdup_7pass(indir,outdir,tarf,addcols=['TARGETID','RA','DEC'],fba=True,tp=pdir,ran='dat')
        
        elif args.famd == 'ab_secondgen':

            print('entering altmtl')
            tarf = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/forFA%d.fits' % mocknum #os.path.join(maindir, 'forFA_Real%d.fits' % mocknum)
            fbadir = os.path.join(maindir, 'altmtl_main_rea00%d' % mocknum, 'Univ000', 'fa', 'MAIN')
            print('entering common.combtiles_wdup_altmtl for FASSIGN')

            asn = common.combtiles_wdup_altmtl('FASSIGN', tiles, fbadir, os.path.join(outdir, 'datcomb_' + pdir + 'assignwdup.fits'), tarf, addcols=['TARGETID','RSDZ','TRUEZ','ZWARN','PRIORITY'])
            #if using alt MTL that should have ZWARN_MTL, put that in here
            asn['ZWARN_MTL'] = np.copy(asn['ZWARN'])
            print('entering common.combtiles_wdup_altmtl for FAVAIL')
            pa = common.combtiles_wdup_altmtl('FAVAIL', tiles, fbadir, os.path.join(outdir, 'datcomb_' + pdir + 'wdup.fits'), tarf, addcols=['TARGETID','RA','DEC','PRIORITY_INIT','DESI_TARGET'])
        else:
            tarf = fbadir+'/targs.fits'
            asn = common.combtiles_assign_wdup(tiles,fbadir,outdir,tarf,tp=pdir)
            #if using alt MTL that should have ZWARN_MTL, put that in here
            asn['ZWARN_MTL'] = np.copy(asn['ZWARN'])
            pa = common.combtiles_pa_wdup(tiles,fbadir,outdir,tarf,addcols=['TARGETID','RA','DEC'],fba=True,tp=pdir,ran='dat')
        
        pa['TILELOCID'] = 10000*pa['TILEID'] +pa['LOCATION']
        tj = join(pa,asn,keys=['TARGETID','LOCATION','TILEID'],join_type='left')
        outfs = os.path.join(lssdir, 'datcomb_'+pdir+'_tarspecwdup_zdone.fits')
        tj.write(outfs,format='fits', overwrite=True)
        print('wrote '+outfs)
        tc = ct.count_tiles_better('dat', pdir, specrel='', survey=args.survey, indir=lssdir, gtl=gtl) 
        outtc =  os.path.join(lssdir, 'Alltiles_'+pdir+'_tilelocs.dat.fits')
        tc.write(outtc,format='fits', overwrite=True)
        print('wrote '+outtc)
    
    if args.combdr == 'y':
        fbadir_data = os.path.join(maindir,'fba'+str(mocknum))
        fbadir_ran = maindir+'random_fba'+str(rannum)
        if args.famd == 'passes':
            fbadir_data = maindir+'/multipass_mock'+str(mocknum)+'_'+pdir+'/faruns/'
            fbadir_ran = maindir+'/ran'+str(rannum)+'_'+pdir+'/faruns/'

        specf = Table(fitsio.read(fbadir_data+'/datcomb_'+pdir+'assignwdup.fits'))
        specf.remove_columns(['TARGETID'])
        fgu = Table(fitsio.read(fbadir_ran+'/rancomb_'+pdir+'wdup.fits'))
        print(len(fgu))
        fgu = join(fgu,specf,keys=['LOCATION','TILEID'],join_type='left')
        fgu['TILELOCID'] = 10000*fgu['TILEID'] +fgu['LOCATION']
        del specf
        print(len(fgu))
        print(fgu.dtype.names)
        fgu.sort('TARGETID')
        outf = lssdir+'/rancomb_'+str(rannum)+pdir+'wdupspec_zdone.fits'
        
        fgu.write(outf,format='fits', overwrite=True)
        print('wrote '+outf)
        del fgu
     
    if args.countran == 'y':
        if mocknum == 0:
            tc = ct.count_tiles_better('ran',pdir,rannum,specrel='',survey=args.survey,indir=lssdir,gtl=gtl)
            print('got counts')
            print(len(tc))
            print(tc.dtype.names)
            print(np.max(tc['NTILE']))
            tc.keep_columns(['TARGETID','NTILE','TILES'])
            #common.write_LSS(tc,lssdir+'/rancomb_'+str(rannum)+pdir+'_Alltilelocinfo.fits')
            tc.write(lssdir+'/rancomb_'+str(rannum)+pdir+'_Alltilelocinfo.fits',format='fits', overwrite=True)
            print('wrote random counts')
        else:
            fn0 = maindir+'mock0/rancomb_'+str(rannum)+pdir+'_Alltilelocinfo.fits'
            fn = lssdir+'/rancomb_'+str(rannum)+pdir+'_Alltilelocinfo.fits'
            os.system('cp '+fn0+' '+fn)
            print('copied '+fn0+' to '+fn)

    #specver = 'mock'    
    imbits = []    
    if args.fulld == 'y':
        
        mainp = main(args.tracer, args.specdata, survey=args.survey)
#inp        imbits = mainp.imbits
        

        #Question, is this necessary? Should be this or the list below?
        maxp = 3400
        if type[:3] == 'LRG' or notqso == 'notqso':
            maxp = 3200
        if type[:3] == 'BGS':
            maxp = 2100

        ftar = None
        dz = os.path.join(lssdir, 'datcomb_'+pdir+'_tarspecwdup_zdone.fits')
        tlf = os.path.join(lssdir, 'Alltiles_'+pdir+'_tilelocs.dat.fits')

        fcoll = os.path.join(lssdir, 'collision_'+pdir+'_mock%d.fits' % mocknum)
        
        if not os.path.isfile(fcoll):
            fin = os.path.join('/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit','mock%d' %mocknum, 'pota-' + pr + '.fits')
            fcoll = mocktools.create_collision_from_pota(fin, fcoll)
        else:
            print('collision file already exist', fcoll)

        ct.mkfulldat(dz, imbits, ftar,args.tracer,bit,os.path.join(dirout, args.tracer+notqso+'_full_noveto.dat.fits'),tlf, survey=args.survey, maxp=maxp, desitarg=desitarg, specver=args.specdata, notqso=notqso, gtl_all=None, mockz=mockz,  mask_coll=fcoll, badfib=mainp.badfib, min_tsnr2=mainp.tsnrcut, mocknum=mocknum)

    maxp = 3400
    pthresh = 3000
    zmin = 0.8
    zmax = 3.5
    P0 = 6000
    dz_step = 0.02

    if args.tracer[:3] == 'LRG':# or notqso == 'notqso':
        maxp = 3200
        P0 = 10000
        dz_step = 0.01
        zmin = 0.4
        zmax = 1.1
    if args.tracer[:3] == 'ELG':
        P0 = 4000
        dz_step = 0.01
        maxp = 3000
        zmin = 0.8
        zmax = 1.6
    if args.tracer[:3] == 'BGS':
        P0 = 7000
        dz_step = 0.01
        maxp = 2100
        pthresh = 2000
        zmin = 0.1
        zmax = 0.5
    if notqso == 'notqso':
        maxp = 3200

    nzmd = 'mock'
        
    if args.fullr == 'y':
        if args.survey != 'Y1':
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
        else:
            #ldata = os.path.join('/global/cfs/cdirs/desi/survey/catalogs', args.survey,'LSS', args.specdata)
            ldata = os.path.join(maindir, 'mock%d'% mocknum, 'datcomb_' + pdir + '_tarspecwdup_zdone.fits')
            specft = fitsio.read(ldata) #Is this from data or mock? 
            #specft = fitsio.read(os.path.join(ldata, 'datcomb_'+args.tracer+notqso+'_tarspecwdup_zdone.fits')) #Is this from data or mock? 
            wg = np.isin(specft['TILELOCID'], gtl)
            specft = Table(specft[wg])
            lznp = common.find_znotposs(specft)
            del specft

            
            ranfile = os.path.join('/global/cfs/cdirs/desi/survey/catalogs', args.survey,'LSS', args.specdata, 'rancomb_%d%swdupspec_zdone.fits' % (rannum, pdir)) 
            alltileloc = os.path.join('/global/cfs/cdirs/desi/survey/catalogs', args.survey,'LSS', args.specdata, 'rancomb_%d%s_Alltilelocinfo.fits' % (rannum, pdir)) 
#            os.system('cp %s %s' %(ranfile, os.path.join(maindir, 'mock'+str(mocknum)+'/.')))
#            os.system('cp %s %s' %(alltileloc, os.path.join(maindir, 'mock'+str(mocknum)+'/.')))

#            ldata = os.path.join('/global/cfs/cdirs/desi/survey/catalogs', args.survey,'LSS', args.specdata)
#            specft = fitsio.read(os.path.join(ldata, 'datcomb_'+pdir+'_spec_zdone.fits'), columns=['TILELOCID','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])

#            mockspec = fitsio.read(os.path.join(maindir, 'fba'+str(mocknum), 'datcomb_' + pdir + 'assignwdup.fits'))
#            mockspec['TILELOCID'] = 10000*mockspec['TILEID'] + mockspec['LOCATION']
#            mockspec = join(mockspec, specft, keys=['TILELOCID'], join_type='left')
#            os.system('mv %s %s' %(os.path.join(maindir, 'fba'+str(mocknum), 'datcomb_' + pdir + 'assignwdup.fits'), os.path.join(maindir, 'fba'+str(mocknum), 'datcomb_' + pdir + 'assignwdup_original.fits'))
#            mockspec.write(os.path.join(maindir, 'fba'+str(mocknum), 'datcomb_' + pdir + 'assignwdup.fits'))
            ranfile, alltileloc = mocktools.createrancomb_wdupspec(lssdir, ranfile, alltileloc, os.path.join(maindir, 'fba'+str(mocknum), 'datcomb_' + pdir + 'assignwdup.fits'), os.path.join('/global/cfs/cdirs/desi/survey/catalogs', args.survey,'LSS', args.specdata, 'datcomb_'+pdir+'_spec_zdone.fits'))
#            print(ranfile, alltileloc)
#            exit()
        outf = os.path.join(dirout, args.tracer+notqso+'_'+str(rannum)+'_full_noveto.ran.fits')
        mainp = main(args.tracer, args.specdata, survey=args.survey)
        imbits = mainp.imbits
        maxp = 3400
        if args.tracer[:3] == 'LRG' or notqso == 'notqso':
            maxp = 3200
        if args.tracer[:3] == 'BGS':
            maxp = 2100
        tsnrcut = mainp.tsnrcut
        ct.mkfullran(gtl,lznp,os.path.join(maindir, 'mock'+str(mocknum)),rannum,imbits,outf,args.tracer,pdir,notqso=notqso,maxp=maxp,min_tsnr2=tsnrcut) 
##        ct.mkfullran(gtlf,lznp,lssdir,rannum,imbits,outf,args.tracer,pdir,notqso=notqso,maxp=maxp,tlid_full=tlid_full)
        

    mainp = main(args.tracer, args.specdata, survey=args.survey)

    if args.apply_veto == 'y':
        print('applying vetos to mock '+str(mocknum))
        fin = os.path.join(dirout, args.tracer+notqso+'_full_noveto.dat.fits')
        if args.tracer == 'LRG' and args.survey=='Y1' and args.famd == 'ab_secondgen':
            lrgf = Table.read(fin)
            lrgmask = Table.read('/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/forFA%d_matched_input_full_lrg_imask.fits' % mocknum)
            lrgf = join(lrgf, lrgmask, keys=['TARGETID'])
            common.write_LSS(lrgf, fin)

        elif args.survey=='Y1' and args.famd == 'ab_secondgen':
            novf = Table.read(fin)
            maskcols = []
            cols = ['NOBS_G','NOBS_R','NOBS_Z','MASKBITS']
            ovw = False
            for colm in cols:
                if colm not in novf.columns:
                    maskcols.append(colm)
            if len(maskcols) > 0:
                maskcols.append('TARGETID')
                targf = Table(fitsio.read('/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/forFA%d.fits' % mocknum, columns=maskcols))
                novf = join(novf, targf, keys=['TARGETID'])
                ovw = True
            if 'PHOTSYS' not in novf.columns:
                novf['PHOTSYS'] = 'N'
                sel = novf['DEC'] < 32.375
                wra = (novf['RA'] > 100-novf['DEC'])
                wra &= (novf['RA'] < 280 +novf['DEC'])
                sel |= ~wra
                novf['PHOTSYS'][sel] = 'S'
                ovw = True
            if ovw:
                common.write_LSS(novf, fin)

        fout = os.path.join(dirout, args.tracer+notqso+'_full.dat.fits')
        common.apply_veto(fin, fout, ebits=mainp.ebits, zmask=False, maxp=maxp, reccircmasks=mainp.reccircmasks)
        print('applying vetos to random '+str(rannum))
        fin = os.path.join(dirout, args.tracer+notqso+'_'+str(rannum)+'_full_noveto.ran.fits')
        fout = os.path.join(dirout, args.tracer+notqso+'_'+str(rannum)+'_full.ran.fits')
        if args.tracer!= 'QSO':
            common.add_veto_col(fin, ran=True, tracer_mask=args.tracer[:3].lower(), rann=rannum)
        common.apply_veto(fin, fout, ebits=mainp.ebits, zmask=False, maxp=maxp, reccircmasks=mainp.reccircmasks)

        #print('random veto '+str(ii)+' done')
    if args.apply_veto_ran == 'y':
        fin = os.path.join(dirout, args.tracer+notqso+'_'+str(rannum)+'_full_noveto.ran.fits')
        fout = os.path.join(dirout, args.tracer+notqso+'_'+str(rannum)+'_full.ran.fits')
        if args.tracer!= 'QSO':
            common.add_veto_col(fin, ran=True, tracer_mask=args.tracer[:3].lower(),rann=rannum)
        common.apply_veto(fin,fout,ebits=mainp.ebits, zmask=False,maxp=maxp, reccircmasks=mainp.reccircmasks)


    regl = ['_N','_S']    


    wzm = ''

    lssmapdirout = os.path.join(dirout,'hpmaps')
    nside = 256
    tracer_clus = args.tracer+notqso+wzm

    if args.mkHPmaps == 'y':
        from LSS.imaging.sky_maps import create_pixweight_file, rancat_names_to_pixweight_name
        if not os.path.exists(lssmapdirout):
            os.mkdir(lssmapdirout)
            print('made '+lssmapdirout)
        new_cols=mainp.new_cols#['STARDENS','HALPHA', 'HALPHA_ERROR', 'CALIB_G', 'CALIB_R', 'CALIB_Z', 'EBV_MPF_Mean_FW15', 'EBV_MPF_Mean_ZptCorr_FW15', 'EBV_MPF_Var_FW15', 'EBV_MPF_VarCorr_FW15', 'EBV_MPF_Mean_FW6P1', 'EBV_MPF_Mean_ZptCorr_FW6P1', 'EBV_MPF_Var_FW6P1', 'EBV_MPF_VarCorr_FW6P1', 'EBV_SGF14', 'BETA_ML', 'BETA_MEAN', 'BETA_RMS', 'HI', 'KAPPA_PLANCK']
        fid_cols=mainp.fid_cols#['EBV','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
        

        


        fieldslist = new_cols+fid_cols
        lssmapdir = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/external_input_maps/'
        rancatname = os.path.join(dirout, tracer_clus+'_*_full.ran.fits')
        rancatlist = sorted(glob.glob(rancatname))
        for temp_rannum in range(0, 18): #rancat in rancatlist:

            rancat = os.path.join(dirout, tracer_clus+'_%d_full.ran.fits'%temp_rannum) 
            
            print('filling randoms with imaging properties for',rancat, temp_rannum)
            common.add_map_cols(rancat, temp_rannum, new_cols=new_cols, fid_cols=fid_cols)
            print('done with ',str(temp_rannum))

        masklist = list(np.zeros(len(fieldslist),dtype=int))
        print('aure', fieldslist) 
        for reg in ['N','S']:
            outfn = os.path.join(lssmapdirout,tracer_clus+'_mapprops_healpix_nested_nside'+str(nside)+'_'+reg+'.fits')
            create_pixweight_file(rancatlist, fieldslist, masklist, nside_out=nside,
                          lssmapdir=lssmapdir, outfn=outfn,reg=reg)    

    if args.apply_map_veto == 'y':
        import healpy as hp
        lssmapdirout = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v0.6/hpmaps'
        mapn = fitsio.read(os.path.join(lssmapdirout, tracer_clus+'_mapprops_healpix_nested_nside'+str(nside)+'_N.fits'))
        maps = fitsio.read(os.path.join(lssmapdirout, tracer_clus+'_mapprops_healpix_nested_nside'+str(nside)+'_S.fits'))
        mapcuts = mainp.mapcuts

#        if args.ranonly != 'y':
        fout = os.path.join(dirout, args.tracer+notqso+'_full.dat.fits')
        fin = fout.replace('global','dvs_ro')
        fout = fout.replace('_full','_full_HPmapcut')
        common.apply_map_veto(fin, fout, mapn, maps, mapcuts)
        print('data veto done, now doing randoms')
        
        def _parfun(rn):
            fout = os.path.join(dirout, args.tracer+notqso+'_'+str(rn)+'_full.ran.fits')
            fin = fout.replace('global', 'dvs_ro')
            fout = fout.replace('_full', '_full_HPmapcut')
            common.apply_map_veto(fin, fout, mapn, maps, mapcuts)
            print('random veto '+str(rn)+' done')
        if args.par == 'n':
            for rn in range(rm,rx):
                _parfun(rn)
        else:
            inds = np.arange(rm,rx)
            from multiprocessing import Pool
            (rx-rm)*2
            nproc = 9 #try this so doesn't run out of memory
            with Pool(processes=nproc) as pool:
                res = pool.map(_parfun, inds)

    
##FKP OVER FULL SAMPLE
    ##fin = fitsio.read(fcd_in)
##        cols = list(fin.dtype.names)
##        nz_in = common.mknz_full(fcd_in, fcr_in, type[:3], bs=dz, zmin=zmin, zmax=zmax, write=wo, randens=randens, md=nzmd)

##        if 'WEIGHT_FKP' not in cols:
##            print('adding FKP weights')
##            common.addFKPfull(fcd_in, nz_in, type[:3], bs=dz, zmin=zmin, zmax=zmax, P0=P0, md=nzmd)

    if args.getFKP == 'y':
        randens = 2500.
        nzmd = 'mock'
        fcd_in = os.path.join(dirout, args.tracer+notqso + '_full'+args.use_map_veto+'.dat.fits')
        fcr_in = os.path.join(dirout, args.tracer+notqso + '_0_full'+args.use_map_veto+'.ran.fits')
        nzf_in = os.path.join(dirout, args.tracer+notqso + '_full_nz.txt')
        wo = 'y'
        if os.path.isfile(nzf_in):
            wo = 'n'
        fin = fitsio.read(fcd_in)
        cols = list(fin.dtype.names)
        nz_in = common.mknz_full(fcd_in, fcr_in, args.tracer[:3], bs=dz_step, zmin=zmin, zmax=zmax, write=wo, randens=randens, md=nzmd)

        if 'WEIGHT_FKP' not in cols:
            print('adding FKP weights')
            common.addFKPfull(fcd_in, nz_in, type[:3], bs=dz_step, zmin=zmin, zmax=zmax, P0=P0, md=nzmd)

    nztl = []
    if args.mkclusdat == 'y':
        nztl.append('')
        fin = os.path.join(dirout, args.tracer+notqso + '_full'+args.use_map_veto+'.dat.fits')
        dz = Table.read(fin)
        if 'PHOTSYS' not in dz.columns:
            dz['PHOTSYS'] = 'N'
            sel = dz['DEC'] < 32.375
            wra = (dz['RA'] > 100-dz['DEC'])
            wra &= (dz['RA'] < 280 +dz['DEC'])
            sel |= ~wra
            dz['PHOTSYS'][sel] = 'S'
            common.write_LSS(dz, fin)

        #ct.mkclusdat(os.path.join(dirout,args.tracer+notqso),tp=args.tracer,dchi2=None,tsnrcut=0,zmin=zmin,zmax=zmax)#,ntilecut=ntile)
        ct.mkclusdat(os.path.join(dirout,args.tracer+notqso),tp=args.tracer,dchi2=None,tsnrcut=0,zmin=zmin,zmax=zmax, use_map_veto=args.use_map_veto)#,ntilecut=ntile,ccut=ccut)


    
    if args.equal_data_dens == 'y':
        data_dir = "/dvs_ro/project/projectdirs/desi/survey/catalogs/edav1/da02/LSScats/clustering"
        
        for i in ['N','S']:
            mocktools.mock_equal_data_density(dirout, data_dir, dirout, args.tracer, i, zmin, zmax, args.nran_clus_data, randens)
        
        

    if args.mkclusran == 'y':
        if len(nztl) == 0:
            nztl.append('')
        rcols=['Z','WEIGHT','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL']
        if args.getFKP == 'y':
            rcols.append('WEIGHT_FKP')
        tsnrcol = 'TSNR2_ELG'
        if args.tracer[:3] == 'BGS':
            tsnrcol = 'TSNR2_BGS'
        fl = os.path.join(dirout, args.tracer+notqso+'_')
        ct.add_tlobs_ran(fl,rannum)
        ct.mkclusran(os.path.join(dirout,args.tracer+notqso+'_'),os.path.join(dirout,args.tracer+notqso+'_'), rannum, rcols=rcols, tsnrcut=0, tsnrcol=tsnrcol, use_map_veto=args.use_map_veto)#,ntilecut=ntile,ccut=ccut)
        #for clustering, make rannum start from 0
        '''
        for reg in regl:
            ranf = dirout+args.tracer+notqso+reg+'_'+str(rannum)+'_clustering.ran.fits'
            ranfm = dirout+args.tracer+notqso+reg+'_'+str(rannum-1)+'_clustering.ran.fits'
            os.system('mv '+ranf+' '+ranfm)
        '''

        ct.clusNStoGC(os.path.join(dirout, args.tracer+notqso+'_'), args.maxr - args.minr)
    
    nran = args.maxr-args.minr
    regions = ['NGC', 'SGC']

    if args.resamp == 'y':
        for reg in regions:
            flin = os.path.join(dirout, tracer_clus + '_'+reg)
        def _parfun(rannum):
            ct.clusran_resamp(flin, rannum, rcols=rcols)#,compmd=args.compmd)#, ntilecut=ntile, ccut=ccut)

        inds = np.arange(nran)
        if args.par == 'y':
            from multiprocessing import Pool
            with Pool(processes=nran*2) as pool:
                res = pool.map(_parfun, inds)
        else:
            for rn in range(rm,rx):
                _parfun(rn)

    allreg = ['N','S','NGC', 'SGC']
    if args.nz == 'y':
        for reg in allreg:
            fb = os.path.join(dirout, tracer_clus+'_'+reg)
            fcr = fb + '_0_clustering.ran.fits'
            fcd = fb + '_clustering.dat.fits'
            fout = fb + '_nz.txt'
            common.mknz(fcd,fcr,fout,bs=dz_step,zmin=zmin,zmax=zmax)
            common.addnbar(fb, bs=dz_step, zmin=zmin, zmax=zmax, P0=P0, nran=nran)

    '''
    if args.FKPfull == 'y':
    
        fb = os.path.join(dirout, tracer_clus)
        fbr = fb
        if type == 'BGS_BRIGHT-21.5':
            fbr = dirout+'BGS_BRIGHT'

        fcr = fbr +' _0_full'+args.use_map_veto+'.ran.fits'
        fcd = fb + '_full'+args.use_map_veto+'.dat.fits'
        nz = common.mknz_full(fcd, fcr, args.tracer[:3], bs=dz_step, zmin=zmin, zmax=zmax)
        common.addFKPfull(fcd, nz, args.tracer[:3], bs=dz, zmin=zmin, zmax=zmax, P0=P0)
    '''

    if args.mkclusdat_allpot == 'y':
        #fbadir = maindir+'fba'+str(mocknum)
        #tarf = fbadir+'/targs.fits'
        #'/global/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit
        tarf = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/forFA%d.fits' % mocknum
        #args.base_output +mockdir+'/forFA'+str(mocknum)+'.fits'
        ztab = Table(fitsio.read(tarf,columns=['TARGETID','RSDZ']))
        ztab.rename_column('RSDZ', 'Z')
        
        mocktools.mkclusdat_allpot(os.path.join(dirout, args.tracer+notqso), ztab, tp=args.tracer, dchi2=None, tsnrcut=0, zmin=zmin, zmax=zmax)#,ntilecut=ntile)



    if args.mkclusran_allpot == 'y':
        nztl.append('_complete')
        rcols=['Z','WEIGHT']
        tsnrcol = 'TSNR2_ELG'
        if args.tracer[:3] == 'BGS':
            tsnrcol = 'TSNR2_BGS'
        fl = os.path.join(dirout, args.tracer+notqso+'_')
        ct.add_tlobs_ran(fl,rannum)
        ct.mkclusran(os.path.join(dirout, args.tracer+notqso+'_'), os.path.join(dirout,args.tracer+notqso+'_complete'+'_'), rannum, rcols=rcols, tsnrcut=0, tsnrcol=tsnrcol, compmd='')#,ntilecut=ntile,ccut=ccut)
        #for clustering, make rannum start from 0
##        for reg in regl:
##            ranf = dirout+args.tracer+notqso+'_complete'+reg+'_'+str(rannum)+'_clustering.ran.fits'
##            ranfm = dirout+args.tracer+notqso+'_complete'+reg+'_'+str(rannum-1)+'_clustering.ran.fits'
##            os.system('mv '+ranf+' '+ranfm)
    #THIS IS TEMP

        ct.clusNStoGC(os.path.join(dirout, args.tracer+notqso+'_complete_'), nran=1)
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

    '''
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
    '''
    

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
#     if args.par == 'y':
#         from multiprocessing import Pool
#         from desitarget.internal import sharedmem
#         
#         N = rx-rm+1
#         inds = []
#         for i in range(rm,rx):
#             inds.append(i)
#         pool = sharedmem.MapReduce(np=N)
#         with pool:
#         
#             def reduce( r):
#                 print('chunk done')
#                 return r
#             pool.map(prep,inds,reduce=reduce)
# 
#         #p.map(doran,inds)
#     else:
    for i in range(rm,rx):
        if args.par == 'y':
            from multiprocessing import Pool
            from itertools import repeat
            from desitarget.internal import sharedmem
            N = mockmax-mockmin+1       
            inds = []
            for mn in range(mockmin,mockmax):
                 inds.append((mn,i))
            #pool = sharedmem.MapReduce(np=N)
            #with pool:
            with Pool(processes=N) as pool:
                #def reduce( r):
                #    print('mock done')
                #    return r
                pool.starmap(docat,inds)#,reduce=reduce)
        for mn in range(mockmin,mockmax):
        
            if args.par != 'y':
                print('processing mock '+str(mn)+' and random '+str(i))
                docat(mn,i)
            if args.split_GC == 'y':
                nztl = ['','_complete']
                lssdir = maindir+'mock'+str(mn)+'/'

                dirout = lssdir+'LSScats/'
            
                for zt in nztl:
                    fb = dirout+args.tracer+notqso+zt+'_'                
                    ct.clusNStoGC(fb,rx-rm)
