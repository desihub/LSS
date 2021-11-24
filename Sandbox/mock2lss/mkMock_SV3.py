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
from desimodel.footprint import is_point_in_desi
from desitarget.sv3 import sv3_targetmask

#from this package
import LSS.mkCat_singletile.fa4lsscat as fa
from LSS.globals import SV3 

import mockcattools as myct

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='everest')
parser.add_argument("--cuttar", help="cut targets to SV3 tiles",default='n')
parser.add_argument("--vis", help="make a plot of data/randoms on tile",default='n')
parser.add_argument("--xi", help="run pair-counting code",default='n')
parser.add_argument("--mockmtl", help="make a mock mtl file for the tile",default='n')
parser.add_argument("--combd", help="combine all the tiles together",default='y')
parser.add_argument("--fulld", help="make the 'full' catalog containing info on everything physically reachable by a fiber",default='y')
parser.add_argument("--clus", help="make the data clustering files; these are cut to a small subset of columns",default='y')
parser.add_argument("--maskz", help="apply sky line mask to redshifts?",default='n')

parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')

parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')

#all random set to n by default since mkCat_SV3_ran.py exists and does it in parallel

args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
version = args.version
specrel = args.verspec

SV3p = SV3(type)

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'


zma = False
if args.maskz == 'y':
    zma = True
ctar = False
if args.cuttar == 'y':
    ctar = True
docatplots = False
if args.vis == 'y':
    docatplots = True
doclus = False
if args.xi == 'y':    
    doclus = True
mkfulld = True #make the 'full' catalog containing info on everything physically reachable by a fiber
if args.fulld == 'n':
    mkfulld = False
mkclus = True #make the data/random clustering files; these are cut to a small subset of columns
mkclusdat = True
if args.clus == 'n':
    mkclus = False
    mkclusdat = False
mknz = False #get n(z) for type and all subtypes
if args.nz == 'y':
    mknz = True

fillNZ = False

combd = True
if args.combd == 'n':
    combd = False
mkmockmtl = False
if args.mockmtl == 'y':
    mkmockmtl = True



if type == 'dark' or type == 'bright':
    #set all type catalog stuff to False in this case
    mkfulld = False
    mkfullr = False
    mkclus = False
    mkclusdat = False
    mkclusran = False
    
def test_dir(value):
    if not os.path.exists(value):
        try:
            os.makedirs(value, 0o755)
            print('made %s'%value)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise    

if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    pr = 'BRIGHT'
    pdir = 'bright'
else:
    pr = 'DARK'
    pdir = 'dark'

pd = pdir

mdir = os.path.join(SV3p.mdir, pdir) #location of ledgers
tdir = os.path.join(SV3p.tdir, pdir) #location of targets
mtld = SV3p.mtld
tiles = SV3p.tiles
imbits = SV3p.imbits #mask bits applied to targeting
ebits = SV3p.ebits #extra mask bits we think should be applied

#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
sv3dir = os.path.join(basedir,'SV3', 'LSS')

#tarbit = int(np.log2(sv3_targetmask.desi_mask[type]))

wp = tiles['PROGRAM'] == pr
tiles = tiles[wp]
print(len(tiles))

wp = np.isin(mtld['TILEID'],tiles['TILEID']) #we want to consider MTL done tiles that correspond to the SV3 tile file
mtld = mtld[wp]
print(len(mtld))


test_dir(os.path.join(sv3dir,'logs'))

ldirspec = os.path.join(sv3dir, specrel)
test_dir(ldirspec)

test_dir(os.path.join(ldirspec,'LSScats'))

dirout = os.path.join(ldirspec,'LSScats', version)
test_dir(dirout)


tiles_ta_file = os.path.join(sv3dir, 'tiles-'+pr+'.fits')
#construct a table with the needed tile information
do_tile_info_again = False

if os.path.isfile(tiles_ta_file) and not do_tile_info_again:
    print('tile info already exist, reading ta from file')
    ta = Table.read(tiles_ta_file)
else:
    print('tile info do not exist, or you ask to run again')
    if len(mtld) > 0:
        tilel = []
        ral = []
        decl = []
        mtlt = []
        fal = []
        obsl = []
        pl = []
        hal = []
        #for tile,pro in zip(mtld['TILEID'],mtld['PROGRAM']):
        for tile in mtld['TILEID']:
            ts = str(tile).zfill(6)
            fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
            tilel.append(tile)
            ral.append(fht['TILERA'])
            decl.append(fht['TILEDEC'])
            mtlt.append(fht['MTLTIME'])
            fal.append(fht['RUNDATE'])
            obsl.append(fht['FIELDROT'])
            hal.append(fht['FA_HA'])
            #pl.append(pro)
            pl.append(pr)
        ta = Table()
        ta['TILEID'] = tilel
        ta['RA'] = ral
        ta['DEC'] = decl
        ta['MTLTIME'] = mtlt
        ta['RUNDATE'] = fal
        ta['FIELDROT'] = obsl
        ta['PROGRAM'] = pl
        ta['FA_HA'] = hal
        #if pd == 'dark':
        ta['OBSCONDITIONS'] = 15
        ta['IN_DESI'] = 1
        ta.write(tiles_ta_file, format='fits', overwrite=True)

    else:
        print('no done tiles in the MTL')

minr = 148
maxr = 274
mind = -2.5
maxd = 68

my_path = '/global/cscratch1/sd/acarnero/fiberassign'
target_file = os.path.join(my_path, 'targets-UNIT-mtlz_SV3_alltracers_sv3bits_v2.fits')

cutsv3_target_file = os.path.join(sv3dir, 'alltilesnofa.fits')

if ctar and not cutsv3_target_file:
    ffile, h = fitsio.read(target_file, header=True)

    print('targets before anything', len(ffile))
    wp = ffile['RA'] > minr
    wp &= ffile['RA'] < maxr
    wp &= ffile['DEC'] > mind
    wp &= ffile['DEC'] < maxd
    ffile = ffile[wp]
    print('targets after selecting by min, max in RA, Dec', len(ffile))
    tilesall = tiles #Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')
    tilesu = unique(tilesall, keys=['RA','DEC'])                
    wi = is_point_in_desi(tilesu, ffile["RA"], ffile["DEC"])
    ffile = ffile[wi]
    print('targets after is point desi', len(ffile))
    fitsio.write(cutsv3_target_file, ffile, clobber=True, header=h)
    print('wrote '+cutsv3_target_file+' from '+target_file)
    target_file = cutsv3_target_file
else:
    print('targets selected in tile already')
    target_file = cutsv3_target_file


if mkmockmtl:
    test_dir('./atest')
    myct.randomtiles_allSV3(ta, target_file, directory_output='./atest')

debug=False
if combd:
    if type == 'dark' or type == 'bright':

        outf = os.path.join(sv3dir,'datcomb_'+type+'_tarwdup_Alltiles.fits')
        if not debug:
            myct.combtiles_wdup(ta, ['./atest', 'tilenofa-{TILE}.fits'], ['./fiberassigment', 'mocks000{TILE}.fits'] , fout=outf)

        print('yeah!')
        tarf = Table.read(outf)
        tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
        remcol = ['PRIORITY','Z','ZWARN','FIBER','SUBPRIORITY'] #subpriority in target files doesn't match what is in fiberassign files
        for col in remcol:
            try:
                tarf.remove_columns([col] )#we get this where relevant from spec file
            except:
                print('column '+col +' was not in tarwdup file')    

        if specrel == 'everest':
            specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-sv3-'+type+'-cumulative.fits')
            wt = np.isin(specf['TILEID'],ta['TILEID']) #cut spec file to dark or bright time tiles
            specf = specf[wt]
            specf.keep_columns(['TARGETID','CHI2','COEFF','Z','ZERR','ZWARN','NPIXELS','SPECTYPE','SUBTYPE','NCOEFF','DELTACHI2'\
            ,'LOCATION','FIBER','COADD_FIBERSTATUS','TILEID','FIBERASSIGN_X','FIBERASSIGN_Y','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT'\
            ,'MEAN_DELTA_X','MEAN_DELTA_Y','RMS_DELTA_X','RMS_DELTA_Y','MEAN_PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B'\
            ,'TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
            'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
            'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
            tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID'],join_type='left')
            specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
        '''            
        elif specrel == 'daily':
            outf = ldirspec+'datcomb_'+type+'_specwdup_Alltiles.fits'
            ct.combtile_spec(mtld,outf,rel=specrel)
            specf = Table.read(ldirspec+'datcomb_'+type+'_specwdup_Alltiles.fits')
            specf.keep_columns(['CHI2','COEFF','Z','ZERR','ZWARN','NPIXELS','SPECTYPE','SUBTYPE','NCOEFF','DELTACHI2'\
            ,'FIBERASSIGN_X','FIBERASSIGN_Y','LOCATION','FIBER','FIBERSTATUS','PRIORITY'\
            ,'DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE','NIGHT','EXPID','MJD','TILEID','INTEG_COADD_FLUX_B',\
            'MEDIAN_COADD_FLUX_B','MEDIAN_COADD_SNR_B','INTEG_COADD_FLUX_R','MEDIAN_COADD_FLUX_R','MEDIAN_COADD_SNR_R','INTEG_COADD_FLUX_Z',\
            'MEDIAN_COADD_FLUX_Z','MEDIAN_COADD_SNR_Z','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
            'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
            'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
            specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
            tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')
        '''
        elif specrel == 'mock':
            outfile_spec = os.path.join(ldirspec, 'datcomb_'+type+'_specwdup_Alltiles.fits')
            if not debug:
                myct.combtile_specmock(ta, ['./fiberassigment', 'mocks000{TILE}.fits'], target_file, outfile_spec)
            specf = Table.read(outfile_spec)
            specf.keep_columns(['FIBER','TARGETID','LOCATION','FIBERSTATUS','LAMBDA_REF','PETAL_LOC','DEVICE_LOC','DEVICE_TYPE','TARGET_RA','TARGET_DEC','FA_TARGET','FA_TYPE','FIBERASSIGN_X','FIBERASSIGN_Y','PLATE_RA','PLATE_DEC','TILEID','PRIORITY','SUBPRIORITY','ZWARN','TRUEZ'])
            specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
            print('targets', tarf.columns)
            print('spec', specf.columns)
            tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')

        try:
            print(np.unique(tj['SV3_DESI_TARGET'], return_counts=True))
        except:
            ftar = Table.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+type+'_targets.fits')
            ftar.keep_columns(['TARGETID','SV3_DESI_TARGET','SV3_BGS_TARGET','SV3_MWS_TARGET'])
            print(len(tj))
            tj = join(tj,ftar,keys=['TARGETID'])  
            print(len(tj))  
        tj.write(os.path.join(ldirspec, 'datcomb_'+type+'_tarspecwdup_Alltiles.fits'), format='fits', overwrite=True)

        tc = myct.count_tiles_better(specf, os.path.join(ldirspec, 'datcomb_'+type+'_tarspecwdup_Alltiles.fits'))
        tc.write(os.path.join(ldirspec, 'Alltiles_'+pdir+'_tilelocs.dat.fits'), format='fits', overwrite=True)
        print('ole')
    else:
        print('nothing to be done for combd, only done for dark/bright now')
        
        
if mkfulld:
    if specrel == 'everest':
        specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-sv3-'+pdir+'-cumulative.fits')
        wt = np.isin(specf['TILEID'],ta['TILEID']) #cut spec file to dark or bright time tiles
        specf = specf[wt]
    elif specrel == 'daily':
        specf = Table.read(ldirspec+'datcomb_'+type+'_specwdup_Alltiles.fits')
    elif specrel == 'mock':
        outfile_spec = os.path.join(ldirspec, 'datcomb_dark_specwdup_Alltiles.fits')
        specf = Table.read(outfile_spec)


    azf=''
    if type[:3] == 'ELG':
        azf = SV3p.elgzf#'/global/cfs/cdirs/desi/users/raichoor/everest/sv3-elg-everest-tiles.fits'
    if type[:3] == 'QSO':
        azf = SV3p.qsozf#'/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/QSO/QSO_catalog_SV3.fits'
    #'/global/homes/r/raichoor/sv3/sv3-elg-daily-thru20210521.fits'
    #/global/homes/r/raichoor/sv3/sv3-elg-daily-thru20210506.fits
    #dz = dirout+'datcomb_'+type+'_Alltiles.fits' old
    dz = os.path.join(ldirspec, 'datcomb_dark_tarspecwdup_Alltiles.fits') #new
    print(dz)
    if type == 'BGS_BRIGHT':
        bit = sv3_targetmask.bgs_mask[type]
        desitarg='SV3_BGS_TARGET'
    else:
        bit = sv3_targetmask.desi_mask[type]
        desitarg='SV3_DESI_TARGET'
    print(desitarg,pdir,bit)
    bitweightfile = None
    if pdir == 'dark':
        bitweightfile = SV3p.darkbitweightfile
        #bitweightfile='/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run64_2/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
        #bitweightfile='/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run64_2/BitweightsRound2/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
    elif pdir == 'bright':
        bitweightfile = SV3p.brightbitweightfile


    myct.mkfulldat(specf,dz,imbits,type,bit,os.path.join(dirout,type+notqso+'_full_noveto.dat.fits'),os.path.join(ldirspec, 'Alltiles_dark_tilelocs.dat.fits'), azf=azf, desitarg=desitarg,specver=specrel,notqso=notqso)
    #get_tilelocweight()
    #logf.write('ran get_tilelocweight\n')
    #print('ran get_tilelocweight\n')


#needs to happen before randoms so randoms can get z and weights
if mkclusdat:
    dchi2 = 9
    tsnrcut = 0
    if type[:3] == 'ELG':
        dchi2 = 0.9 #This is actually the OII cut criteria for ELGs
        tsnrcut = 80
    if type == 'LRG':
        dchi2 = 16  
        tsnrcut = 80  
    if type[:3] == 'BGS':
        dchi2 = 40
        tsnrcut = 800
    myct.mkclusdat(os.path.join(dirout,type+notqso+'_full_noveto.dat.fits'), zmask=zma, tp=type, dchi2=dchi2, tsnrcut=tsnrcut, ebits=None)
    #logf.write('ran mkclusdat\n')
    #print('ran mkclusdat\n')

    
if mknz:
    wzm = ''
#     if zmask:
#         wzm = 'zmask_'
    if rcut is not None:
        wzm += '_rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntile > 0:
        wzm += '_ntileg'+str(ntilecut)+'_'    
    if args.ccut is not None:
        wzm += '_'+args.ccut #you could change this to however you want the file names to turn out

    regl = ['','_N','_S']
    
    for reg in regl:
        fb = dirout+type+wzm+reg
        fcr = fb+'_0_clustering.ran.fits'
        fcd = fb+'_clustering.dat.fits'
        fout = fb+'_nz.dat'
        if type == 'QSO':
            zmin = 0.6
            zmax = 4.5
            dz = 0.05
            
        else:    
            dz = 0.02
            zmin = 0.01
            zmax = 1.61
        myct.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
        myct.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax)
        
