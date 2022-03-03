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
import random
#from this package
from LSS.globals import SV3 

import mockcattools as myct

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='everest')
parser.add_argument("--cutran", help="cut randoms to SV3 tiles",default='n')
parser.add_argument("--vis", help="make a plot of data/randoms on tile",default='n')
parser.add_argument("--xi", help="run pair-counting code",default='n')
parser.add_argument("--ranmtl", help="make a random mtl file for the tile",default='n')
parser.add_argument("--combr", help="combine the random tiles together",default='n')
parser.add_argument("--fullr", help="make the random files associated with the full data files",default='n')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--maskz", help="apply sky line mask to redshifts?",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=1) 
parser.add_argument("--id", help="Mock id",default=0)

parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')

parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')

#all random set to n by default since mkCat_SV3_ran.py exists and does it in parallel

args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
version = args.version
specrel = args.verspec

id_ = "%03d"%int(args.id)

SV3p = SV3(type)

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'


zma = False
if args.maskz == 'y':
    zma = True
cran = False
if args.cutran == 'y':
    cran = True
docatplots = False
if args.vis == 'y':
    docatplots = True
mkranmtl = False
if args.ranmtl == 'y':
    mkranmtl = True
mkfullr = True #make the random files associated with the full data files
if args.fullr == 'n':
    mkfullr = False
mkclusran = False
if args.clusran == 'y':
    mkclusran = True
mknz = False #get n(z) for type and all subtypes
if args.nz == 'y':
    mknz = True

fillNZ = False

combr = True
if args.combr == 'n':
    combr = False   


if type == 'dark' or type == 'bright':
    #set all type catalog stuff to False in this case
    mkfullr = False
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

randir = os.path.join(sv3dir, 'random')
test_dir(randir)

#randir = sv3dir+'random'
rm = int(args.minr)
rx = int(args.maxr)
print(rm,rx)
for i in range(rm,rx):
    test_dir(os.path.join(sv3dir, 'random'+str(i)))
    print('made '+str(i)+' random directory')

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

my_path = os.path.join(basedir,'SV3') #'/global/cscratch1/sd/acarnero/codes/LSS/Sandbox/mock2lss'

target_file = os.path.join(my_path, 'mockTargets_{ID}_FirstGen_CutSky_alltracers_sv3bits.fits'.format(ID=id_))

list_randoms = np.linspace(0, 9, num=10, dtype=np.int)
ranid = random.choice(list_randoms)

random_file = os.path.join(my_path, 'mockRandom_5X_{RANID}_FirstGen_CutSky_alltracers_sv3bits.fits'.format(RANID=ranid))

cutsv3_target_file = os.path.join(sv3dir, 'alltilesnofa_{ID}.fits'.format(ID=id_))
cutsv3_random_file = os.path.join(randir, 'alltilesnofa_random_{ID}.fits'.format(ID=id_))

file_info = open(os.path.join(basedir,'SV3','info_mock_{ID}_randomchoice'.format(ID=id_)), 'w')
file_info.write('random used '+ 'mockRandom_5X_{RANID}_FirstGen_CutSky_alltracers_sv3bits.fits'.format(RANID=ranid))
file_info.close()

if cran:# and not os.path.isfile(cutsv3_random_file):
#    dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
#Make for one random, then change for several ones    for ii in range(rm,rx):
##        ranf = fitsio.read(dirrt+'/randoms-1-'+str(ii)+'.fits')
        ranf, h = fitsio.read(random_file, header=True)
        print(len(ranf))
        #if ctar:
        wp = ranf['RA'] > minr
        wp &= ranf['RA'] < maxr
        wp &= ranf['DEC'] > mind
        wp &= ranf['DEC'] < maxd
        ranf = ranf[wp]
        print(len(ranf))
        tilesall = tiles #Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')
        tilesu = unique(tilesall,keys=['RA','DEC'])                
        wi = is_point_in_desi(tilesu, ranf["RA"], ranf["DEC"])
        ranf = ranf[wi]
        fitsio.write(cutsv3_random_file, ranf, clobber=True, header=h)
        print('wrote ' + cutsv3_random_file +' from ' + random_file)
        random_file = cutsv3_random_file

else:
        print('already there!')
        random_file = cutsv3_random_file

if mkranmtl:
    test_dir(os.path.join(basedir,'SV3','mtlran{ID}'.format(ID=id_)))
    myct.randomtiles_allSV3_parallel(ta, random_file, directory_output=os.path.join(basedir,'SV3','mtlran{ID}'.format(ID=id_)))

    

print('end here')
if combr:
    if type == 'dark' or type == 'bright':
        spec_file = os.path.join(ldirspec, 'datcomb_'+type+'_specwdup_Alltiles_{ID}.fits'.format(ID=id_))
        kc = ['ZWARN','FIBER','LOCATION','TILEID','TILELOCID','FIBERSTATUS','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY']#,'DELTA_X','DELTA_Y','EXPTIME','PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B','TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z','TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG']
        out_rancomb = os.path.join(ldirspec, 'rancomb_' + type + 'wdupspec_Alltiles_{ID}.fits'.format(ID=id_))
        myct.combran_wdup(ta,randir,type,ldirspec, Table.read(spec_file), [os.path.join(basedir,'SV3','fba_randoms'), 'ran_5X_0{RANID}_'.format(RANID=ranid)], [os.path.join(basedir,'SV3','mtlran{ID}'.format(ID=id_)), 'tilenofa-'], keepcols=kc, outf=[os.path.join(randir, 'rancomb_'+type+'wdup_Alltiles_{ID}.fits'.format(ID=id_)), out_rancomb])
        tc = myct.count_tiles_better(Table.read(spec_file), out_rancomb)
        tc.write(os.path.join(randir, 'rancomb_'+pdir+'_Alltilelocinfo_{ID}.fits'.format(ID=id_)), format='fits', overwrite=True)
    else:
        print('nothing to be done for combr, only done for dark/bright now')
        
        
if mkfullr:
    if specrel == 'everest':
        specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-sv3-'+pdir+'-cumulative.fits')
        fbcol = 'COADD_FIBERSTATUS'
    if specrel == 'daily':
        specf = Table.read(ldirspec+'datcomb_'+pdir+'_specwdup_Alltiles.fits')
        fbcol = 'FIBERSTATUS'
    if specrel == 'mock':
        specf = Table.read(os.path.join(ldirspec, 'datcomb_dark_specwdup_Alltiles_{ID}.fits'.format(ID=id_)))
        fbcol = 'FIBERSTATUS'

    outf = os.path.join(dirout, type+notqso+'_full_noveto_{ID}.ran.fits'.format(ID=id_))

    if type == 'BGS_BRIGHT':
        bit = sv3_targetmask.bgs_mask[type]
        desitarg='SV3_BGS_TARGET'
    else:
        bit = sv3_targetmask.desi_mask[type]    
        desitarg='SV3_DESI_TARGET'

    myct.mkfullran(specf, ldirspec, randir, imbits,outf,type,pdir,bit,desitarg=desitarg,tsnr='LOCATION_ASSIGNED',fbcol=fbcol,notqso=notqso, id_=id_)
#    for ii in range(rm,rx):
#        outf = dirout+type+'_'+str(ii)+'_full_noveto.ran.fits'
#        ct.mkfullran(randir,ii,imbits,outf,type,pdir,sv3_targetmask.desi_mask[type])
    #logf.write('ran mkfullran\n')
    #print('ran mkfullran\n')


if mkclusran:
    tsnrcol = 'TSNR2_ELG'
    tsnrcut = 0
    if type[:3] == 'ELG':
        #dchi2 = 0.9 #This is actually the OII cut criteria for ELGs
        tsnrcut = 80
    if type == 'LRG':
        #dchi2 = 16
        tsnrcut = 80
    if type[:3] == 'BGS':
        tsnrcol = 'TSNR2_BGS'
        dchi2 = 40
        tsnrcut = 1000

    myct.mkclusran(dirout,type,notqso,zmask=zma,tsnrcut=tsnrcut,tsnrcol=tsnrcol,ebits=None, id_=id_)

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
        ct.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
        ct.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax)
        
