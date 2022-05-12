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

#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.SV3.cattools as ct
import LSS.common_tools as common
import LSS.mkCat_singletile.fa4lsscat as fa
from LSS.globals import SV3 
import myfa
import mocktools as mt
import targettools as tt
from desitarget.sv3 import sv3_targetmask

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='fuji')
parser.add_argument("--cutran", help="cut randoms to SV3 tiles",default='n')
parser.add_argument("--ranmtl", help="make a random mtl file for the tile",default='n')
parser.add_argument("--rfa", help="run randoms through fiberassign",default='n')
parser.add_argument("--combr", help="combine the random tiles together",default='n')
parser.add_argument("--fullr", help="make the random files associated with the full data files",default='n')
parser.add_argument("--apply_veto", help="make the random files associated with the full data files",default='n')
parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')
parser.add_argument("--maskz", help="apply sky line mask to redshifts?",default='n')
parser.add_argument("--faver", help="version of fiberassign code to use for random; versions for SV3 are '2.3.0' '2.4.0' '2.5.0' '2.5.1' '3.0.0' '4.0.0'",default='2.3.0')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=20) 
parser.add_argument("--par", help="run different random number in parallel?",default='y')
parser.add_argument("--clus", help="make the data/random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
parser.add_argument("--univ", help="Which AltMTL realization?",default=1)
parser.add_argument("--mockrea", help="Which mock realization?",default=0)

args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
version = args.version
faver = args.faver
specrel = args.verspec
rm = int(args.minr)
rx = int(args.maxr)
par = False
if args.par == 'y':
    par = True


zma = False
if args.maskz == 'y':
    zma = True


cran = False
if args.cutran == 'y':
    cran = True

mkranmtl = False
if args.ranmtl == 'y':
    mkranmtl = True
runrfa = True#run randoms through fiberassign
if args.rfa == 'n':
    runrfa = False
combr = True
if args.combr == 'n':
    combr = False   


mkfullr = True #make the random files associated with the full data files
if args.fullr == 'n':
    mkfullr = False

mkclusran = True
if args.clus == 'n':
    mkclusran = False

if type == 'bright' or type == 'dark':
    #don't do any of the actual catalogs in this case
    mkclusran = False
    mkfullr = False

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'


if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    pr = 'BRIGHT'
    pdir = 'bright'
else:
    pr = 'DARK'
    pdir = 'dark'

pd = pdir

SV3p = SV3(type, specver=specrel)

id_ = "%03d"%int(args.univ)
mockrea = "%03d"%int(args.mockrea)

release = 'sv3'

pmdir = '/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_256dirs_rea{MOCKREA}/Univ{UNIV}'.format(MOCKREA=mockrea, UNIV=id_)
mdir = os.path.join(pmdir, release, pd)

tdir = '/global/cscratch1/sd/acarnero/SV3/mockTargets_{MOCKREA}_FirstGen_CutSky_alltracers_sv3bits.fits'.format(MOCKREA=mockrea)
mtld = SV3p.mtld
tiles = SV3p.tiles
imbits = SV3p.imbits #mask bits applied to targeting
ebits = SV3p.ebits #extra mask bits we think should be applied



#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
sv3dir = os.path.join(basedir,'SV3', 'LSS_rea{MOCKREA}_univ{UNIV}'.format(MOCKREA=mockrea, UNIV=args.univ))
mt.test_dir(sv3dir)


#tarbit = int(np.log2(sv3_targetmask.desi_mask[type]))

wp = tiles['PROGRAM'] == pr
tiles = tiles[wp]
print(len(tiles))

wp = np.isin(mtld['TILEID'],tiles['TILEID']) #we want to consider MTL done tiles that correspond to the SV3 tile file
mtld = mtld[wp]
print(len(mtld))


mt.test_dir(os.path.join(sv3dir,'logs'))

ldirspec = os.path.join(sv3dir, specrel)
mt.test_dir(ldirspec)

mt.test_dir(os.path.join(ldirspec,'LSScats'))

dirout = os.path.join(ldirspec,'LSScats', version)
mt.test_dir(dirout)

randir = os.path.join(sv3dir,'random')

for i in range(rm,rx):
    mt.test_dir(randir+str(i))
    print('made '+str(i)+' random directory')

#READ INFORMATION ABOUT TILES,
tilef = os.path.join(sv3dir,'tiles-'+pr+'.fits')
ta = mt.read_info_tiles(tilef, mtld, pr)

''' I LEAVE THIS HERE FOR THE FAVER VERSION WHICH IS DIFFERENT THAN IN DATA
#construct a table with the needed tile information
if len(mtld) > 0:
    tilel = []
    ral = []
    decl = []
    mtlt = []
    fal = []
    obsl = []
    pl = []
    fver = []
    fahal = []
    #for tile,pro in zip(mtld['TILEID'],mtld['PROGRAM']):
    for tile in mtld['TILEID']:
        ts = str(tile).zfill(6)
        fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
        tilel.append(tile)
        ral.append(fht['TILERA'])
        decl.append(fht['TILEDEC'])
        mtlt.append(fht['MTLTIME'])
        fal.append(fht['FA_RUN'])
        obsl.append(fht['OBSCON'])
        fav = fht['FA_VER']
        if np.isin(fav,['2.2.0.dev2811','2.3.0','2.3.0.dev2838']):#2.3.0 confirmed to work for these
            fver.append('2.3.0')
        else:
            fver.append(fav)    
        #try:
        #    faha = fht['FA_HA']
        #except:
        #    faha = 0
        #    print(tile,'no FA_HA in this tile header')        
        #pl.append(pro)
        pl.append(pr)
    ta = Table()
    ta['TILEID'] = tilel
    ta['RA'] = ral
    ta['DEC'] = decl
    ta['MTLTIME'] = mtlt
    ta['FA_RUN'] = fal
    ta['OBSCON'] = obsl
    ta['PROGRAM'] = pl
    #ta['FA_HA'] = fahal
    #ta['FA_VER'] = fver
    print(np.unique(fver))
    wfv = (np.array(fver) == faver)
    #mtld =  mtld[wfv]
    #ta = ta[wfv]
else:
    print('no done tiles in the MTL')
'''


ran_ids = np.linspace(100,5000,50)

#CREATE A DICTIONARY WITH TILEID AND THE DIRECTORY OF THE ALTMTL FBA RUN
##############################################################################################
list_runFA = mt.create_tile_altmtldir(mockrea, id_, ta)


def doran(ii):
    dirrt='/global/cscratch1/sd/acarnero/SV3'
    if cran:
        ranfile = os.path.join(dirrt, 'mockRandom_{RANID}_FirstGen_CutSky_alltracers_sv3bits.fits'.format(RANID=int(ran_ids[ii])))
        ranf = fitsio.read(ranfile)

        print('size parent catalog', len(ranf))

        tilesall = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')
        tilesu = unique(tilesall,keys=['RA','DEC'])                
        wi = is_point_in_desi(tilesu, ranf["RA"], ranf["DEC"])
        ranf = ranf[wi]
        fitsio.write(os.path.join(randir+str(ii),'alltilesnofa.fits'), ranf, clobber=True)
        print('wrote '+sv3dir+'random'+str(ii)+'/alltilesnofa.fits')

    if mkranmtl:
        mt.randomtiles_allSV3(ta, os.path.join(randir+str(ii),'alltilesnofa.fits'), directory_output=randir+str(ii))

    if runrfa:
        for it in range(0,len(mtld)):
            #print(it,len(mtld))
            tile = mtld['TILEID'][it]
            ts = str(tile).zfill(6)

            testfbaf = os.path.join(randir+str(ii),'fba-'+str(tile).zfill(6)+'.fits')
            if os.path.isfile(testfbaf):
                print('fba file already made')
            else:
                stamp = list_runFA[tile]
                myfa.dofa(os.path.join(randir+str(ii),'tilenofa-'+str(tile)+'.fits'),ts,stamp,randir+str(ii),id_,mockrea)
    '''
    if runrfa:
        print('DID YOU DELETE THE OLD FILES!!!')
        for it in range(0,len(mtld[wfv])):
            #print(it,len(mtld))    
            tile = mtld['TILEID'][it]
            ts = str(tile).zfill(6)
            fbah = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
            dt = fbah['RUNDATE'][:19]
            fav = fbah['FA_VER']
            if np.isin(fav,['2.2.0.dev2811','2.3.0','2.3.0.dev2838']):#2.3.0 confirmed to work for these
                fav = '2.3.0'
            if fav == faver:
                ttemp = Table(ta[wfv][it])
                ttemp['OBSCONDITIONS'] = 516
                ttemp['IN_DESI'] = 1
                try:
                    ttemp['FA_PLAN'] = fbah['FA_PLAN']
                    ttemp['FA_HA'] = fbah['FA_HA']
                    ttemp['FIELDROT'] = fbah['FIELDROT']
                except:
                    print('did not add FA_PLAN and FIELDROT')
                #for i in range(rm,rx):
                testfbaf = randir+str(ii)+'/fba-'+str(tile).zfill(6)+'.fits'
                if os.path.isfile(testfbaf):
                    print('fba file already made')
                else:                   
                    print(ttemp)
                    print(fav,dt)
                    ttemp.write('tiletemp'+str(ii)+'.fits',format='fits', overwrite=True)
                    fa.getfatiles(randir+str(ii)+'/tilenofa-'+str(tile)+'.fits','tiletemp'+str(ii)+'.fits',dirout=randir+str(ii)+'/',dt = dt,faver=faver)

    '''

    if combr:
        print(len(mtld['TILEID']))
        if type == 'dark' or type == 'bright':
            if specrel == 'everest' or specrel == 'fuji':
                specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/'+specrel+'/zcatalog/ztile-sv3-'+type+'-cumulative.fits')
                wt = np.isin(specf['TILEID'],ta['TILEID']) #cut spec file to dark or bright time tiles
                specf = specf[wt]
                specf = mt.add_tilelocid(specf)
                
                kc = ['ZWARN','LOCATION','FIBER','COADD_FIBERSTATUS','TILEID','TILELOCID','FIBERASSIGN_X','FIBERASSIGN_Y','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT'\
                ,'MEAN_DELTA_X','MEAN_DELTA_Y','RMS_DELTA_X','RMS_DELTA_Y','MEAN_PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B'\
                ,'TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
                'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
                'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','PRIORITY']

            if specrel == 'daily':
                specf = Table.read(ldirspec+'datcomb_'+type+'_specwdup_Alltiles.fits')
                kc = ['ZWARN','LOCATION','TILEID','TILELOCID','FIBERSTATUS','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY','DELTA_X','DELTA_Y','EXPTIME','PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B','TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z','TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG']

            mt.combran_wdup(mtld,ii,randir,type,ldirspec,specf,keepcols=kc)
            tc = mt.count_tiles_better(specf, os.path.join(ldirspec,'rancomb_'+str(ii)+type+'wdupspec_Alltiles.fits'), type,rann=ii,specrel=specrel)  
            ##AUREtc = mt.count_tiles_better(specf,'ran',type,rann=ii,specrel=specrel)  
            tc.write(os.path.join(ldirspec,'rancomb_'+str(ii)+type+'_Alltilelocinfo.fits'), format='fits', overwrite=True)


        
    if mkfullr:
        if specrel == 'everest' or specrel == 'fuji':
            specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/'+specrel+'/zcatalog/ztile-sv3-'+pdir+'-cumulative.fits')
            fbcol = 'COADD_FIBERSTATUS'
        if specrel == 'daily':
            specf = Table.read(ldirspec+'datcomb_'+pdir+'_specwdup_Alltiles.fits')
            fbcol = 'FIBERSTATUS'

        outf = os.path.join(dirout,type+notqso+'_'+str(ii)+'_full_noveto.ran.fits')
        if type == 'BGS_BRIGHT':
            bit = sv3_targetmask.bgs_mask[type]
            desitarg='SV3_BGS_TARGET'
        else:
            bit = sv3_targetmask.desi_mask[type]    
            desitarg='SV3_DESI_TARGET'
        maxp = 3400
        #maxp = 103400
        if type[:3] == 'LRG' or notqso == 'notqso':
            maxp = 3200
            #maxp = 103200
        if type[:3] == 'ELG' and notqso == 'notqso':
            maxp = 3100
            #maxp = 103100
        if type[:3] == 'BGS':
            maxp = 102100

        tt.mkfullran(specf, ldirspec, randir+str(ii), ii, imbits, outf, type, pdir, bit, desitarg=desitarg, fbcol=fbcol, notqso=notqso, maxp=maxp)
    #logf.write('ran mkfullran\n')
    #print('ran mkfullran\n')
    if args.apply_veto == 'y':
        print('applying vetos')
#        maxp = 3400
        maxp = 103400
        if type[:3] == 'LRG' or notqso == 'notqso':
            #maxp = 3200
            maxp = 103200
        if type[:3] == 'ELG' and notqso == 'notqso':
            #maxp = 3100
            maxp = 103100
        if type[:3] == 'BGS':
            maxp = 102100
        fin = os.path.join(dirout,type+notqso+'_'+str(ii)+'_full_noveto.ran.fits')
        fout = os.path.join(dirout,type+notqso+'_'+str(ii)+'_full.ran.fits')
        tt.apply_veto(fin,fout,ebits=ebits,zmask=False,maxp=maxp)


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
        rcols=['Z','WEIGHT']
        if type[:3] == 'BGS':
            rcols.append('flux_r_dered')

        tt.mkclusran(os.path.join(dirout,type+notqso+'_'), ii, zmask=zma, tsnrcut=tsnrcut, tsnrcol=tsnrcol, ebits=ebits, rcols=rcols)
        #ct.mkclusran(dirout+type+'Alltiles_',ii,zmask=zma)
    #logf.write('ran mkclusran\n')
    #print('ran mkclusran\n')
    
if __name__ == '__main__':
    if par:
        from multiprocessing import Pool
#        import sys
        #N = int(sys.argv[2])
        N = rx
        p = Pool(N)
        inds = []
        for i in range(0,N):
            inds.append(i)
        p.map(doran,inds)
        p.close()
        p.join()
    else:
        for i in range(rm,rx):
            doran(i)
