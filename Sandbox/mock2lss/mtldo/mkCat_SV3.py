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
import LSS.SV3.cattools as ct
import LSS.mkCat_singletile.fa4lsscat as fa
from LSS.globals import SV3 
import mockcattools as mt 

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='everest')
parser.add_argument("--cuttar", help="cut targets to SV3 tiles",default='n')
parser.add_argument("--vis", help="make a plot of data/randoms on tile",default='n')
parser.add_argument("--xi", help="run pair-counting code",default='n')
parser.add_argument("--combd", help="combine all the tiles together",default='y')
parser.add_argument("--dodt", help="process individual tiles; not really necessary anymore",default='n')
parser.add_argument("--redodt", help="remake already done data tiles",default='n')
parser.add_argument("--fulld", help="make the 'full' catalog containing info on everything physically reachable by a fiber",default='y')
parser.add_argument("--clus", help="make the data clustering files; these are cut to a small subset of columns",default='y')
parser.add_argument("--maskz", help="apply sky line mask to redshifts?",default='n')
parser.add_argument("--univ", help="Which AltMTL realization?",default=1)
parser.add_argument("--isoMTL", help="isodate for initial ledger",default='2022-03-10T16:32:15.000')

parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')

parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')

#all random set to n by default since mkCat_SV3_ran.py exists and does it in parallel


#READING CONFIGURATION
##########################################################################
args = parser.parse_args()
print(args)
id_ = "%03d"%int(args.univ)
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
remake_dtile = True
if args.redodt == 'n':
    remake_dtile = False

mkdtiles = False #not really necessary anymore
if args.dodt == 'y':
    mkdtiles = True

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

if type == 'dark' or type == 'bright':
    #set all type catalog stuff to False in this case
    mkfulld = False
    mkclus = False
    mkclusdat = False
    
    

if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    pr = 'BRIGHT'
    pdir = 'bright'
else:
    pr = 'DARK'
    pdir = 'dark'

pd = pdir

mdir = '/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_016dirs/Univ{UNIV}/sv3/dark'.format(UNIV=id_)  #SV3p.mdir+pdir+'/' #location of ledgers
tdir = '/global/cscratch1/sd/acarnero/SV3/mockTargets_000_FirstGen_CutSky_alltracers_sv3bits.fits' #SV3p.tdir+pdir+'/' #location of targets
#tdir = '/global/cscratch1/sd/acarnero/SV3/atest000' #SV3p.tdir+pdir+'/' #location of targets
mtld = SV3p.mtld
tiles = SV3p.tiles
imbits = SV3p.imbits #mask bits applied to targeting
ebits = SV3p.ebits #extra mask bits we think should be applied

print('mdir',mdir)
print('tdir',tdir)


#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'

def test_dir(value):
    if not os.path.exists(value):
        try:
            os.makedirs(value, 0o755)
            print('made %s'%value)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


#sv3dir = basedir +'/SV3/LSS/'
sv3dir = os.path.join(basedir,'SV3', 'LSS_MTL_{UNIV}'.format(UNIV=args.univ))
test_dir(sv3dir)

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




#construct a table with the needed tile information
##################################################################################
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
    ttf = Table() #to write out to use for fiberassign all at once
    ttf['TILEID'] = tilel
    ttf['RA'] = ral
    ttf['DEC'] = decl
    ttf['OBSCONDITIONS'] = 15
    ttf['IN_DESI'] = 1
    ttf['PROGRAM'] = 'SV3'
    ta.write(os.path.join(sv3dir,'tiles-'+pr+'.fits'),format='fits', overwrite=True)

else:
    print('no done tiles in the MTL')


#CREATE INITIAL LEDGER FOR NEW RUN. 
#THIS IS OBSOLETE
####################################################################################
minr = 148
maxr = 274
mind = -2.5
maxd = 68
if ctar:
    tard = read_targets_in_tiles(mdir,tiles,mtl=True)#,isodate='2021-04-06T00:00:00') #this date should be after initial creation and before 1st update
    print('read in mtl targets')
    print('should be 0 '+str(np.unique(tard['NUMOBS'])))
    minr = np.min(tard['RA'])-1
    maxr = np.max(tard['RA'])+1
    mind = np.min(tard['DEC'])-1
    maxd = np.max(tard['DEC'])+1
    print('should fail here')
    tardi = inflate_ledger(tard,tdir)
    tardi = Table(tardi)
    tardi.write(os.path.join(sv3dir,pdir+'_targets.fits'),overwrite=True,format='fits')
    print('wrote inflated ledger target file to '+sv3dir+pdir+'_targets.fits')
    del tardi
    del tard

#THIS MATCH TARGET TO REDSHIFTS
#THIS IS OBSOLETE
##############################################################################################################################
if mkdtiles:
    #for tile,zdate in zip(mtld['TILEID'],mtld['ZDATE']):
    for tile,zdate in zip(mtld['TILEID'],mtld['LASTNIGHT']): 
        ffd = os.path.join(dirout,'ALL'+str(tile)+'_full.dat.fits')
        if os.path.isfile(ffd) and remake_dtile == False:
            print(ffd +' file already made and remakes not requested')
        else:
            zdate = str(zdate)
            tspec = ct.combspecdata(tile,zdate)
            pdict,goodloc = ct.goodlocdict(tspec)
            wloc = (np.isin(tspec['LOCATION'],goodloc))
            tspec = tspec[wloc]
            ts = str(tile).zfill(6)
            fbaf = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
            wt = ta['TILEID'] == tile
            tars = read_targets_in_tiles(mdir,ta[wt],mtl=True)
            
            ftar = Table.read(sv3dir+pdir+'_targets.fits')
            ftar.keep_columns(['TARGETID','EBV','FLUX_G','FLUX_R','FLUX_Z','FLUX_IVAR_G','FLUX_IVAR_R','FLUX_IVAR_Z','MW_TRANSMISSION_G','MW_TRANSMISSION_R',\
            'MW_TRANSMISSION_Z','FRACFLUX_G','FRACFLUX_R','FRACFLUX_Z','FRACMASKED_G','FRACMASKED_R','FRACMASKED_Z','FRACIN_G','FRACIN_R',\
            'FRACIN_Z','NOBS_G','NOBS_R','NOBS_Z','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','FLUX_W1',\
            'FLUX_W2','FLUX_IVAR_W1','FLUX_IVAR_W2','MW_TRANSMISSION_W1','MW_TRANSMISSION_W2','ALLMASK_G','ALLMASK_R','ALLMASK_Z','FIBERFLUX_G',\
            'FIBERFLUX_R','FIBERFLUX_Z','FIBERTOTFLUX_G','FIBERTOTFLUX_R','FIBERTOTFLUX_Z','WISEMASK_W1','WISEMASK_W2','MASKBITS',\
            'RELEASE','BRICKID','BRICKNAME','BRICK_OBJID','MORPHTYPE','PHOTSYS'])
            ol = len(tars)
            tars = join(tars,ftar,keys=['TARGETID'])
            print('lengths after join:'+str(ol),len(tars))
            #tars = inflate_ledger(tars,tdir) #need to specify columns here or MTL updates will be reversed to original state
            tars = tars[[b for b in list(tars.dtype.names) if b != 'Z']]
            tars = tars[[b for b in list(tars.dtype.names) if b != 'ZWARN']]
            tars = tars[[b for b in list(tars.dtype.names) if b != 'PRIORITY']]
            tars = join(tars,tspec,keys=['TARGETID'],join_type='left')
            tout = ct.gettarinfo_type(fbaf,tars,goodloc,pdict)
            #tout = join(tfa,tspec,keys=['TARGETID','LOCATION'],join_type='left') #targetid should be enough, but all three are in both and should be the same
            print(tout.dtype.names)
            wz = tout['ZWARN']*0 == 0
            wzg = tout['ZWARN'] == 0
            print('there are '+str(len(tout[wz]))+' rows with spec obs redshifts and '+str(len(tout[wzg]))+' with zwarn=0')
        
            tout.write(ffd,format='fits', overwrite=True) 
            print('wrote matched targets/redshifts to '+ffd)
            #logf.write('made full data files\n')

#CREATE A DICTIONARY WITH TILEID AND THE DIRECTORY OF THE ALTMTL FBA RUN
##############################################################################################
list_runFA = {}
infp = Table.read('/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_016dirs/Univ{UNIV}/mtl-done-tiles.ecsv'.format(UNIV=id_))
for tile in ta['TILEID']:
    ts = str(tile).zfill(6)
    faf_d = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
    fht = fitsio.read_header(faf_d)
    stamp = fht['RUNDATE'].split('T')[0].replace('-','')
    list_runFA[tile] = stamp


#isodate is retrieve with check_isomax.py from ALTMTL newly created files, in
#/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_016dirs/Univ000/sv3/dark/orig     2022-02-11T16:42:57.000
run_tarwdup = True
run_specwdup = True

if combd:
    if type == 'dark' or type == 'bright':
#CREATE tarwdup FILE READING FROM ALTMTL RESULTS, USING FAVAIL hdu
##########################################################################################################################
        outf = os.path.join(sv3dir,'datcomb_'+type+'_tarwdup_Alltiles.fits')
        #AUREct.combtiles_wdup(ta,mdir,outf)
        if run_tarwdup:
            #Univ001 2022-02-14T19:37:05.000
            #Univ000 2022-02-11T16:42:58.000
            mt.combtiles_wdup_mtl(ta, mdir, outf, mtl_done='/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_016dirs/Univ{UNIV}/mtl-done-tiles.ecsv'.format(UNIV=id_), univ=id_, isodate=args.isoMTL)

        tarf = Table.read(outf)
        tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
        remcol = ['PRIORITY','Z','ZWARN','SUBPRIORITY'] #subpriority in target files doesn't match what is in fiberassign files
#AURE        remcol = ['PRIORITY','Z','ZWARN','FIBER','SUBPRIORITY'] #subpriority in target files doesn't match what is in fiberassign files
        for col in remcol:
            try:
                tarf.remove_columns([col] )#we get this where relevant from spec file
            except:
                print('column '+col +' was not in tarwdup file')    
        print('hasta aquio')
        
#FOR NEWLY CREATED tarwdup FILE, JOIN WITH REAL DATA TO RETRIEVE HARDWARE AND TEMPLATE INFO BASED ON TILEID AND LOCATION, including coadd_fiberstatus
################################################################################################################
        data_specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-sv3-'+type+'-cumulative.fits')
        wt = np.isin(data_specf['TILEID'],ta['TILEID']) #cut spec file to dark or bright time tiles
        data_specf = data_specf[wt]

        data_specf.keep_columns(['LOCATION','COADD_FIBERSTATUS','TILEID','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT',\
            'TSNR2_ELG_B','TSNR2_LYA_B',\
            'TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
            'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
            'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])

        tarf = join(tarf, data_specf, keys=['LOCATION','TILEID'], join_type='left')
        
        if specrel == 'everest':

#CREATE specwdup FILE READING fba RESULT, USING FASSIGN hdu, MATCH TO masterTarget file
#######################################################################################################################################
            outfile_spec = os.path.join(ldirspec, 'datcomb_'+type+'_specwdup_Alltiles.fits')

            namecomb = os.path.join('/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_016dirs/Univ{UNIV}/fa/SV3'.format(UNIV=id_),'{stamp}','fba-{ts}.fits') 
            if run_specwdup:
                mt.combtile_specmock_mtl(ta, namecomb, list_runFA, tdir, outfile_spec)

            specf = Table.read(outfile_spec)
            specf.keep_columns(['FIBER','TARGETID','LOCATION','FIBERSTATUS','LAMBDA_REF','PETAL_LOC','DEVICE_LOC','DEVICE_TYPE','TARGET_RA','TARGET_DEC','FA_TARGET','FA_TYPE','FIBERASSIGN_X','FIBERASSIGN_Y','TILEID','PRIORITY','SUBPRIORITY','ZWARN','TRUEZ','RSDZ'])
            specf['TILELOCID'] = 10000*specf['TILEID'] + specf['LOCATION']

            ##AUREtj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID'],join_type='left')
#JOIN TARGET WITH DUP TO SPEC WITH DUP 
####################################################################################################################
            tj = join(tarf, specf, keys=['TARGETID','LOCATION','TILEID','TILELOCID'], join_type='left')

            tj.write(os.path.join(ldirspec,'datcomb_'+type+'_tarspecwdup_Alltiles.fits'),format='fits', overwrite=True)
        
        print(np.unique(tj['SV3_DESI_TARGET'],return_counts=True))

#Make tilelocs file with NTILE per TARGETID
######################################################################################
        tc = mt.count_tiles_better_mtl(data_specf, os.path.join(ldirspec,'datcomb_'+type+'_tarspecwdup_Alltiles.fits'), pdir, specrel=specrel)
#AURE        tc = ct.count_tiles_better(specf, 'dat', pdir, specrel=specrel)
        tc.write(os.path.join(ldirspec,'Alltiles_'+pdir+'_tilelocs.dat.fits'), format='fits', overwrite=True)
    else:
        print('nothing to be done for combd, only done for dark/bright now')


        
        
if mkfulld:
    if specrel == 'everest':
        outfile_spec = os.path.join(ldirspec, 'datcomb_'+pdir+'_specwdup_Alltiles.fits')
        specf = Table.read(outfile_spec)
##AURE        specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-sv3-'+pdir+'-cumulative.fits')
##AURE        wt = np.isin(specf['TILEID'],ta['TILEID']) #cut spec file to dark or bright time tiles
##AURE        specf = specf[wt]
    if specrel == 'daily':
        specf = Table.read(ldirspec+'datcomb_'+type+'_specwdup_Alltiles.fits')

    #ct.mkfulldat(dirout+'ALLAlltiles_'+pd+'_full.dat.fits',imbits,tdir,'SV3_DESI_TARGET',sv3_targetmask.desi_mask[type],dirout+type+'Alltiles_full.dat.fits')
    azf=''
    if type[:3] == 'ELG':
        azf = SV3p.elgzf#'/global/cfs/cdirs/desi/users/raichoor/everest/sv3-elg-everest-tiles.fits'
    if type[:3] == 'QSO':
        azf = SV3p.qsozf#'/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/QSO/QSO_catalog_SV3.fits'
    #'/global/homes/r/raichoor/sv3/sv3-elg-daily-thru20210521.fits'
    #/global/homes/r/raichoor/sv3/sv3-elg-daily-thru20210506.fits
    #dz = dirout+'datcomb_'+type+'_Alltiles.fits' old
    dz = os.path.join(ldirspec,'datcomb_'+pdir+'_tarspecwdup_Alltiles.fits') #new
    print(dz)
    if type == 'BGS_BRIGHT':
        bit = sv3_targetmask.bgs_mask[type]
        desitarg='SV3_BGS_TARGET'
    else:
        bit = sv3_targetmask.desi_mask[type]
        desitarg='SV3_DESI_TARGET'
    print(desitarg,pdir,bit)
    bitweightfile = None
    bitweightfile = '/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_016dirs/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
    '''AURE
    if pdir == 'dark':
        bitweightfile = SV3p.darkbitweightfile
        #bitweightfile='/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run64_2/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
        #bitweightfile='/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run64_2/BitweightsRound2/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
    if pdir == 'bright':
        bitweightfile = SV3p.brightbitweightfile
    '''

##    mt.mkfulldat_mtl(specf,dz,imbits,tdir,type,bit,os.path.join(dirout,type+notqso+'_full_noveto.dat.fits'),os.path.join(ldirspec,'Alltiles_'+pdir+'_tilelocs.dat.fits'),azf=azf,desitarg=desitarg,specver=specrel,notqso=notqso,bitweightfile=bitweightfile,fbcol='FIBERSTATUS')
    mt.mkfulldat_mtl(specf,dz,imbits,tdir,type,bit,os.path.join(dirout,type+notqso+'_full_noveto.dat.fits'),os.path.join(ldirspec,'Alltiles_'+pdir+'_tilelocs.dat.fits'),azf=azf,desitarg=desitarg,specver=specrel,notqso=notqso,bitweightfile=bitweightfile,otherspec='/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-sv3-'+pdir+'-cumulative.fits')#,fbcol='FIBERSTATUS')
    #get_tilelocweight()
    #logf.write('ran get_tilelocweight\n')
    #print('ran get_tilelocweight\n')

###weightmd='probobs'
weightmd='tileloc'
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
    mt.mkclusdat_mtl(os.path.join(dirout,type+'_'), weightmd='tileloc', zmask=zma,tp=type,dchi2=dchi2,tsnrcut=tsnrcut,ebits=None)
##AURE    ct.mkclusdat(os.path.join(dirout,type+'_'),zmask=zma,tp=type,dchi2=dchi2,tsnrcut=tsnrcut,ebits=ebits)
    #logf.write('ran mkclusdat\n')
    #print('ran mkclusdat\n')

rcut = None 
if mknz:
    wzm = ''
#     if zmask:
#         wzm = 'zmask_'
    if rcut is not None:
        wzm += '_rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
#    if ntile > 0:
#        wzm += '_ntileg'+str(ntilecut)+'_'    
#    if args.ccut is not None:
#        wzm += '_'+args.ccut #you could change this to however you want the file names to turn out

    regl = ['']#,'_N','_S']
    
    for reg in regl:
        fb = os.path.join(dirout,type+wzm+reg)
        
        fbr = '/global/cscratch1/sd/acarnero/SV3/LSS_MTL/everest/LSScats/test/'+type+wzm+reg
##with data real        fbr = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/LSScats/2.1/'+type+wzm+reg
        
        fcr = fbr+'_0_clustering.ran.fits'
#AURE        fcr = fb+'_0_clustering.ran.fits'
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
        mt.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
#        mt.addnbar(fb,fbr,bs=dz,zmin=zmin,zmax=zmax)
        
