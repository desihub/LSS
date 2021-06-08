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

sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.SV3.cattools as ct
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   
import LSS.mkCat_singletile.fa4lsscat as fa

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='daily')
parser.add_argument("--cuttar", help="cut targets to SV3 tiles",default='n')
parser.add_argument("--cutran", help="cut randoms to SV3 tiles",default='n')
parser.add_argument("--vis", help="make a plot of data/randoms on tile",default='n')
parser.add_argument("--xi", help="run pair-counting code",default='n')
parser.add_argument("--ranmtl", help="make a random mtl file for the tile",default='n')
parser.add_argument("--rfa", help="run randoms through fiberassign",default='n')
parser.add_argument("--combd", help="combine all the tiles together",default='y')
parser.add_argument("--combr", help="combine the random tiles together",default='n')
parser.add_argument("--dodt", help="process individual tiles; not really necessary anymore",default='n')
parser.add_argument("--redodt", help="remake already done data tiles",default='n')
parser.add_argument("--fulld", help="make the 'full' catalog containing info on everything physically reachable by a fiber",default='y')
parser.add_argument("--fullr", help="make the random files associated with the full data files",default='n')
parser.add_argument("--clus", help="make the data clustering files; these are cut to a small subset of columns",default='y')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--maskz", help="apply sky line mask to redshifts?",default='n')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=1) 

parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')

#all random set to n by default since mkCat_SV3_ran.py exists and does it in parallel

args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
version = args.version
specrel = args.verspec

zma = False
if args.maskz == 'y':
    zma = True
ctar = False
if args.cuttar == 'y':
    ctar = True
cran = False
if args.cutran == 'y':
    cran = True
docatplots = False
if args.vis == 'y':
    docatplots = True
doclus = False
if args.xi == 'y':    
    doclus = True
mkranmtl = False
if args.ranmtl == 'y':
    mkranmtl = True
runrfa = True#run randoms through fiberassign
if args.rfa == 'n':
    runrfa = False
remake_dtile = True
if args.redodt == 'n':
    remake_dtile = False

mkdtiles = False #not really necessary anymore
if args.dodt == 'y':
    mkdtiles = True

mkfulld = True #make the 'full' catalog containing info on everything physically reachable by a fiber
if args.fulld == 'n':
    mkfulld = False
mkfullr = True #make the random files associated with the full data files
if args.fullr == 'n':
    mkfullr = False
mkclus = True #make the data/random clustering files; these are cut to a small subset of columns
mkclusdat = True
mkclusran = False
if args.clus == 'n':
    mkclus = False
    mkclusdat = False
if args.clusran == 'y':
    mkclusran = True
mknz = False #get n(z) for type and all subtypes
if args.nz == 'y':
    mknz = True

fillNZ = False

combd = True
if args.combd == 'n':
    combd = False
combr = True
if args.combr == 'n':
    combr = False   

if type == 'dark' or type == 'bright':
    #set all type catalog stuff to False in this case
    mkfulld = False
    mkfullr = False
    mkclus = False
    mkclusdat = False
    mkclusran = False
    
    

if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    pr = 'BRIGHT'
    pdir = 'bright'
else:
    pr = 'DARK'
    pdir = 'dark'

pd = pdir

mdir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/sv3/'+pdir+'/' #location of ledgers
tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/targets/sv3/resolve/'+pdir+'/' #location of targets
#mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/mtl-done-tiles.ecsv') #log of tiles completed for mtl
mtld = Table.read('/global/cfs/cdirs/desi/spectro/redux/daily/tiles.csv')
wdone = mtld['ZDONE'] == 'true'
mtld = mtld[wdone]
tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')
imbits = [1,5,6,7,8,9,11,12,13]

#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
sv3dir = basedir +'/SV3/LSS/'



#tarbit = int(np.log2(sv3_targetmask.desi_mask[type]))

wp = tiles['PROGRAM'] == pr
tiles = tiles[wp]
print(len(tiles))

wp = np.isin(mtld['TILEID'],tiles['TILEID']) #we want to consider MTL done tiles that correspond to the SV3 tile file
mtld = mtld[wp]
print(len(mtld))



if not os.path.exists(sv3dir+'/logs'):
    os.mkdir(sv3dir+'/logs')
    print('made '+sv3dir+'/logs')

ldirspec = sv3dir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    print('made '+ldirspec)


if not os.path.exists(ldirspec+'LSScats'):
    os.mkdir(ldirspec+'LSScats')
    print('made '+ldirspec+'LSScats')

dirout = ldirspec+'LSScats/'+version+'/'
if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)


randir = sv3dir+'random'
rm = int(args.minr)
rx = int(args.maxr)
print(rm,rx)
#logf.write('using random files '+str(rm)+ ' through '+str(rx)+' (this is python, so max is not inclusive)\n')
for i in range(rm,rx):
    if not os.path.exists(sv3dir+'random'+str(i)):
        os.mkdir(sv3dir+'random'+str(i))
        print('made '+str(i)+' random directory')


#construct a table with the needed tile information
if len(mtld) > 0:
    tilel = []
    ral = []
    decl = []
    mtlt = []
    fal = []
    obsl = []
    pl = []
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
else:
    print('no done tiles in the MTL')


minr = 148
maxr = 274
mind = -2.5
maxd = 68
if ctar:
    tard = read_targets_in_tiles(mdir,tiles,mtl=True,isodate='2021-04-06T00:00:00') #this date should be after initial creation and before 1st update
    print('read in mtl targets')
    print('should be 0 '+str(np.unique(tard['NUMOBS'])))
    minr = np.min(tard['RA'])-1
    maxr = np.max(tard['RA'])+1
    mind = np.min(tard['DEC'])-1
    maxd = np.max(tard['DEC'])+1

    tardi = inflate_ledger(tard,tdir)
    tardi = Table(tardi)
    tardi.write(sv3dir+pdir+'_targets.fits',overwrite=True,format='fits')
    print('wrote inflated ledger target file to '+sv3dir+pdir+'_targets.fits')
    del tardi
    del tard

if cran:
    dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
    for ii in range(rm,rx):
        ranf = fitsio.read(dirrt+'/randoms-1-'+str(ii)+'.fits')
        print(len(ranf))
        #if ctar:
        wp = ranf['RA'] > minr
        wp &= ranf['RA'] < maxr
        wp &= ranf['DEC'] > mind
        wp &= ranf['DEC'] < maxd
        ranf = ranf[wp]
        print(len(ranf))
        tilesall = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')
        tilesu = unique(tilesall,keys=['RA','DEC'])                
        wi = is_point_in_desi(tilesu, ranf["RA"], ranf["DEC"])
        ranf = ranf[wi]
        fitsio.write(sv3dir+'random'+str(ii)+'/alltilesnofa.fits',ranf,clobber=True)
        print('wrote '+sv3dir+'random'+str(ii)+'/alltilesnofa.fits')

if mkranmtl:
    ct.randomtiles_allSV3(ta,imin=rm,imax=rx)
    
if runrfa:
    print('DID YOU DELETE THE OLD FILES!!!')
    for ii in range(0,len(mtld)):
        tile = mtld['TILEID'][ii]
        ts = str(tile).zfill(6)
        fbah = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
        dt = fbah['RUNDATE']
        ttemp = Table(ta[ii])
        ttemp['OBSCONDITIONS'] = 516
        ttemp['IN_DESI'] = 1
        ttemp.write('tiletemp.fits',format='fits', overwrite=True)
        for i in range(rm,rx):
            testfbaf = randir+str(i)+'/fba-'+str(tile).zfill(6)+'.fits'
            if os.path.isfile(testfbaf):
                print('fba file already made')
            else:                   
                fa.getfatiles(randir+str(i)+'/tilenofa-'+str(tile)+'.fits','tiletemp.fits',dirout=randir+str(i)+'/',dt = dt)

if mkdtiles:
    #for tile,zdate in zip(mtld['TILEID'],mtld['ZDATE']):
    for tile,zdate in zip(mtld['TILEID'],mtld['LASTNIGHT']): 
        ffd = dirout+'ALL'+str(tile)+'_full.dat.fits'
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

if combd:
    #print(len(mtld['TILEID']))
    #if type == 'dark' or type == 'bright':
    #    ct.count_tiles(mtld['TILEID'],dirout,type)
    #ct.combtiles(mtld['TILEID'],dirout,type,sv3_targetmask.desi_mask)   
    #dd = ct.combtarinfo_all(ta,mdir=mdir)
    if type == 'dark' or type == 'bright':
        outf = sv3dir+'datcomb_'+type+'_tarwdup_Alltiles.fits'
        ct.combtiles_wdup(ta,mdir,outf)
        outf = ldirspec+'datcomb_'+type+'_specwdup_Alltiles.fits'
        ct.combtile_spec(mtld,outf,rel=specrel)
        tarf = Table.read(sv3dir+'datcomb_'+type+'_tarwdup_Alltiles.fits')
        tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
        remcol = ['PRIORITY','Z','ZWARN','FIBER']
        for col in remcol:
            try:
                tarf.remove_columns([col] )#we get this where relevant from spec file
            except:
                print('column '+col +' was not in tarwdup file')    
        specf = Table.read(ldirspec+'datcomb_'+type+'_specwdup_Alltiles.fits')
        specf.keep_columns(['CHI2','COEFF','Z','ZERR','ZWARN','NPIXELS','SPECTYPE','SUBTYPE','NCOEFF','DELTACHI2'\
        ,'FIBERASSIGN_X','FIBERASSIGN_Y','TARGETID','LOCATION','FIBER','FIBERSTATUS','PRIORITY'\
        ,'DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE','NIGHT','EXPID','MJD','TILEID','INTEG_COADD_FLUX_B',\
        'MEDIAN_COADD_FLUX_B','MEDIAN_COADD_SNR_B','INTEG_COADD_FLUX_R','MEDIAN_COADD_FLUX_R','MEDIAN_COADD_SNR_R','INTEG_COADD_FLUX_Z',\
        'MEDIAN_COADD_FLUX_Z','MEDIAN_COADD_SNR_Z','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
        'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
        'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
        specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
        tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')
        try:
            print(np.unique(tj['SV3_DESI_TARGET'],return_counts=True))
        except:
            ftar = Table.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+type+'_targets.fits')
            ftar.keep_columns(['TARGETID','SV3_DESI_TARGET','SV3_BGS_TARGET','SV3_MWS_TARGET'])
            print(len(tj))
            tj = join(tj,ftar,keys=['TARGETID'])  
            print(len(tj))  
        tj.write(ldirspec+'datcomb_'+type+'_tarspecwdup_Alltiles.fits',format='fits', overwrite=True)
        tc = ct.count_tiles_better('dat',pdir,specrel=specrel)
        tc.write(ldirspec+'Alltiles_'+pdir+'_tilelocs.dat.fits',format='fits', overwrite=True)
    else:
        print('nothing to be done for combd, only done for dark/bright now')
    #outf = dirout+'datcomb_'+type+'wdup_Alltiles.fits'
    #dd.write(outf,format='fits', overwrite=True)


if combr:
    #print(len(mtld['TILEID']))
    if type == 'dark' or type == 'bright':
        for i in range(rm,rx):
            #ct.combran(mtld,i,randir,dirout,type,sv3_targetmask.desi_mask)
            ct.combran_wdup(mtld,i,randir,type,sv3dir)
            tc = ct.count_tiles_better('ran',pdir,i)
            tc.write(randir+str(i)+'/rancomb_'+pdir+'_Alltilelocinfo.fits',format='fits', overwrite=True)
    else:
        print('nothing to be done for combr, only done for dark/bright now')
        
        
        
if mkfulld:
    #ct.mkfulldat(dirout+'ALLAlltiles_'+pd+'_full.dat.fits',imbits,tdir,'SV3_DESI_TARGET',sv3_targetmask.desi_mask[type],dirout+type+'Alltiles_full.dat.fits')
    azf = '/global/homes/r/raichoor/sv3/sv3-elg-daily-thru20210521.fits'
    #/global/homes/r/raichoor/sv3/sv3-elg-daily-thru20210506.fits
    #dz = dirout+'datcomb_'+type+'_Alltiles.fits' old
    dz = ldirspec+'datcomb_'+pdir+'_tarspecwdup_Alltiles.fits' #new
    if type == 'BGS_BRIGHT':
        bit = sv3_targetmask.bgs_mask[type]
        desitarg='SV3_BGS_TARGET'
    else:
        bit = sv3_targetmask.desi_mask[type]
        desitarg='SV3_DESI_TARGET'
    print(desitarg,pdir,bit)
    ct.mkfulldat(dz,imbits,tdir,type,bit,dirout+type+'Alltiles_full.dat.fits',ldirspec+'Alltiles_'+pdir+'_tilelocs.dat.fits',azf=azf,desitarg=desitarg,specver=specrel)
    #get_tilelocweight()
    #logf.write('ran get_tilelocweight\n')
    #print('ran get_tilelocweight\n')

if mkfullr:
    for ii in range(rm,rx):
        outf = dirout+type+'Alltiles_'+str(ii)+'_full.ran.fits'
        ct.mkfullran(randir,ii,imbits,outf,type,pdir,sv3_targetmask.desi_mask[type])
    #logf.write('ran mkfullran\n')
    #print('ran mkfullran\n')

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
        tsnrcut = 1000
    ct.mkclusdat(dirout+type+'Alltiles_',zmask=zma,tp=type,dchi2=dchi2,tsnrcut=tsnrcut)
    #logf.write('ran mkclusdat\n')
    #print('ran mkclusdat\n')

if mkclusran:
    print('doing clustering randoms')
    tsnrcol = 'TSNR2_ELG'
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

    for ii in range(rm,rx):
        ct.mkclusran(dirout+type+'Alltiles_',ii,zmask=zma,tsnrcut=tsnrcut,tsnrcol=tsnrcol)
    #logf.write('ran mkclusran\n')
    #print('ran mkclusran\n')
    
if mknz:
    regl = ['','_N','_S']
    for reg in regl:
        if zma:
            reg = '_zmask'+reg
        fcr = dirout+type+'Alltiles'+reg+'_0_clustering.ran.fits'
        fcd = dirout+type+'Alltiles'+reg+'_clustering.dat.fits'
        fout = dirout+type+reg+'_nz.dat'
        if type == 'QSO':
            zmin = 0.6
            zmax = 4.5
            dz = 0.05
            ct.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
        else:    
            ct.mknz(fcd,fcr,fout,bs=0.02)

if fillNZ:
    e2e.fillNZ(target_type,program,P0=P0,truez=truez)   
    logf.write('put NZ and weight_fkp into clustering catalogs\n')    
    print('put NZ and weight_fkp into clustering catalogs\n')
        