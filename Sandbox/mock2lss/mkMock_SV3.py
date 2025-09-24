#standard python
import os
import numpy as np
import argparse
from astropy.table import Table,join,unique,vstack
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desimodel.footprint import is_point_in_desi
from desitarget.sv3 import sv3_targetmask

#from this package
import LSS.common_tools as common
from LSS.globals import SV3 
import mocktools as mt
import targettools as tt

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='fuji')
parser.add_argument("--cuttar", help="cut targets to SV3 tiles",default='n')
parser.add_argument("--combd", help="combine all the tiles together",default='n')
parser.add_argument("--fulld", help="make the 'full' catalog containing info on everything physically reachable by a fiber",default='n')
parser.add_argument("--apply_veto", help="apply vetos for imaging, priorities, and hardware failures",default='n')
parser.add_argument("--clus", help="make the data clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--maskz", help="apply sky line mask to redshifts?",default='n')

parser.add_argument("--nz", help="get n(z) for type and all subtypes",default='n')

parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
parser.add_argument("--ntile",help="add any constraint on the number of overlapping tiles",default=0,type=int)
parser.add_argument("--rcut",help="add any cut on the rosette radius, use string like rmin,rmax",default=None)
parser.add_argument("--ccut",help="add some extra cut based on target info; should be string that tells cattools what to ",default=None)

#SV3 mock specific
parser.add_argument("--univ", help="Which AltMTL realization?",default=1)
parser.add_argument("--mockrea", help="Which Mock realization",default=0)
parser.add_argument("--isoMTL", help="isodate for initial ledger", default='2022-03-10T16:32:15.000') 

#all random set to n by default since mkCat_SV3_ran.py exists and does it in parallel

args = parser.parse_args()
print(args)

release = 'sv3'

type = args.type
basedir = args.basedir
version = args.version
specrel = args.verspec

ntile = args.ntile
rcut = args.rcut
if rcut is not None:
    rcutstr = rcut.split(',')
    rcut = []
    rcut.append(float(rcutstr[0]))
    rcut.append(float(rcutstr[1]))

ccut = args.ccut

SV3p = SV3(type,specver=specrel)

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'

zma = False
if args.maskz == 'y':
    zma = True

ctar = False
if args.cuttar == 'y':
    ctar = True

mkfulld = True #make the 'full' catalog containing info on everything physically reachable by a fiber
if args.fulld == 'n':
    mkfulld = False

mkclusdat = True

if args.clus == 'n':
    mkclusdat = False

mknz = False #get n(z) for type and all subtypes
if args.nz == 'y':
    mknz = True

combd = True
if args.combd == 'n':
    combd = False

if type == 'dark' or type == 'bright':
    #set all type catalog stuff to False in this case
    mkfulld = False
    mkclusdat = False
    
if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    pr = 'BRIGHT'
    pdir = 'bright'
else:
    pr = 'DARK'
    pdir = 'dark'

pd = pdir

id_ = "%03d"%int(args.univ)
mockrea = "%03d"%int(args.mockrea)

pmdir = '/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_256dirs_rea{MOCKREA}/Univ{UNIV}'.format(MOCKREA=mockrea, UNIV=id_)
mdir = os.path.join(pmdir, release, pd)
tdir = '/global/cscratch1/sd/acarnero/SV3/mockTargets_{MOCKREA}_FirstGen_CutSky_alltracers_sv3bits.fits'.format(MOCKREA=mockrea) 
mtld = SV3p.mtld
tiles = SV3p.tiles
imbits = SV3p.imbits #mask bits applied to targeting
ebits = SV3p.ebits #extra mask bits we think should be applied


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

#READ INFORMATION ABOUT TILES,
tilef = os.path.join(sv3dir,'tiles-'+pr+'.fits')

ta = mt.read_info_tiles(tilef, mtld, pr)


#CREATE A DICTIONARY WITH TILEID AND THE DIRECTORY OF THE ALTMTL FBA RUN
##############################################################################################
list_runFA = mt.create_tile_altmtldir(mockrea, id_, ta)

print('Running for mock ', mockrea,' Universe ', id_)
print('Location of ledgers ', mdir)
print('Location of input target', tdir)

filename_tarspecwdup = os.path.join(ldirspec, 'datcomb_'+type+'_tarspecwdup_Alltiles.fits')
outfile_spec = os.path.join(ldirspec, 'datcomb_'+type+'_specwdup_Alltiles.fits')

if combd:
    
    if type == 'dark' or type == 'bright':

        print('Have you corrected the isodate to be read in comb data?', args.isoMTL)
        print('First create datcomb_'+type+'_specwdup_Alltiles.fits')
        namecomb = os.path.join(pmdir,'fa', release.upper(), '{stamp}','fba-{ts}.fits')
        if os.path.isfile(outfile_spec):
            mock_fassign = Table.read(outfile_spec)
        else:
            mock_fassign = mt.combtile_specmock(ta, namecomb, list_runFA, tdir, outfile_spec)

    
        print('Second create datcomb_'+type+'_tarwdup_Alltiles.fits')
        outf = os.path.join(sv3dir,'datcomb_'+type+'_tarwdup_Alltiles.fits')
        if os.path.isfile(outf):
            tarf = Table.read(outf)
        else:
            tarf = mt.combtiles_wdup(ta, mdir=mdir, fout=outf, mtl_done=os.path.join(pmdir,'mtl-done-tiles.ecsv'), univ=id_, isodate=args.isoMTL, mockrea=mockrea)
        


        remcol = ['PRIORITY','Z', 'ZTILEID', 'FIBER','SUBPRIORITY'] #subpriority in target files doesn't match what is in fiberassign files
        for col in remcol:
            try:
                tarf.remove_columns([col] )#we get this where relevant from spec file
            except:
                print('column '+col +' was not in tarwdup file')    

        #JOIN TO GET HARDWARE STATUS FROM REAL OBSERVATIONS
        real_specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/'+specrel+'/zcatalog/ztile-sv3-'+type+'-cumulative.fits')
        wt = np.isin(real_specf['TILEID'], ta['TILEID']) #cut spec file to dark or bright time tiles
        real_specf = real_specf[wt]

        real_specf.keep_columns(['LOCATION','COADD_FIBERSTATUS','TILEID','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT',\
            'TSNR2_ELG_B','TSNR2_LYA_B',\
            'TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
            'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
            'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])

        if os.path.isfile(filename_tarspecwdup):
            pass
        else:
            tarf = join(tarf, real_specf, keys=['LOCATION','TILEID'], join_type='left')


        #JOIN TO GET tarspecwdup joining tarf with mock_fassign
            mock_fassign.keep_columns(['TARGETID','FIBER','LOCATION','TILEID','TILELOCID',
                'TRUEZ','RSDZ','PRIORITY','SUBPRIORITY','ZWARN','FIBERSTATUS'])

            tj = join(tarf, mock_fassign, keys=['TARGETID','LOCATION','TILEID','TILELOCID'], join_type='left')
            print(np.unique(tj['SV3_DESI_TARGET'],return_counts=True))

            tj.write(filename_tarspecwdup, format='fits', overwrite=True)
        
        #Count completeness on real data
        tc = mt.count_tiles_better(real_specf, filename_tarspecwdup, pdir, specrel=specrel)
        tc.write(os.path.join(ldirspec, 'Alltiles_'+pdir+'_tilelocs.dat.fits'), format='fits', overwrite=True)


    else:
        print('nothing to be done for combd, only done for dark/bright now')


       
        
        
if mkfulld:
    #    mock_fassign = Table.read(outfile_spec)  #This is old specf
    #specf = Table.read(outfile_spec)

    '''
    if specrel == 'everest' or specrel == 'fuji':
        specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/'+specrel+'/zcatalog/ztile-sv3-'+pdir+'-cumulative.fits')
        wt = np.isin(specf['TILEID'],ta['TILEID']) #cut spec file to dark or bright time tiles
        specf = specf[wt]
    if specrel == 'daily':
        specf = Table.read(ldirspec+'datcomb_'+type+'_specwdup_Alltiles.fits')
    '''
    azf=''
    if type[:3] == 'ELG':
        azf = SV3p.elgzf#'/global/cfs/cdirs/desi/users/raichoor/everest/sv3-elg-everest-tiles.fits'
    if type[:3] == 'QSO':
        azf = SV3p.qsozf#'/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/QSO/QSO_catalog_SV3.fits'
    #'/global/homes/r/raichoor/sv3/sv3-elg-daily-thru20210521.fits'
    #/global/homes/r/raichoor/sv3/sv3-elg-daily-thru20210506.fits
    #dz = dirout+'datcomb_'+type+'_Alltiles.fits' old
#    dz = ldirspec+'datcomb_'+pdir+'_tarspecwdup_Alltiles.fits' #new
    print(filename_tarspecwdup)
    if type == 'BGS_BRIGHT':
        bit = sv3_targetmask.bgs_mask[type]
        desitarg='SV3_BGS_TARGET'
    else:
        bit = sv3_targetmask.desi_mask[type]
        desitarg='SV3_DESI_TARGET'
    print(desitarg,pdir,bit)
    bitweightfile = os.path.join('/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_256dirs_rea{MOCKREA}/BitweightFiles'.format(MOCKREA=mockrea), release, pdir, 'sv3bw-dark-AllTiles.fits')

    filename_forfull_tarspecwdup = filename_tarspecwdup.replace(type, pdir)
    tt.mkfulldat(filename_forfull_tarspecwdup, imbits, tdir, type, bit, os.path.join(dirout,type+notqso+'_full_noveto.dat.fits'), os.path.join(ldirspec,'Alltiles_'+pdir+'_tilelocs.dat.fits'), azf=azf, desitarg=desitarg, specver=specrel, notqso=notqso, bitweightfile=bitweightfile)
##    mt.mkfulldat(mock_fassign, filename_tarspecwdup, imbits,tdir,type,bit,os.path.join(dirout,type+notqso+'_full_noveto.dat.fits'),os.path.join(ldirspec,'Alltiles_'+pdir+'_tilelocs.dat.fits'),azf=azf,desitarg=desitarg,specver=specrel,notqso=notqso,bitweightfile=bitweightfile)


if args.apply_veto == 'y':
    print('applying vetos')
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
    fin = os.path.join(dirout, type+notqso+'_full_noveto.dat.fits')
    fout = os.path.join(dirout, type+notqso+'_full.dat.fits')
    tt.apply_veto(fin, fout, ebits=ebits, zmask=False, maxp=maxp)

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
    tt.mkclusdat(os.path.join(dirout, type+'_'), tp=type, dchi2=dchi2, tsnrcut=tsnrcut, rcut=rcut, ntilecut=ntile, ccut=ccut, weightmd='probobs', ebits=ebits)
    #ct.mkclusdat(os.path.join(dirout, type+'_'), tp=type, dchi2=dchi2, tsnrcut=tsnrcut, rcut=rcut, ntilecut=ntile, ccut=ccut, weightmd=SV3p.weightmode, ebits=ebits)
    #logf.write('ran mkclusdat\n')
    #print('ran mkclusdat\n')

    
#changed to be done at same time as clustering catalogs within mkclusdat
if mknz:
    wzm = ''
#     if zmask:
#         wzm = 'zmask_'
    if rcut is not None:
        wzm += '_rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntile > 0:
        wzm += '_ntileg'+str(ntilecut)+'_'    
    if ccut is not None:
        wzm += '_'+ccut #you could change this to however you want the file names to turn out

    regl = ['','_N','_S']
    
    
    if type[:3] == 'QSO':
        zmin = 0.6
        zmax = 4.5
        dz = 0.05
        P0 = 6000
        
    else:    
        dz = 0.02
        zmin = 0.01
        zmax = 1.61
    
    if type[:3] == 'LRG':
        P0 = 10000
    if type[:3] == 'ELG':
        P0 = 4000
    if type[:3] == 'BGS':
        P0 = 7000
    
    for reg in regl:
        fb = dirout+type+notqso+wzm+reg
        fcr = fb+'_0_clustering.ran.fits'
        fcd = fb+'_clustering.dat.fits'
        fout = fb+'_nz.dat'
        common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
        common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0)


