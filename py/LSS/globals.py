from astropy.table import Table,join
import numpy as np

class SV3:
    def __init__(self,tp,weightmode='probobs',specver='fuji'):
        self.mdir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/sv3/' #location of ledgers
        self.tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/targets/sv3/resolve/'#location of targets
        self.mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
        self.tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')

        ebits = None
        #these are all the ELG parameters, they get changed for other tracers below
        self.tsnrcut = 80
        self.dchi2 = 0.9 #used for the ELG OII criteria
        self.tsnrcol = 'TSNR2_ELG'
        self.zmin = 0.6
        self.zmax = 1.6
        if tp[:3] == 'BGS':
            self.imbits = [1,13]
            self.tsnrcut = 1000
            self.dchi2 = 40
            self.tsnrcol = 'TSNR2_BGS'
            self.zmin = 0.01
            self.zmax = 0.6
        else:
            self.imbits = [1,12,13]
        self.ebits = None
        if tp[:3] == 'QSO':
            self.ebits = [8,9,11]    
            self.tsnrcut = 0
            self.dchi2 = 0
            self.zmin = 0.6
            self.zmax = 3.5
            #self.tsnrcol = 'TSNR2_QSO'
        if tp[:3] == 'LRG':
            self.ebits = 'lrg_mask'
            self.dchi2 = 15
            self.zmin = 0.4
            self.zmax = 1.1
        if tp[:3] == 'ELG' or tp[:3] == 'BGS':
            self.ebits = [11]    
        if specver == 'everest':
            self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/ELG/sv3-elg-everest-tiles.fits'
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/QSO/QSO_catalog_SV3.fits'
        if specver == 'fuji':
            self.elgzf = '/global/cfs/cdirs/desi/users/raichoor/spectro/fuji/sv3-elg-fuji-tiles.fits'
            self.elgzfhp = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/fuji/emline_darkallhealpix.fits'
            self.qsozf = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/fuji/QSO_cat_fuji_cumulative.fits'
            self.qsozfhp = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/fuji/healpix/QSO_cat_fuji_sv3_dark_healpix_only_qso_targets.fits'
            self.badfib = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/fuji/lrg+bgs_3sig_bad_fibers.txt')
        if specver == 'newQSOtemp_tagged':
            self.qsozf = '/global/cfs/cdirs/desi/users/abrodze/QSOtemplates/rroutput_ALL-fuji/catalogs/QSO_cat_fuji_cumulative_only_qso_targets.fits'
            self.badfib = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/fuji/lrg+bgs_3sig_bad_fibers.txt')
        self.darkbitweightfile = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/DESI_EDA_SV3AltMTLs/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
        #self.darkbitweightfile = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run128/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
        self.brightbitweightfile = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/DESI_EDA_SV3AltMTLs/BitweightFiles/sv3/bright/sv3bw-bright-AllTiles.fits'
        #self.brightbitweightfile = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run128/BitweightFiles/sv3/bright/sv3bw-bright-AllTiles.fits'
        self.weightmode = weightmode
        #'/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run64_2/BitweightsRound2/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
        
class main:
    def __init__(self,tp,specver='guadalupe',survey='main'):
        self.mdir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/main/' #location of ledgers
        self.tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/'#location of targets
        ss = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
        md = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/mtl-done-tiles.ecsv')
        self.mtld = join(ss,md,keys=['TILEID','ARCHIVEDATE'])
        self.tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv')
        self.ebits = None
        self.badfib = None
        self.tsnrcol = 'TSNR2_ELG'
        self.tsnrcut = 0
        self.dchi2 = 0
        self.zmin = 0
        self.zmax = 4.5
        if tp[:3] == 'BGS':
            self.imbits = [1,13]
            self.tsnrcut = 1000
            self.tsnrcol = 'TSNR2_BGS'
            self.dchi2 = 40
            self.zmin = 0.1
            self.zmax = 0.5
            self.ebits = [11,12] 
        else:
            self.imbits = [1,12,13]
        if tp[:3] == 'QSO':
            self.ebits = [8,9,11]    
            self.tsnrcut = 0
            self.dchi2 = 0
            self.zmin = 0.8
            self.zmax = 3.5
            #self.tsnrcol = 'TSNR2_QSO'
        if tp[:3] == 'LRG':
            self.ebits = 'lrg_mask'
            self.tsnrcut = 80
            self.dchi2 = 15
            self.zmin = 0.4
            self.zmax = 1.1
        if tp[:3] == 'ELG':
            self.tsnrcut = 80
            self.dchi2 = 0.9
            self.zmin = 0.8
            self.zmax = 1.6
        if tp[:3] == 'ELG':# or tp[:3] == 'BGS':
            self.ebits = [11]    
        if specver == 'everest':
            self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/everest/ELG/main-elg-everest-tiles.fits'
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/everest/QSO/QSO_catalog_MAIN.fits'
        if specver == 'guadalupe':
            #originally '/global/cfs/cdirs/desi/users/raichoor/spectro/guadalupe/main-elg-guadalupe-tiles.fits'
            #self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/main-elg-guadalupe-tiles.fits'
            self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/dark_emlin_catalog.fits'
            #self.qsozf = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/guadalupe/QSO_cat_guadalupe_cumulative.fits'
            #originally self.qsozf = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/guadalupe/QSO_cat_guadalupe_healpix.fits'
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/QSO_cat_guadalupe_healpix.fits'
            self.badfib = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/lrg+bgs_3sig_bad_fibers.txt')
        if specver == 'newQSOtemp_tagged':
            #originally '/global/cfs/cdirs/desi/users/raichoor/spectro/guadalupe/main-elg-guadalupe-tiles.fits'
            #self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/main-elg-guadalupe-tiles.fits'
            self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/newQSOtemp_tagged/dark_emlin_catalog.fits'
            #self.qsozf = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/guadalupe/QSO_cat_guadalupe_cumulative.fits'
            #originally self.qsozf = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/guadalupe/QSO_cat_guadalupe_healpix.fits'
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/QSO_cat_guadalupe_healpix.fits'
            self.badfib = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/lrg+bgs_3sig_bad_fibers.txt')

        if specver == 'daily':
            self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/daily/emlin_catalog.fits'
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/daily/QSO_catalog.fits'
            self.badfib = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/lrg+bgs_3sig_bad_fibers.txt')
