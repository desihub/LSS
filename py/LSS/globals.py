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
    def __init__(self,tp,specver='iron',survey='main',relax_zbounds='n'):
        print(tp)
        self.mdir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/main/' #location of ledgers
        
        self.tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/'#location of targets
        ss = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
        md = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/mtl-done-tiles.ecsv')
        self.mtld = join(ss,md,keys=['TILEID','ARCHIVEDATE'])
        self.tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv')
        self.ebits = None
        self.badfib = None
        self.badfib_td = None
        self.badfib_status = None
        self.tsnrcol = 'TSNR2_ELG' 
        #self.tsnrcut = 0 #better to throw an error having this undefined rather than having code use a bad value
        #self.dchi2 = 0
        self.zmin = 0
        self.zmax = 4.5
        zfloor = 0.002
        self.reccircmasks=None
        self.mapcuts = {'EBV':0.15,'STARDENS':4.4,'PSFSIZE_G':2.4,'PSFSIZE_R':2.3,'PSFSIZE_Z':2,'GALDEPTH_G':250,'GALDEPTH_R':80,'GALDEPTH_Z':30,'PSFDEPTH_W1':2}
        if tp[:3] == 'BGS':
            self.fit_maps_all = ['STARDENS','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','HI'] #used up until v0.6
            self.fit_maps = ['STARDENS','GALDEPTH_R','HI']
            self.fit_maps_allebv = ['STARDENS','GALDEPTH_R','HI','EBV_DIFF_GR','EBV_DIFF_RZ']
            self.fit_maps_allebvcmb = ['STARDENS','GALDEPTH_R','HI','EBV_DIFF_GR','EBV_DIFF_RZ','ZCMB']

            self.tsnrcut = 1000
            self.tsnrcol = 'TSNR2_BGS'
            self.dchi2 = 40
            self.zmin = zfloor
            self.zmax = 0.6
            if tp == 'BGS_BRIGHT-21.5':
                self.zmin = 0.1
                self.zmax = 0.4
            #self.zmin = 0.1
            #self.zmax = 0.5
            #if survey == 'Y1':
            #    self.zmax = 0.4
            self.ebits = [11] 
            self.imbits = [1,13]
        else:
            self.imbits = [1,12,13]
        if tp[:3] == 'QSO':
            self.fit_maps = ['PSFDEPTH_W1','PSFDEPTH_W2','STARDENS','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','EBV_DIFF_GR','EBV_DIFF_RZ','HI']
            self.fit_maps_all = self.fit_maps
            self.ebits = [8,9,11]    
            self.tsnrcut = 80
            self.dchi2 = 0
            self.zmin = 0.8
            if relax_zbounds == 'y':
                self.zmin = zfloor
            self.zmax = 3.5
            self.reccircmasks=['/global/cfs/cdirs/desi/users/rongpu/desi_mask/desi_custom_mask_v1.txt']
            #self.tsnrcol = 'TSNR2_QSO'
        if tp[:3] == 'LRG':
            self.fit_maps_all = ['STARDENS','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','HI','PSFDEPTH_W1'] #used up until v0.6
            self.fit_maps_allebv = ['STARDENS','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','HI','PSFDEPTH_W1','EBV_DIFF_GR','EBV_DIFF_RZ']
            self.fit_maps_allebvcmb = ['STARDENS','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','HI','PSFDEPTH_W1','EBV_DIFF_GR','EBV_DIFF_RZ','ZCMB']
            self.fit_maps = ['STARDENS','PSFSIZE_R','GALDEPTH_Z','HI','PSFDEPTH_W1']
            self.fit_maps46s = ['STARDENS','PSFSIZE_R','GALDEPTH_Z','HI','PSFDEPTH_W1','GALDEPTH_R']
            self.fit_maps68s = ['STARDENS','PSFSIZE_R','GALDEPTH_Z','HI','PSFDEPTH_W1','GALDEPTH_G']
            self.fit_maps81s = ['STARDENS','PSFSIZE_R','GALDEPTH_Z','HI','PSFDEPTH_W1','PSFSIZE_Z']
            self.ebits = 'lrg_mask'
            self.tsnrcut = 80
            self.dchi2 = 15
            self.zmin = 0.4
            self.zmax = 1.1
            if relax_zbounds == 'y':
                self.zmin = zfloor
                self.zmax = 1.5
            
        if tp[:3] == 'ELG':
            self.fit_maps = ['STARDENS','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','EBV_DIFF_GR','EBV_DIFF_RZ','HI']#,'EBV_DIFF_MPF']
            self.fit_maps_all = self.fit_maps
            self.tsnrcut = 80
            self.dchi2 = 0.9
            self.zmin = 0.8
            self.zmax = 1.6
            if relax_zbounds == 'y':
                self.zmin = zfloor

            self.reccircmasks=['/global/cfs/cdirs/desi/users/rongpu/desi_mask/desi_custom_mask_v1.txt','/global/cfs/cdirs/desi/users/rongpu/desi_mask/elg_custom_mask_v1.1_draft.txt']

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
            self.badfib = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/unique_badfibers.txt')
            #self.badfib = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/lrg+bgs_3sig_bad_fibers.txt')
        
        if specver == 'himalayas':
            self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/'+specver+'/emlin_catalog.fits'
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/'+specver+'/QSO_catalog.fits'
            self.badfib = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/lrg+bgs_3sig_bad_fibers.txt')

        if specver == 'iron':
            self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/'+specver+'/emlin_catalog.fits'
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/Y1/QSO/'+specver+'/QSO_cat_iron_cumulative_v0.fits'
            self.badfib = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/unique_badfibers.txt')

        if specver == 'jura-v1':
            self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/'+specver+'/emlin_catalog.fits'
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/DA2/QSO/jura/QSO_cat_jura_cumulative_v1.fits'
            self.badfib = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/unique_badfibers.txt')

        if specver == 'kibo-v1' or specver == 'loa-v1':
            self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/'+specver+'/emlin_catalog.fits'
            
            
            self.badfib = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/kibo-v1/unique_badfibers.txt')
            if specver == 'loa-v1':
                self.badfib_td = open('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/unique_badfibers_time-dependent.txt').readlines()
            self.badfib_status  = [13,14]
        if specver == 'kibo-v1':
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs//DA2/QSO/kibo/QSO_cat_kibo_cumulative_v1.fits'
        if specver == 'loa-v1':
            #self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs//DA2/QSO/loa/QSO_cat_loa_cumulative_v0.fits' #used for v1
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs//DA2/QSO/loa/QSO_cat_loa_cumulative_v2.fits' #used for v1.1 onward


        #self.darkbitweightfile = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/mainbw-dark-allTiles_v1.fits'
        #self.brightbitweightfile = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/mainbw-bright-allTiles_v1.fits'
        
        #properties for maps
        self.new_cols=['EBV_CHIANG_SFDcorr','STARDENS','HALPHA', 'HALPHA_ERROR', 'CALIB_G', 'CALIB_R', 'CALIB_Z', 'EBV_MPF_Mean_FW15', 'EBV_MPF_Mean_ZptCorr_FW15', 'EBV_MPF_Var_FW15', 'EBV_MPF_VarCorr_FW15', 'EBV_MPF_Mean_FW6P1', 'EBV_MPF_Mean_ZptCorr_FW6P1', 'EBV_MPF_Var_FW6P1', 'EBV_MPF_VarCorr_FW6P1', 'EBV_SGF14', 'BETA_ML', 'BETA_MEAN', 'BETA_RMS', 'HI', 'KAPPA_PLANCK']
        self.fid_cols=['EBV','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
