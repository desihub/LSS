from astropy.table import Table

class SV3:
    def __init__(self,tp,weightmode='probobs',specver='fuji'):
        self.mdir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/sv3/' #location of ledgers
        self.tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/targets/sv3/resolve/'#location of targets
        self.mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
        self.tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')
        ebits = None
        if tp[:3] == 'BGS':
            self.imbits = [1,13]
        else:
            self.imbits = [1,12,13]
        self.ebits = None
        if tp[:3] == 'QSO':
            self.ebits = [8,9,11]    
        if tp[:3] == 'LRG':
            self.ebits = 'lrg_mask'
        if tp[:3] == 'ELG' or tp[:3] == 'BGS':
            self.ebits = [11]    
        if specver == 'everest':
            self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/ELG/sv3-elg-everest-tiles.fits'
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/QSO/QSO_catalog_SV3.fits'
        if specver == 'fuji':
            self.elgzf = '/global/cfs/cdirs/desi/users/raichoor/spectro/fuji/sv3-elg-fuji-tiles.fits'
            self.qsozf = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/fuji/QSO_cat_fuji_cumulative.fits'
        
        self.darkbitweightfile = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run128/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
        self.brightbitweightfile = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run128/BitweightFiles/sv3/bright/sv3bw-bright-AllTiles.fits'
        self.weightmode = weightmode
        #'/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run64_2/BitweightsRound2/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
        
class main:
    def __init__(self,tp,specver='guadalupe'):
        self.mdir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/main/' #location of ledgers
        self.tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/'#location of targets
        self.mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
        self.tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv')
        self.ebits = None
        
        if tp[:3] == 'BGS':
            self.imbits = [1,13]
        else:
            self.imbits = [1,12,13]
        if tp[:3] == 'QSO':
            self.ebits = [8,9,11]    
        if tp[:3] == 'LRG':
            self.ebits = 'lrg_mask'
        if tp[:3] == 'ELG' or tp[:3] == 'BGS':
            self.ebits = [11]    
        if specver == 'everest':
            self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/everest/ELG/main-elg-everest-tiles.fits'
            self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/everest/QSO/QSO_catalog_MAIN.fits'
        if specver == 'guadalupe':
            self.elgzf = '/global/cfs/cdirs/desi/users/raichoor/spectro/guadalupe/main-elg-guadalupe-tiles.fits'
            #self.qsozf = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/guadalupe/QSO_cat_guadalupe_cumulative.fits'
            self.qsozf = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/guadalupe/QSO_cat_guadalupe_healpix.fits'
        
        #recon parameters
        self.om = 0.31519
        if tp[:3] == 'LRG':
            self.bias = 1.8
            self.ff = 0.4*self.bias
