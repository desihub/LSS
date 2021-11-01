from astropy.table import Table

class SV3:
    def __init__(self,tp,weightmode='probobs'):
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
        if tp[:3] == 'LRG' or tp[:3] == 'QSO':
            self.ebits = [8,9,11]    
        if tp[:3] == 'ELG' or tp[:3] == 'BGS':
            self.ebits = [11]    
        self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/ELG/sv3-elg-everest-tiles.fits'
        self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/QSO/QSO_catalog_SV3.fits'
        self.darkbitweightfile = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run128/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
        self.brightbitweightfile = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run128/BitweightFiles/sv3/bright/sv3bw-bright-AllTiles.fits'
        self.weightmode = weightmode
        #'/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_run64_2/BitweightsRound2/BitweightFiles/sv3/dark/sv3bw-dark-AllTiles.fits'
        
class main:
    def __init__(self,tp):
        self.mdir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/main/' #location of ledgers
        self.tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/'#location of targets
        self.mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
        self.tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv')
        self.ebits = None
        
        if tp[:3] == 'BGS':
            self.imbits = [1,13]
        else:
            self.imbits = [1,12,13]
        if tp[:3] == 'LRG' or tp[:3] == 'QSO':
            self.ebits = [8,9,11]    
        if tp[:3] == 'ELG' or tp[:3] == 'BGS':
            self.ebits = [11]    
        self.elgzf = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/everest/ELG/main-elg-everest-tiles.fits'
        self.qsozf = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/everest/QSO/QSO_catalog_MAIN.fits'
