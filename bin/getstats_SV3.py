import numpy as np
import fitsio
from astropy.table import Table,unique
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')

args = parser.parse_args()
print(args)

version = args.version


mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/mtl-done-tiles.ecsv') #log of tiles completed for mtl
tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')

wp = tiles['PROGRAM'] == 'DARK'
tiles = tiles[wp]

wp = np.isin(mtld['TILEID'],tiles['TILEID']) #we want to consider MTL done tiles that correspond to the SV3 tile file
mtld = mtld[wp]
print('number of completed dark tiles:',len(mtld))

mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/mtl-done-tiles.ecsv') #log of tiles completed for mtl
tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')

wp = tiles['PROGRAM'] == 'BRIGHT'
tiles = tiles[wp]

wp = np.isin(mtld['TILEID'],tiles['TILEID']) #we want to consider MTL done tiles that correspond to the SV3 tile file
mtld = mtld[wp]
print('number of completed bright tiles:',len(mtld))


#ran = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random0/rancomb_dark_Alltiles.fits')
ran = Table.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random0/rancomb_dark_Alltilelocinfo.fits')
ran = unique(ran,keys=['TARGETID'])

print('#area covered on DARK tiles\n#>N_tiles area(deg2)')
for nt in np.unique(ran['NTILE']):
    wt = ran['NTILE'] >= nt
    print(nt,len(ran[wt])/2500)
    
print('#')
#ran = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random0/rancomb_bright_Alltiles.fits')
ran = Table.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random0/rancomb_bright_Alltilelocinfo.fits')
ran = unique(ran,keys=['TARGETID'])


print('#area covered on BRIGHT tiles\n#>N_tiles area(deg2)')
for nt in np.unique(ran['NTILE']):
    wt = ran['NTILE'] >= nt
    print(nt,len(ran[wt])/2500)
    
print('#')
dat = Table.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/datcomb_dark_tarspecwdup_Alltiles.fits')
wz = dat['ZWARN']*0 == 0
wz &= dat['ZWARN'] != 999999
datz = unique(dat[wz],keys=['TARGETID'])
wzg = dat['ZWARN'] == 0
datzg = unique(dat[wzg],keys=['TARGETID'])
print('number of unique dark time targets observed (good hardware):',len(datz))
print('number of unique dark time targets with a good observation (ZWARN==0):',len(datzg))

print('#')
dat = Table.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/datcomb_bright_tarspecwdup_Alltiles.fits')
wz = dat['ZWARN']*0 == 0
wz &= dat['ZWARN'] != 999999
datz = unique(dat[wz],keys=['TARGETID'])
wzg = dat['ZWARN'] == 0
datzg = unique(dat[wzg],keys=['TARGETID'])

print('number of unique bright time targets observed (good hardware):',len(datz))
print('number of unique bright time targets with a good observation (ZWARN==0):',len(datzg))

print('#')
print('splitting by type, numbers are after all veto masks:')

tpl = ['ELG','LRG','QSO','BGS_ANY','MWS_ANY']
for tp in tpl:
    dat = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/LSScats/'+version+'/'+tp+'Alltiles_full.dat.fits')
    wz = dat['ZWARN']*0 == 0
    wz &= dat['ZWARN'] != 999999
    wzg = dat['ZWARN'] == 0
    print('number of unique '+tp+' targets observed (good hardware):',len(dat[wz]))
    print('number of unique '+tp+' targets with a good observation (ZWARN==0):',len(dat[wzg]))
    ran = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/LSScats/'+version+'/'+tp+'Alltiles_0_clustering.ran.fits')
    print('effective '+tp+' area, after vetoing higher-priority positions, and imaging: ',str(len(ran)/2500))
    print('#')
    

    
    
