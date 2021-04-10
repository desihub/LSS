import numpy as np
import fitsio

ran = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random0/rancomb_dark_Alltiles.fits')

print('#area covered on DARK tiles\n#>N_tiles area(deg2)')
for nt in np.unique(ran['NTILE']):
    wt = ran['NTILE'] >= nt
    print(nt,len(ran[wt])/2500)
    
print('\n')
ran = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random0/rancomb_bright_Alltiles.fits')

print('#area covered on BRIGHT tiles\n#>N_tiles area(deg2)')
for nt in np.unique(ran['NTILE']):
    wt = ran['NTILE'] >= nt
    print(nt,len(ran[wt])/2500)
    
print('\n')
dat = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/LSScats/test/ALLAlltiles_dark_full.dat.fits')
wz = dat['ZWARN']*0 == 0
wz &= dat['ZWARN'] != 999999
wzg = dat['ZWARN'] == 0
print('number of dark time observations (good hardware):',len(dat[wz]))
print('number of dark time good observations (ZWARN==0):',len(dat[wzg]))

print('\n')
dat = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/LSScats/test/ALLAlltiles_bright_full.dat.fits')
wz = dat['ZWARN']*0 == 0
wz &= dat['ZWARN'] != 999999
wzg = dat['ZWARN'] == 0
print('number of bright time observations (good hardware):',len(dat[wz]))
print('number of bright time good observations (ZWARN==0):',len(dat[wzg]))

print('\n')
print('splitting by type, numbers are after all veto masks:')

tpl = ['ELG','LRG','QSO','BGS_ANY']
for tp in tpl:
    dat = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/LSScats/test/'+tp+'Alltiles_full.dat.fits')
	wz = dat['ZWARN']*0 == 0
	wz &= dat['ZWARN'] != 999999
	wzg = dat['ZWARN'] == 0
	print('number of '+tp+' observations (good hardware):',len(dat[wz]))
	print('number of '+tp+' time good observations (ZWARN==0):',len(dat[wzg]))
	ran = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/LSScats/test/'+tp+'Alltiles_clustering.dat.fits')
	print('effective '+tp+' area, after vetoing higher-priority positions, positions with redshift failures, and imaging: ',str(len(ran)/2500))
	print('\n')
	

    
    
