import fitsio
from time import time
t0 = time()
ranfl = [] #making a list of file names and then reading all at once, takes forever otherwise for whatever reason
tpstr = 'QSO'
use_map_veto = '_HPmapcut'
dirout = '/global/cfs/cdirs/desi/survey/catalogs//DA2/LSS/loa-v1/LSScats/v2'
for i in range(0,5):
    ranf = os.path.join(dirout, tpstr+'_'+str(i)+'_full'+use_map_veto+'.ran.fits'.replace('global','dvs_ro'))
    ranfl.append(ranf)
ranl = [fitsio.read(ranfi, columns=['RA', 'DEC','PHOTSYS']) for ranfi in ranfl]
print('took '+str(time()-t0)+' seconds')
