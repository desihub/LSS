import fitsio
import numpy as np
from desitarget import targetmask
import sys

def comp_forFAdens_dark(mockdir,dataver='loa-v1/LSScats/v2/',survey='DA2',rootdir='/global/cfs/cdirs/desi/survey/catalogs/'):
    ranf = rootdir+survey+'/LSS/rands_intiles_DARK_nomask_0.fits'
    nran = fitsio.read_header(ranf,ext=1)['NAXIS2']
    area = nran/2500

    fadata = fitsio.read(mockdir+'forFA0_noimagingmask_applied_withcontaminants.fits')
    target_types = ['ELG','ELG_LOP','LRG','QSO']
    for tp in target_types:
        nq = ''
        if 'ELG' in tp:
            nq = 'notqso'
        obsarea = fitsio.read_header(rootdir+survey+'/LSS/'+dataver+tp+nq+'_0_full_HPmapcut.ran.fits',ext=1)['NAXIS2']/2500
        ndat = fitsio.read_header(rootdir+survey+'/LSS/'+dataver+tp+nq+'_full_HPmapcut.dat.fits',ext=1)['NAXIS2']
        seltar = fadata['DESI_TARGET'] & targetmask.desi_mask[tp] > 0
        datdens = ndat/obsarea
        mockdens = np.sum(seltar)/area
        print('The '+tp+' data and mock target densities are '+str(datdens)+' and '+str(mockdens)+'; the ratio is '+str(round(datdens/mockdens,3)))
        if 'ELG' in tp:
            selhip = fadata['DESI_TARGET'] & targetmask.desi_mask['ELG_HIP'] > 0
            print('The fraction that are ELG_HIP is '+str(round(np.sum(selhip&seltar)/np.sum(seltar),3))+'; should be 0.1')
        if 'QSO' in tp:
            selhiz = fadata['RSDZ'] > 2.1
            dat = fitsio.read(rootdir+survey+'/LSS/'+dataver+tp+nq+'_full_HPmapcut.dat.fits',columns=['Z_not4clus','ZWARN'])
            selhiz_dat = dat['Z_not4clus'] > 2.1
            selhiz_dat &= dat['Z_not4clus'] != 999999
            selobs = dat['ZWARN'] != 999999
            print('The fractions of data and mock QSO targets that have z > 2.1 are '+str(round(np.sum(selhiz)/np.sum(seltar),3))+','+str(round(np.sum(selhiz_dat)/np.sum(selobs),3)))

comp_forFAdens_dark(sys.argv[1])
