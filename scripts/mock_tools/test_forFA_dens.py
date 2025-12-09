import sys,os
import numpy as np
import fitsio
from desitarget import targetmask

def comp_forFAdens_dark(mockdir,realn=0,dataver='loa-v1/LSScats/v2/',survey='DA2',rootdir='/global/cfs/cdirs/desi/survey/catalogs/'):
    #ranf = rootdir+survey+'/LSS/rands_intiles_DARK_nomask_0.fits'
    ranf = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/randoms/imaging_mask_applied/rands_intiles_DARK_with_imagingmask_0.fits'
    nran = fitsio.read_header(ranf,ext=1)['NAXIS2']
    area = nran/2500

    fadata = fitsio.read(mockdir+'forFA'+str(realn)+'.fits')
    sel_contam = fadata['TARGETID'] > 1e13
    target_types = ['ELG','ELG_LOP','LRG','QSO']
    obsarea_noveto = fitsio.read_header(rootdir+survey+'/LSS/'+dataver+'dark_0_full_noveto.ran.fits',ext=1)['NAXIS2']/2500
    for tp in target_types:
        nq = ''
        if 'ELG' in tp:
            nq = 'notqso'
        obsarea = fitsio.read_header(rootdir+survey+'/LSS/'+dataver+tp+nq+'_0_full_HPmapcut.ran.fits',ext=1)['NAXIS2']/2500
        
        ndat = fitsio.read_header(rootdir+survey+'/LSS/'+dataver+tp+nq+'_full_HPmapcut.dat.fits',ext=1)['NAXIS2']
        seltar = fadata['DESI_TARGET'] & targetmask.desi_mask[tp] > 0
        datdens = ndat/obsarea
        mockdens = np.sum(seltar&~sel_contam)/area
        mockdens_contam = np.sum(seltar&sel_contam)/obsarea_noveto
        mockdenst = mockdens+mockdens_contam
        print('The '+tp+' mock contaminant density is '+str(mockdens_contam)+' and the uncontaminated density is '+str(mockdens))
        print('The '+tp+' data and mock target densities are '+str(datdens)+' and '+str(mockdenst)+'; the ratio is '+str(round(datdens/mockdenst,3)))
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

mockdir = sys.argv[1]
realn = sys.argv[2]

comp_forFAdens_dark(mockdir,realn)
