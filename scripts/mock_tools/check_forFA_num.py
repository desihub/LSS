import fitsio
import numpy as np
import argparse
import os

ref_real = 0
real_min = 1
real_max = 200
mock = 'holi_v3'
indir = os.getenv('SCRATCH')+'/DA2/mocks/'+mock+'/'

ref_dat = fitsio.read(indir+'forFA'+str(ref_real)+'.fits',
                      columns=['TARGETID', 'DESI_TARGET'])
sel_lrg = (ref_dat['DESI_TARGET'] & 1) != 0
nlrg_ref = np.sum(sel_lrg)
sel_qso = (ref_dat['DESI_TARGET'] & 4) != 0
nqso_ref = np.sum(sel_qso)
sel_elg = (ref_dat['DESI_TARGET'] & 2) != 0
nelg_ref = np.sum(sel_elg)
sel_contam = ref_dat['TARGETID'] < 419430400000000
nqso_ref_nocontam = np.sum(sel_qso & sel_contam)
nelg_ref_nocontam = np.sum(sel_elg & sel_contam)
del ref_dat
for real in range(real_min, real_max):
    dat = fitsio.read(indir+'forFA'+str(real)+'.fits',
                      columns=['TARGETID', 'DESI_TARGET'])
    sel_lrg = (dat['DESI_TARGET'] & 1) != 0
    nlrg = np.sum(sel_lrg)
    sel_qso = (dat['DESI_TARGET'] & 4) != 0
    nqso = np.sum(sel_qso)
    sel_elg = (dat['DESI_TARGET'] & 2) != 0
    nelg = np.sum(sel_elg)
    sel_contam = dat['TARGETID'] < 419430400000000
    nqso_nocontam = np.sum(sel_qso & sel_contam)
    nelg_nocontam = np.sum(sel_elg & sel_contam)
    print('realization', real, 'lrg', nlrg, 'qso', nqso, 'elg', nelg, 'qso_nocontam', nqso_nocontam, 'elg_nocontam', nelg_nocontam,
          'lrg_ratio', nlrg/nlrg_ref, 'qso_ratio', nqso /
          nqso_ref, 'elg_ratio', nelg/nelg_ref,
          'qso_nocontam_ratio', nqso_nocontam/nqso_ref_nocontam, 'elg_nocontam_ratio', nelg_nocontam/nelg_ref_nocontam)
