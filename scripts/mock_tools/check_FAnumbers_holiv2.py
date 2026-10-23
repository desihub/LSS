import fitsio
import os
import numpy as np
mockdir = '/pscratch/sd/d/desica/DA2/mocks/holi_v2/'
outfn = mockdir+'tracer_num.txt'
fo = open(outfn,'w')
fo.write('#realization num_LRG num_ELG num_QSO\n')
mini = 0
maxi = 1001
for realn in range(mini,maxi):
    fname = mockdir+'forFA'+str(realn)+'.fits'
    if os.path.isfile(fname):
        f = fitsio.read(fname,columns=['DESI_TARGET'])
        #print(len(f),f.dtype.names)
        sel_LRG = f['DESI_TARGET'] & 1 > 0
        sel_ELG = f['DESI_TARGET'] & 2 > 0
        sel_QSO = f['DESI_TARGET'] & 4 > 0
        fo.write(str(realn)+' '+str(np.sum(sel_LRG))+' '+str(np.sum(sel_ELG))+' '+str(np.sum(sel_QSO))+'\n')
        print(str(realn)+' '+str(np.sum(sel_LRG))+' '+str(np.sum(sel_ELG))+' '+str(np.sum(sel_QSO)))
fo.close()
