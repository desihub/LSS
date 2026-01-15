import os
import sys

indir = '/global/cfs/cdirs/desi/mocks/cai/LSS/'+sys.argv[1]+'/mocks/'+sys.argv[2]
out3 = open(indir +'/has_QSOLRGELG_LOP.txt','w')
out4 = open(indir +'/has_QSOLRGELGELG_LOP.txt','w')
missing3 = open(indir +'/missing_atleast_oneof_QSOLRGELG_LOP.txt','w')
missinge = open(indir +'/missing_ELG.txt','w')
imin = 1
imax = 1001
tracers = ['QSO','LRG','ELG_LOPnotqso']
for i in range(imin,imax):
    lssdir = indir+'altmtl'+str(i)+'/loa-v1/mock'+str(i)+'/LSScats/'
    h3 = 1
    he = 1
    ml = []
    for tr in tracers:        
        fd = lssdir + tr+'_NGC_clustering.dat.fits'
        if not os.path.isfile(fd):
            ml.append(fd)
            h3 = 0
        fr = lssdir + tr+'_0_NGC_clustering.ran.fits'
        if not os.path.isfile(fr):
            ml.append(fr)
            h3 = 0
    if h3 == 1:
        out3.write(str(i)+'\n')
        me = []
        tr = 'ELGnotqso'
        fd = lssdir + tr+'_NGC_clustering.dat.fits'
        if not os.path.isfile(fd):
            ml.append(fd)
            he = 0
        fr = lssdir + tr+'_0_NGC_clustering.ran.fits'
        if not os.path.isfile(fr):
            ml.append(fr)
            he = 0
        if he == 1:
            out4.write(str(i)+'\n')
        else:
            missinge.write(str(i)+'\n')
    else:
        missing3.write(str(i)+' ')
        for j in range(0,len(ml)):
            missing3.write(ml[j]+' ')
            missing3.write('\n')
  
missing3.close()
out3.close()
missinge.close()
out4.close()
