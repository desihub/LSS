import sys,os
import argparse

#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

import LSS.common_tools as common

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='everest')


args = parser.parse_args()
print(args)

basedir = args.basedir
version = args.version
specrel = args.verspec

ntile = args.ntile
rcut = args.rcut
if rcut is not None:
    rcutstr = rcut.split(',')
    rcut = []
    rcut.append(float(rcutstr[0]))
    rcut.append(float(rcutstr[1]))
ccut = args.ccut

dadir = basedir +'/DA02/LSS/'



ldirspec = dadir+specrel+'/'

dirout = ldirspec+'LSScats/'+version+'/'


types = ['ELG','ELG_LOP','LRG','ELG_LOPnoqso','QSO','BGS_ANY','BGS_BRIGHT']

for type in types:
    wzm = ''
#     if zmask:
#         wzm = 'zmask_'
    if rcut is not None:
        wzm += '_rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntile > 0:
        wzm += '_ntileg'+str(ntilecut)+'_'    
    if ccut is not None:
        wzm += '_'+ccut #you could change this to however you want the file names to turn out

    regl = ['_DN','_DS','','_N','_S']
    
    for reg in regl:
        fb = dirout+type+wzm+reg
        fcr = fb+'_0_clustering.ran.fits'
        fcd = fb+'_clustering.dat.fits'
        fout = fb+'_nz.dat'
        if type == 'QSO':
            zmin = 0.6
            zmax = 4.5
            dz = 0.05
            
        else:    
            dz = 0.02
            zmin = 0.01
            zmax = 1.61
        common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
        common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax)
