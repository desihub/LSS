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


dadir = basedir +'/main/LSS/'



ldirspec = dadir+specrel+'/'

dirout = ldirspec+'LSScats/'+version+'/'


types = ['ELG','ELG_LOP','LRG','ELG_LOPnoqso','QSO','BGS_ANY','BGS_BRIGHT']

for type in types:
    wzm = 'zdone'

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
