import glob
import argparse

sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.SV3.cattools as ct

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='daily')

sv3dir = basedir +'/SV3/LSS/'



ldirspec = sv3dir+specrel+'/'

dirout = ldirspec+'LSScats/'+version+'/'


types = ['ELG','ELG_HIP','LRG','QSO','BGS_ANY','BGS_BRIGHT']

for type in types:
	regl = ['','_N','_S']
	for reg in regl:
		if zma:
			reg = '_zmask'+reg
		fcr = dirout+type+'Alltiles'+reg+'_0_clustering.ran.fits'
		fcd = dirout+type+'Alltiles'+reg+'_clustering.dat.fits'
		fout = dirout+type+reg+'_nz.dat'
		if type == 'QSO':
			zmin = 0.6
			zmax = 4.5
			dz = 0.05
			ct.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax)
		else:    
			ct.mknz(fcd,fcr,fout,bs=0.02)
