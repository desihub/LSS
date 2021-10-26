
from astropy.table import Table

def fix_dupcol(fn,col='FIBER'):
    of = Table.read(fn)
    of.remove_columns([col+'_2'])
    of[col+'_1'].name = col
    of.write(fn,overwrite=True,format='fits')
    
indir = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/'

#fix_dupcol(indir+'datcomb_bright_tarwdup_Alltiles.fits')
fix_dupcol(indir+'datcomb_dark_tarwdupspec_Alltiles.fits')

types = ['ELG','ELG_HIP','LRG','LRG_main','QSO','BGS_ANY','BGS_BRIGHT']

nran = 18
for tp in types:
    print(tp)
    for ii in range(0,nran):
        fix_dupcol(indir+'LSScats/test/'+tp+'_'+str(ii)+'_full.ran.fits')
        print(str(ii))