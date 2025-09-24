
from astropy.table import Table

def fix_dupcol(fn,col='FIBER'):
    of = Table.read(fn)
    try:
        print(len(of[col+'_2']))
    except:
        print('did not find column '+col+'_2 for file '+fn)
        return True
    
    of.remove_columns([col+'_2'])
    of[col+'_1'].name = col
    of.write(fn,overwrite=True,format='fits')
    
indir = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/'


types = ['ELG','ELG_HIP','LRG','QSO','BGS_ANY','BGS_BRIGHT']

nran = 18
for ii in range(0,nran):
    print(str(ii))
    fix_dupcol(indir+'rancomb_'+str(ii)+'brightwdupspec_Alltiles.fits')
    fix_dupcol(indir+'rancomb_'+str(ii)+'darkwdupspec_Alltiles.fits')
    for tp in types:
        print(tp)
        fix_dupcol(indir+'LSScats/test/'+tp+'_'+str(ii)+'_full_noveto.ran.fits')
        fix_dupcol(indir+'LSScats/test/'+tp+'_'+str(ii)+'_full.ran.fits')
        