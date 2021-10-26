
from astropy.table import Table

def fix_nullcol(fn,collist=[]):
    of = Table.read(fn)
    try:    
        of.remove_columns(collist)
        of.write(fn,overwrite=True,format='fits')
    except:
        print('columns for removal not in '+fn)    
   
indir = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'


types = ['ELG','ELG_HIP','LRG','QSO','BGS_ANY','BGS_BRIGHT']


colstar = ['REF_EPOCH','PARALLAX','PMRA','PMDEC','OBSCONDITIONS','NUMOBS_INIT',\
'SV3_SCND_TARGET','NUMOBS_MORE','NUMOBS','Z','ZTILEID','VERSION']

fix_nullcol(indir+'datcomb_dark_tarwdup_Alltiles.fits',colstar)
fix_nullcol(indir+'datcomb_bright_tarwdup_Alltiles.fits',colstar)

colspec = ['REF_EPOCH',
 'PARALLAX',
 'PMRA',
 'PMDEC',
 'OBSCONDITIONS',
 'NUMOBS_INIT',
 'NUMOBS_MORE',
 'NUMOBS',
 'ZTILEID',
 'VERSION']
 
fix_nullcol(indir+'everest/datcomb_dark_tarspecwdup_Alltiles.fits',colspec)
fix_nullcol(indir+'everest/datcomb_bright_tarspecwdup_Alltiles.fits',colspec)
 
for tp in types:
    print(tp)
    fix_nullcol(indir+'everest/LSScats/test/'+tp+'_full_noveto.dat.fits',colspec)
    fix_nullcol(indir+'everest/LSScats/test/'+tp+'_full.dat.fits',colspec)
        