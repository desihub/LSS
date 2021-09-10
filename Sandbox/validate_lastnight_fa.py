#right now, requires source /project/projectdirs/desi/software/desi_environment.sh master
#and export PYTHONPATH=/global/homes/a/ajross/desicode/fiberassign_rerun/py:$PYTHONPATH
#the 2nd likley only works for AJR at the moment, once fba_rerun_io is in master branch 2nd shouldn't be necessary
from astropy.table import Table
import os

#open exposures file
exps = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/exposures.ecsv')
#remove the rows that don't have night so that next lines work
sel = exps['NIGHT'] != 'None'
exps = exps[sel]
#find the most recent night
maxn = np.max(exps['NIGHT'].astype(int))
seln = exps['NIGHT'].astype(int) == maxn
#get the list of tileids observed on the last night
tidl = np.unique(exps[seln]['TILEID'])

print('looking at fiberassign files used on the night '+str(maxn))
print('the tileids are:')
print(tidl)

if not os.path.exists(os.environ['CSCRATCH']+'/fiberassign'):
    os.mkdir(os.environ['CSCRATCH']+'/fiberassign')
    print('made '+os.environ['CSCRATCH']+'/fiberassign')

outdir = os.environ['CSCRATCH']+'/fiberassign/rerun'
if not os.path.exists(outdir):
    os.mkdir(outdir)
    print('made '+outdir)

fol = [] #empty list to contain orginal fiberassign file names
fnl = [] #empty list to contain orginal fiberassign file names

for tid in tidl:
    sel = exps['TILEID'] == tid
    sel &= exps['NIGHT'].astype(int) == maxn
    expid = exps[sel]['EXPID'][0] #just select one expid to get the fiberassign file
    fn = '/global/cfs/cdirs/desi/spectro/data/'+str(maxn)+'/'+str(expid).zfill(8)+'/fiberassign-'+str(tid).zfill(6)+'.fits.gz'
    print('reproducing data for fiberassign file '+fn)
    fol.append(fn)
    fnl.append(outdir+'/'+fn[12:15]+'/fiberassign-'+str(tid).zfill(6)+'.fits.gz')
    #system call run fiberassign
    os.system('fba_rerun --infiberassign '+fn+' --outdir '+outdir+' --dtver 1.1.1 --nosteps qa')

try:
    #edit out the fiberassign_rerun piece to be fiberassign once fba_rerun_io is merged into master fiberassign 
    from fiberassign_rerun.fba_rerun_io import fba_rerun_check
    docheck = True
except:
    print('import failed, not doing check')

if docheck:
    for ii in range(0,len(fnl)):
        fba_rerun_check(fol[ii], fnl[ii], str(tidl[ii])+'.diff')   
