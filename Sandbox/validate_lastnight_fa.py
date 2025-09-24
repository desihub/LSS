#right now, requires source /project/projectdirs/desi/software/desi_environment.sh master
#and export PYTHONPATH=/global/homes/a/ajross/desicode/fiberassign_rerun/py:$PYTHONPATH
#the 2nd likley only works for AJR at the moment, once fba_rerun_io is in master branch 2nd shouldn't be necessary
from astropy.table import Table
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--night", help="use this if you want to specify the night, rather than just use the last one",default=None)
parser.add_argument("--clean", help="remove all files in the output directory prior to start, True/Fase",default=False)
args = parser.parse_args()

#open exposures file
exps = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/exposures.ecsv')
#remove the rows that don't have night so that next lines work
sel = exps['NIGHT'] != 'None'
#check if tileid are in main
tlm = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv')
sel &= np.isin(exps['TILEID'],tlm['TILEID'])
exps = exps[sel]

if args.night == None:
    #find the most recent night
    maxn = np.max(exps['NIGHT'].astype(int))
else:
    maxn = int(args.night)  
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

if args.clean:
    os.system('rm -R '+outdir+'/*')

fol = [] #empty list to contain orginal fiberassign file names
fnl = [] #empty list to contain new fiberassign file names

for tid in tidl:
    sel = exps['TILEID'] == tid
    sel &= exps['NIGHT'].astype(int) == maxn
    expid = exps[sel]['EXPID'][0] #just select one expid to get the fiberassign file
    ff = 'fiberassign-'+str(tid).zfill(6)+'.fits.gz'
    fn = '/global/cfs/cdirs/desi/spectro/data/'+str(maxn)+'/'+str(expid).zfill(8)+'/'+ff
    print('reproducing data for fiberassign file '+fn)
    fol.append(fn)
    fnl.append(outdir+'/'+str(tid).zfill(6)[:3]+'/fiberassign-'+str(tid).zfill(6)+'.fits.gz')
    #system call run fiberassign
    os.system('fba_rerun --infiberassign '+fn+' --outdir '+outdir+' --nosteps qa') #--dtver 1.1.1 

#try:
from fiberassign.fba_rerun_io import fba_rerun_check
docheck = True
#except:
#    print('import failed, not doing check')

tids_passl = []
if docheck:
    for ii in range(0,len(fnl)):
        dfn = outdir+'/'+str(tidl[ii])+'.diff'
        if os.path.isfile(fnl[ii]):
            fba_rerun_check(fol[ii], fnl[ii],dfn )  
            tids = np.genfromtxt(dfn,usecols = (3))
            if len(tids) > 0:
                #tids = dd[3]
                sel = tids > 0
                if len(tids[sel]) > 0:
                    print('found '+str(len(tids[sel]))+' positive targetid that are different')
                    print('FOLLOW-UP NEEDED, DO NOT ALLOW ZDONE FOR TILEID '+str(tidl[ii])+'!!!')
                else:
                    print('TILEID '+str(tidl[ii])+' PASSED') 
                    tids_passl.append(tidl[ii])
            else:
                print('TILEID '+str(tidl[ii])+' PASSED')          
                tids_passl.append(tidl[ii])
        else:
            print('WHY IS THERE NO NEW FIBERASSIGN FILE FOR '+str(tidl[ii])+'!?!? (check above output for clues)')

#move intermediate files for tiles that pass

intermediate_dir = '/global/cfs/cdirs/desi/survey/fiberassign/main/test' #remove the test once we are happy
for tid in tids_passl:
    mv_tiddir = os.path.join(intermediate_dir, str(tid).zfill(6)[:3])
    if not os.path.isdir(mv_tiddir):
        print("create {}".format(mv_tiddir))
        os.mkdir(mv_tiddir)
    for name in ["tiles", "sky", "gfa", "targ", "scnd", "too"]:        
        #fn = os.path.join(outdir, str(tid).zfill(6)[:3], "{:06d}-{}.fits".format(tid))
        fn = outdir+'/'+str(tid).zfill(6)[:3]+'/'+str(tid).zfill(6)+'-'+name+'.fits'
        if os.path.isfile(fn):
            mv_fn =os.path.join(mv_tiddir, os.path.basename(fn))
            print("moving {} to {}".format(fn, mv_fn))
            os.system("mv {} {}".format(fn, mv_fn))

