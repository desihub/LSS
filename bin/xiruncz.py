#make sure to type these two commands:
#export OMP_NUM_THREADS=64
#module load gsl
import subprocess
import sys
import argparse
sys.path.append('../py')
import LSS.mkCat_singletile.xitools as xt


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')

args = parser.parse_args()

type = args.type
basedir = args.basedir
version = args.version


lssdir = basedir+'/SV1/LSS/LSScats/'
dirout = svdir+'LSScats/'+version+'/'


if type == 'LRG':
    zmin=0.6
    zmax=1.

xt.prep4czxi(type,zmin,zmax,nran=10,indir=lssdir,ver=version,outdir=os.environ['CSCRATCH']+'/cz/')
subprocess.run(['chmod','+x','czpc.sh'])
subprocess.run('./czpc.sh')
xt.calcxi_dataCZ(type,zmin,zmax,ver=version)
