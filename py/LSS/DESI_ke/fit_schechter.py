import  os
import  argparse
import  numpy as np
import  getdist.plots as gdplt

from    schechter import schechter
from    astropy.table import Table
from    cobaya.run import run
from    scipy import stats
from    getdist.mcsamples import MCSamplesFromCobaya

# Run me on interactive:
# desienv master
# srun -N 1 -n 1 python3 fit_schechter.py

parser  = argparse.ArgumentParser(description='MCMC fitting of measured Schechter function.')
parser.add_argument('--known',  action='store_true', help='Run for toy luminosity function.')

args  = parser.parse_args()
known = args.known

root  = os.environ['GOLD_DIR']
fpath = root + '/gama_gold_lumfn.fits'
lumfn = Table.read(fpath)

if known:
    lumfn['PHI_N'] = schechter(lumfn['MEDIAN_M'], -2.01, -20.89, -1.25)
    lumfn['PHI_N_ERROR'] = 1.e-2 * lumfn['PHI_N']

else:
    # Bins with no objects currently results in Nans. 
    lumfn = lumfn[lumfn['PHI_N'] > 0.]


lumfn = lumfn[lumfn['PHI_N'] > 0.]
 
lumfn.pprint()

def chi2(log10phistar, Mstar, alpha):
    phistar = 10.**log10phistar
    
    res  = (lumfn['PHI_IVMAX'] - schechter(lumfn['MEDIAN_M'], phistar, Mstar, alpha))
    res /= lumfn['PHI_IVMAX_ERROR']
    res  = res * res

    return np.sum(res)

def lnlike(log10phistar, Mstar, alpha):
    return -chi2(log10phistar, Mstar, alpha) / 2.

x2 = chi2(-2.01, -20.89, -1.25)

print('Initial Chi sq. evaluation: {:.4f}'.format(x2))


info = {"likelihood": {"schechter": lnlike}}

# https://cobaya.readthedocs.io/en/latest/params_prior.html
info["params"] = {
    "log10phistar": {"prior": {"dist": "norm", "loc": -2.00, "scale": 0.25},   "ref": -2.00,   "proposal": 0.01},  # TMR ref.  -2.01
     "Mstar":        {"prior": {"dist": "norm", "loc": -20.89, "scale": 0.15}, "ref": -20.89, "proposal": 0.01},  # TMR ref.  -20.89 
     "alpha":        {"prior": {"dist": "norm", "loc": -1.25, "scale": 0.05},  "ref": -1.25,  "proposal": 0.01}}  # TMR ref.  -1.25

# https://cobaya.readthedocs.io/en/latest/sampler_mcmc.html
info["sampler"] = {"mcmc": {"Rminus1_stop": 0.001, "max_tries": 5000}}

updated_info, sampler = run(info, output='{}/cobaya/schechter_chain'.format(root))

print('Written to {}/cobaya/schechter_chain*'.format(root))
