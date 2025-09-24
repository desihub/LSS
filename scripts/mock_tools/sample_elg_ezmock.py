import numpy as np
import time

elg_frac_fn = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/y1_elg_lop_fraction.txt'

def fraction_func(x):
    start_time = time.time()
    red,bin_values = np.loadtxt(elg_frac_fn, unpack=True)
    red_inf = red-0.005
    bin_limits = np.append(red_inf,1.6)


# Convert inputs to NumPy arrays for efficient calculations
    x_array = np.array(x)
    bin_limits = np.array(bin_limits)
    bin_values = np.array(bin_values)

    # Find the indices of the bins that contain the given values of x
    bin_indices = np.digitize(x_array, bin_limits, right=True) - 1

    print("--- fraction function %s seconds ---" % (time.time() - start_time))
    return bin_values[bin_indices]


def create_subsample(main_sample):
    x_values = main_sample['Z']
    fractions = fraction_func(x_values)
    start_time = time.time()
    mask = np.random.rand(len(main_sample)) < fractions
    df = main_sample.to_pandas()
    subsample_lop = df.loc[mask]
    subsample_vlo = df.loc[~mask]
    print("--- create_subsample %s seconds ---" % (time.time() - start_time))
    print(len(subsample_lop), len(subsample_vlo))
    return subsample_lop, subsample_vlo


'''
# Create the subsample
lop, vlo = create_subsample(da, fraction_func)

print("LOP Sample:", len(lop))
print("VLO sample:", len(vlo))
import matplotlib.pyplot as plt
n_main,bins,_ = plt.hist(da['Z'],bins=100)
n_lop,bins,_ = plt.hist(lop['Z'],bins=bins)
n_vlo,bins,_ = plt.hist(vlo['Z'],bins=bins)
xcenters = (bins[:-1] + bins[1:]) / 2
np.savetxt('nz_main_lop_vlo.txt', np.array([xcenters, n_main, n_lop, n_vlo]).T)
from astropy.table import Table

t_main = Table(da)
t_lop = Table(lop)
t_vlo = Table(vlo)
t_main.write('main.fits')
t_vlo.write('vlo.fits')
t_lop.write('lop.fits')
'''
