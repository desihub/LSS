from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as Spline

data = fits.open('cosmos-elg-match.fits')[1].data

plt.ion()
plt.show()

plt.figure()

counts, bins = np.histogram(data['COSMOS2020_ZPHOT'][~data['GOOD_ELG']],range=(0,3),bins=31)

plt.plot(0.5 * (bins[1:] + bins[:-1]), counts, color='b',label='ELG Failure')

counts, bins = np.histogram(data['COSMOS2020_ZPHOT'][data['GOOD_ELG']],range=(0,3),bins=31)

plt.plot(0.5 * (bins[1:] + bins[:-1]), counts, color='r',label='Good ELG')

counts, bins = np.histogram(data['COSMOS2020_ZPHOT'],range=(0,3),bins=31)

plt.plot(0.5 * (bins[1:] + bins[:-1]), counts, color='g',label='All ELG')


#plt.hist(data['COSMOS2020_ZPHOT'][data['GOOD_ELG']],range=(0,3),bins=31,histtype='step',color='r',label='Good ELG')
#plt.hist(data['COSMOS2020_ZPHOT'],range=(0,3),bins=31,histtype='step',color='g',label='All ELG')

plt.legend(fontsize=15)


plt.xlabel('Redshift',fontsize=15)
plt.ylabel('dN/dz',fontsize=15)

plt.savefig('cosmos-elg-dndz.png')

plt.figure()

plt.hist(data['COSMOS2020_ZPHOT'][~data['GOOD_ELG']],range=(0,3),bins=61,histtype='step',color='b',label='ELG Failure')
plt.hist(data['COSMOS2020_ZPHOT'][data['GOOD_ELG']],range=(0,3),bins=61,histtype='step',color='r',label='Good ELG')


counts, bins = np.histogram(data['COSMOS2020_ZPHOT'][~data['GOOD_ELG']],range=(0,3),bins=61)
bc = 0.5 * (bins[1:] + bins[:-1])

spl = Spline(bc, counts, s=0,k=1)

a, b = np.polyfit(bc[(bc > 0.5) & (bc < 1.49)],counts[(bc > 0.5) & (bc < 1.49)],1)

plt.plot(bc, spl(bc),color='k',label='Spline fit')
plt.plot(bc, a * bc + b, color='k', linestyle='--', label='Linear fit')

plt.legend(fontsize=15)

elg_density = 2404 # In per square degrees
failure_rate = 0.241

len_in_range = np.shape(np.where(data['COSMOS2020_ZPHOT'][~data['GOOD_ELG'] & (data['COSMOS2020_ZPHOT'] > 0.5) & (data['COSMOS2020_ZPHOT'] < 1.6)]))[-1]

len_failures = np.shape(np.where(data['COSMOS2020_ZPHOT'][~data['GOOD_ELG']]))[-1]

density_in_range = elg_density * failure_rate * (len_in_range/len_failures)
print(density_in_range)

elg_data = np.loadtxt('ELGnotqso_full_HPmapcut_nz.txt')

plt.figure()

plt.plot(elg_data[:,0],elg_data[:,4]/np.sum(elg_data[:,4] * np.gradient(elg_data[:,0])),label='loa-v1')

def splfunc(x, a, b):
	out = np.zeros_like(x)
	#out = spl(1.49) * np.ones_like(x) - 5 * (x - 1.49)
	out = spl(x)
	out[x < 1.49] = a * x[x < 1.49] + b
	out[x < 0.5] = 0
	out[x > 1.6] = 0
	out[out < 0] = 0
	return out/np.sum(out * np.gradient(x))

plt.plot(elg_data[:,0], splfunc(elg_data[:,0], a, b) * failure_rate * (len_in_range/len_failures) / (1-failure_rate),label='Failures, from COSMOS')
	
#plt.plot(elg_data[:,0],elg_data[:,4]/np.sum(elg_data[:,4] * np.gradient(elg_data[:,0])) + 
#	splfunc(elg_data[:,0], a, b) * failure_rate * (len_in_range/len_failures) / (1-failure_rate))
plt.grid()

total_nz = (elg_data[:,4]/np.sum(elg_data[:,4] * np.gradient(elg_data[:,0])) + 
	splfunc(elg_data[:,0], a, b) * failure_rate * (len_in_range/len_failures) / (1-failure_rate))
	
c, d, = np.polyfit(elg_data[:,0][(elg_data[:,0] > 1.4) & (elg_data[:,0] < 1.5)], 
	total_nz[(elg_data[:,0] > 1.4) & (elg_data[:,0] < 1.5)], 1)
	
new_nz = np.zeros_like(elg_data[:,0])
new_nz[elg_data[:,0] < 1.49] = total_nz[elg_data[:,0] < 1.49]
new_nz[elg_data[:,0] > 1.49] = c * elg_data[:,0][elg_data[:,0] > 1.49] + d
	
plt.plot(elg_data[:,0], new_nz, label='Proposal')

plt.plot(elg_data[:,0], new_nz - (elg_data[:,4]/np.sum(elg_data[:,4] * np.gradient(elg_data[:,0]))), label='Failures, from proposal')

plt.xlabel('Redshift',fontsize=15)
plt.ylabel('dN/dz',fontsize=15)

plt.legend(fontsize=15)

plt.savefig('proposal-dndz.png')

new_Nbin = new_nz * np.sum(elg_data[:,4] * np.gradient(elg_data[:,0]))
new_nbar = new_Nbin / elg_data[:,5]

plt.figure()
plt.plot(elg_data[:,0],elg_data[:,3])
plt.plot(elg_data[:,0], new_nbar)

np.savetxt('elg_nz_including_contaminants.txt',np.array([elg_data[:,0],elg_data[:,3],new_nbar]).T,
	header='zmid n(z) from loa-v1 n(z) with contaminant')
