from astropy.table import Table
from scipy import special

tracers = ['LRG','QSO','ELG_LOPnotqso','BGS_BRIGHT']

for tp in tracers:
    #data1 = np.loadtxt('/global/cfs/cdirs/desi/users/akrolew/summaryscale/'+tp+'fibersim.txt.3',usecols=range(10))
    data2 = np.loadtxt('/global/cfs/cdirs/desi/users/akrolew/summaryscale/'+tp+'fibersim_kibo-v1.txt',usecols=range(10))
    #data3 = np.loadtxt('/global/cfs/cdirs/desicollab/users/akrolew/summary/ELG_LOPnotqso/0.txt',usecols=range(10)).reshape((1,10))
    data = data2
    
    col_names = ['FIBER', 'p_value', 'var', 'nsig', 'n', 's', 'mean_mod_suc', 'frac_suc', 'mean_x','mean_y']
    exec(tp+'summary = Table(data, names=col_names)')
    exec(tp+'summary ["n_sig_sim"] = np.around((-np.sqrt(2) * special.erfcinv(2 * '+tp+'summary["p_value"])),5)')
    
#np.savetxt('/global/cfs/cdirs/desi/users/akrolew/summaryscale/'+tp+'fibersim.txt',data)

badfib_LRG = LRGsummary['FIBER'][np.where(LRGsummary['n_sig_sim'] <= -4)]
badfib_QSO = QSOsummary['FIBER'][np.where(QSOsummary['n_sig_sim'] <= -4)]
badfib_ELG = ELG_LOPnotqsosummary['FIBER'][np.where(ELG_LOPnotqsosummary['n_sig_sim'] <= -4)]
badfib_BGS = BGS_BRIGHTsummary['FIBER'][np.where(BGS_BRIGHTsummary['n_sig_sim'] <= -4)]


np.savetxt('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/kibo-v1/LRG_badfibers.txt',badfib_LRG[np.argsort(badfib_LRG)].astype('int'),fmt='%i')
np.savetxt('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/kibo-v1/QSO_badfibers.txt',badfib_QSO[np.argsort(badfib_QSO)].astype('int'),fmt='%i')
np.savetxt('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/kibo-v1/BGS_BRIGHT_badfibers.txt',badfib_BGS[np.argsort(badfib_BGS)].astype('int'),fmt='%i')
np.savetxt('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/kibo-v1/ELG_LOPnotqso_badfibers.txt',badfib_ELG[np.argsort(badfib_ELG)].astype('int'),fmt='%i')

badfib_all = np.array(list(set(np.concatenate((badfib_LRG,badfib_QSO,badfib_ELG,badfib_BGS)))))

# Add 7 fibers failing the n(z) KS test
badfib_nz = np.array([916, 917, 1261, 1269, 1295, 2246, 3359])
badfib_all = np.concatenate((badfib_all, badfib_nz))

np.savetxt('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/kibo-v1/unique_badfibers.txt',badfib_all[np.argsort(badfib_all)].astype('int'),fmt='%i')
