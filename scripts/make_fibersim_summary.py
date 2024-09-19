import numpy as np
import os

tracers = ['QSO']
for tp in tracers:
    path = "/global/cfs/cdirs/desi/users/akrolew/summary/"+tp + "/" #path of text files
    concatenated_text = ''
    for file_name in os.listdir(path):
        if file_name.endswith('.txt'):
            if os.path.getmtime(os.path.join(path, file_name)) >= 1726728251.3433545:
                with open(os.path.join(path, file_name), 'r') as f:
                    file_contents = f.read()
                    # Add the file contents to the concatenated text
                    concatenated_text += file_contents
    with open('/global/cfs/cdirs/desi/users/akrolew/summaryscale/'+tp+'fibersim_kibo-v1.txt', 'w') as f2:
        f2.write(concatenated_text)
        
from astropy.table import Table
from scipy import special

for tp in tracers:
    #data1 = np.loadtxt('/global/cfs/cdirs/desi/users/akrolew/summaryscale/'+tp+'fibersim.txt.3',usecols=range(10))
    data2 = np.loadtxt('/global/cfs/cdirs/desi/users/akrolew/summaryscale/'+tp+'fibersim_kibo-v1.txt',usecols=range(10))
    #data3 = np.loadtxt('/global/cfs/cdirs/desicollab/users/akrolew/summary/ELG_LOPnotqso/0.txt',usecols=range(10)).reshape((1,10))
    data = data2
    
    col_names = ['FIBER', 'p_value', 'var', 'nsig', 'n', 's', 'mean_mod_suc', 'frac_suc', 'mean_x','mean_y']
    exec(tp+'summary = Table(data, names=col_names)')
    exec(tp+'summary ["n_sig_sim"] = np.around((-np.sqrt(2) * special.erfcinv(2 * '+tp+'summary["p_value"])),5)')
    
#np.savetxt('/global/cfs/cdirs/desi/users/akrolew/summaryscale/'+tp+'fibersim.txt',data)