import numpy as np
import os

tracers = ['QSO','LRG','BGS_BRIGHT','ELG_LOPnotqso']
for tp in tracers:
    path = "/global/cfs/cdirs/desi/users/akrolew/summary/"+tp + "/" #path of text files
    concatenated_text = ''
    for file_name in os.listdir(path):
        if file_name.endswith('.txt'):
            if os.path.getmtime(os.path.join(path, file_name)) >= 1723504931.3894317:
                with open(os.path.join(path, file_name), 'r') as f:
                    file_contents = f.read()
                    # Add the file contents to the concatenated text
                    concatenated_text += file_contents
    with open('/global/cfs/cdirs/desi/users/akrolew/summaryscale/'+tp+'fibersim.txt', 'w') as f2:
        f2.write(concatenated_text)
        
from astropy.table import Table
from scipy import special

#tracers = ['QSO']
for tp in tracers:
    data = np.loadtxt('/global/cfs/cdirs/desi/users/akrolew/summaryscale/'+tp+'fibersim.txt',usecols=range(10))
    col_names = ['FIBER', 'p_value', 'var', 'nsig', 'n', 's', 'mean_mod_suc', 'frac_suc', 'mean_x','mean_y']
    exec(tp+'summary = Table(data, names=col_names)')
    exec(tp+'summary ["n_sig_sim"] = np.around((-np.sqrt(2) * special.erfcinv(2 * '+tp+'summary["p_value"])),5)')