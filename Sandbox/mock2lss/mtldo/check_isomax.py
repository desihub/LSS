import glob
import astropy.table as tb
import numpy as np
import os

path = '/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_256dirs_rea024/Univ025/sv3/dark/orig' 
times = []
for f in glob.glob(os.path.join(path,'*.ecsv')):
        tab = tb.Table.read(f)
        #times.append(np.unique(tab['TIMESTAMP']).data)#[0])
        times.extend(np.unique(tab['TIMESTAMP']).data)#[0])
#        print(np.unique(tab['TIMESTAMP']).data)#[0])
        times = list(np.unique(times))
#        print(times)

print(times)
from astropy.time import Time
Times = Time(times, format='isot')
tt = Times.to_value('mjd')
print(tt)
print(np.max(tt))

a = Time(np.max(tt), format='mjd')

print (a.isot)

