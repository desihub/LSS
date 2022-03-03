import glob
import astropy.table as tb
import numpy as np
times = []
for f in glob.glob('*.ecsv'):
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

