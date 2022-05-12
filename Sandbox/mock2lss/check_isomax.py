import glob
import astropy.table as tb
import numpy as np
import os

from astropy.time import Time

path = '/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_256dirs_rea{REA}/Univ{UNIV}/sv3/dark/orig' 

for i in range(1):
    id_ = "%03d"%int(i)
    un = "%03d"%int(i+1)
    times = []
    pa = path.format(REA=id_, UNIV=un)
    for f in glob.glob(os.path.join(pa,'*.ecsv')):
            tab = tb.Table.read(f)
            times.extend(np.unique(tab['TIMESTAMP']).data)#[0])
            times = list(np.unique(times))

    Times = Time(times, format='isot')
    tt = Times.to_value('mjd')

    a = Time(np.max(tt), format='mjd')

    print (id_, a.isot)

