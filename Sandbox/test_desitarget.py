from desitarget.io import read_targets_in_tiles
import time

import desimodel.io as dmio
tiles = dmio.load_tiles()

t0 = time.time()
tr = read_targets_in_tiles(
    '/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/2.4.0/randoms/resolve/randoms-1-0', tiles[:10])
print(str(time.time()-t0))

t0 = time.time()
tr = read_targets_in_tiles(
    '/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/2.4.0/randoms/resolve/randoms-1-0', tiles[:10], use_concatenate=True)
print(str(time.time()-t0))
