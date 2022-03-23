#Edmond CHAUSSIDON (CEA Saclay)

import os
import logging
import numpy as np
import healpy as hp

class RFWeight(object):

    logger = logging.getLogger('RFWeight')

    def __init__(self, tracer="LRG", survey='SV3', Nside=None, use_stream=False, dir_weight="/global/cfs/cdirs/desi/users/edmondc/Imaging_weight/"):
        """
        Survey is eihter SV3 or MAIN or DA02.
        The list of tracer available can be found here : /global/cfs/cdirs/desi/users/edmondc/Imaging_weight/

        Remark :
            * SV3 / MAIN are for target selection
            * DA02 is for clustering object (see README.md in the sub-directory for more info)
        """
        if Nside is not None:
            self.nside = Nside
        elif tracer in ["QSO", "ELG_VLO"]:
            self.nside=256
        elif survey == 'DA02':
            self.nside=128
        else:
            self.nside=512

        suffixe=''
        if use_stream:
            suffixe='_with_stream'

        weight_file = os.path.join(dir_weight, f"{survey}/{tracer}_imaging_weight_{self.nside}{suffixe}.npy")
        logging.info(f"Read imaging weight: {weight_file}")
        self.map = np.load(weight_file)

    def __call__(self, ra, dec):
        pix = hp.ang2pix(self.nside, ra, dec, nest=True, lonlat=True)
        return self.map[pix]

if __name__ == '__main__':
    import fitsio

    logger = logging.getLogger()
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    weight = RFWeight(tracer="LRG", survey="DA02", Nside=128)

    catalog_file = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/everest/LSScats/test/LRGzdone_clustering.dat.fits'
    logger.info(f"Read catalog: {catalog_file}")
    catalog = fitsio.read(catalog_file)

    w = weight(catalog['RA'],catalog['DEC'])

    print(f"Info on w: mean = {w.mean()}, std = {w.std()}")
    print(f"Info: #of objects in clustering catalog: {w.size}")
