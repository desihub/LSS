"""
LSS.masks
=====================

Use desiutil utilities to parse LSS masks in the data directory
Ported from desitarget.targetmask
"""
import os.path
from desiutil.bitmask import BitMask
import yaml
from pkg_resources import resource_filename


def load_mask_bits():
    """Load bit definitions from yaml file.
    """
    # ADM Assumes all masks are stored in one standard location.
    ender = os.path.join('data', 'masks.yaml')
    fn = resource_filename('LSS', ender)
    with open(fn) as fx:
        bitdefs = yaml.safe_load(fx)

    return bitdefs

# -convert to BitMask objects
# if bitdefs is None:
#     load_bits()
bitdefs = load_mask_bits()
try:
    elg_mask = BitMask('elg_mask', bitdefs)
    lrg_mask = BitMask('lrg_mask', bitdefs)
    skymap_mask = BitMask('skymap_mask', bitdefs)
except TypeError:
    elg_mask = object()
    lrg_mask = object()
    skymap_mask = object()
