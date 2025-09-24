import numpy as np

# Table 3 of https://arxiv.org/pdf/1409.4681.pdf
d8_limits = [[-1.0, -0.75],\
             [-0.75, -0.55],\
             [-0.55, -0.40],\
             [-0.4, 0.0],\
             [0.0, 0.7],\
             [0.7, 1.60],\
             [1.60, 2.90],\
             [2.90, 4.00],\
             [4.0, 1.e4]]

d8_limits = np.array(d8_limits)

def delta8_tier(delta8):
    result = -99 * np.ones(len(delta8), dtype=np.int)

    for i, lims in enumerate(d8_limits):
        result[(delta8 >= lims[0]) & (delta8 < lims[1])] = i

    return  result
