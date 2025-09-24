import numpy as np

from   cosmo import cosmo
from   scipy.spatial.transform import Rotation as R


def rotate(ras, decs, pos):
    phi        = np.radians(ras)
    theta      = np.pi/2. - np.radians(decs)

    mean_phi   = np.median(phi)
    mean_theta = np.median(theta)

    rot        = R.from_rotvec(-mean_phi * np.array([0, 0, 1]))

    res        = rot.apply(pos)

    rot        = R.from_rotvec((np.pi/2. - mean_theta) * np.array([0, 1, 0]))

    resres     = rot.apply(res)

    return  resres 

def cartesian(ras, decs, zs, rotate=False):
    phi        = np.radians(ras)
    theta      = np.pi/2. - np.radians(decs)

    mean_phi   = np.median(phi)
    mean_theta = np.median(theta)
    
    chis       = cosmo.comoving_distance(zs).value # [Mpc/h].

    zs         = chis * np.cos(theta)
    ys         = chis * np.sin(theta) * np.sin(phi)
    xs         = chis * np.sin(theta) * np.cos(phi)

    pos        = np.c_[xs, ys, zs]

    if rotate:
        pos    = rotate(ras, decs, pos)

    return  pos


if __name__ == '__main__':
    ras  = np.random.uniform(10., 20., 100)
    decs = np.random.uniform(2., 3., 100)
    zs   = np.random.uniform(2., 2.4, 100)

    pos  = cartesian(ras, decs, zs)

    print(pos)
