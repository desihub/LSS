import numpy as np

def bn(n):
    # https://en.wikipedia.org/wiki/Sérsic_profile
    return 2. * n - 1./3. + 4./405./n + 46./25515./n/n # + ...

def sersic()
