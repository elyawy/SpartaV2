import numpy as np
import pandas as pd


from scipy.optimize import fsolve


def calc_zip_mom(a, mom=1, truncation=50, cutoff=None):
    """
    Calculation of Zipfian moment
    inputs:
    a - a parameter of the Zipfian distribution
    mom - (optional, by default 1) the moment you would like to calculate (for example mom=1 -> average)
    truncation - (optional, by default 150) maximal value of the Zipfian distribution
    cutoff - (optional, by default None) if not none it will give the sum for the distribution until cutoff
    output:
    moment of the distribution (if cutoff!=None it will return the partial sum
    """
    if a <= 1:
        return 10000000
    a = float(a)
    mom = int(mom)
    truncation = int(truncation)
    z = 0
    out = 0
    for i in range(1, truncation + 1, 1):
        z += i ** -a
        if (cutoff is None) or (cutoff is not None and i <= cutoff):
            out += (i ** -a) * (i ** mom)
    return out / z


def zip_mom_to_a(mom, truncation=50, init_guess=1.00000001, epsfcn=0.0000001, xtol=1.49012e-10):
    """
    Returns the "a" parameter of the truncated zipfian function with the input moment
    mom - value of the moment
    mom_order - order of the moment (1 is the average)
    truncation -  where the distribution is truncated
    """
    return fsolve(lambda x: (mom - calc_zip_mom(x, mom=1, truncation=truncation)) \
                  , init_guess, epsfcn=epsfcn, xtol=xtol)[0]