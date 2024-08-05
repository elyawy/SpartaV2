import numpy as np
from scipy.optimize import least_squares
from scipy.stats import geom



def get_moment_geoemtric(p, mom, truncation):
    """
    Calculation of geoemtric moment
    inputs:
    p - p parameter of the geoemtric distribution
    mom - (optional, by default 1) the moment you would like to calculate (for example mom=1 -> average)
    n_trunc - (optional, by default 150) maximal value of the geoemtric distribution
    cutoff - (optional, by default None) if not none it will give the sum for the distribution until cutoff
    output:
    moment of the distribution (if cutoff!=None it will return the partial sum
    """
    CDF = lambda x: geom.cdf(x,p)

    norm_factor = CDF(truncation)



    moment = 0
    for i in range(1, truncation+1):
        prob_i = geom.pmf(i,p)
        # print(prob_i)

        moment += (prob_i/norm_factor)*(i**mom)

    return moment


def geo_moment_to_p(mom, truncation):
    solution = least_squares(lambda x: (mom - get_moment_geoemtric(x, 1, truncation)),
                        x0=0.5, bounds=(0.0001, 1)).x[0]

    return solution

# m  = get_moment_geoemtric(0.001, 1, 150)
# print(m)