import numpy as np
from scipy.optimize import least_squares
from scipy.stats import poisson


def get_moment_poisson(mu, mom, truncation):
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
    CDF = lambda x: poisson.cdf(x, mu)
    
    norm_factor = CDF(truncation)



    moment = 0
    for i in range(1, truncation+1):
        prob_i = poisson.pmf(i, mu)
        moment += (prob_i/norm_factor)*(i**mom)

    return moment


def pois_moment_to_mu(mom, truncation):
    solution = least_squares(lambda x: (mom - get_moment_poisson(x, 1, truncation)),
                        x0=mom, bounds=(1, 100)).x[0]

    return solution


# m = get_moment_poisson(1.5, 1, 150)
# print(m)