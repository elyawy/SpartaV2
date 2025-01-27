from spartaabc.getting_priors import geometric as geo
from spartaabc.getting_priors import poisson as pois
from spartaabc.getting_priors import zipf


length_distributions = ["zipf", "geometric", "poisson"]


params_search_map = {
    "geometric" :  geo.geo_moment_to_p,
    "zipf":  zipf.zip_mom_to_a,
    "poisson": pois.pois_moment_to_mu
}

mean_check_map = {
    "geometric" :  geo.get_moment_geoemtric,
    "zipf":  zipf.calc_zip_mom,
    "poisson": pois.get_moment_poisson
}

TRUNCATION = 150

lower_mean = 1.5
upper_mean = 25.0

lower_prior_lims = {}
upper_prior_lims = {}


epsilon = 0.001

for len_dist in length_distributions:
    # calculate prior limits (params) for given means:
    lower_prior_lims[len_dist] = params_search_map[len_dist](lower_mean, TRUNCATION)    
    upper_prior_lims[len_dist] = params_search_map[len_dist](upper_mean, TRUNCATION)    

    # get the resulting means from calculated params:
    lower_mean_test = mean_check_map[len_dist](lower_prior_lims[len_dist], 1, TRUNCATION)
    upper_mean_test = mean_check_map[len_dist](upper_prior_lims[len_dist], 1, TRUNCATION)

    # test if calculated param is correct within epsilon threshold:
    if (lower_mean_test < lower_mean - epsilon) or (lower_mean > lower_mean + epsilon):
        print(f"Could not find lower mean param for {len_dist}")

    if (upper_mean_test < upper_mean - epsilon) or (upper_mean_test > upper_mean + epsilon):
        print(f"Could not find upper mean param for {len_dist}")

# organize the results for printing:
priors = tuple(zip(lower_prior_lims.values(), upper_prior_lims.values()))
final_priors = dict(zip(length_distributions, priors))

# print(final_priors)