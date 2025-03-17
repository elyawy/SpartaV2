import pickle, pathlib, os
import numpy as np
import random

from msasim import sailfish as sf

from spartaabc.getting_priors import get_means

TRUNCATION = 150

length_distribution_priors = {
    "zipf": {
        "insertion": sorted(get_means.final_priors["zipf"]),
        "deletion": sorted(get_means.final_priors["zipf"])
    },
    "geometric": {
        "insertion": sorted(get_means.final_priors["geometric"]),
        "deletion": sorted(get_means.final_priors["geometric"])
    },
    "poisson": {
        "insertion": sorted(get_means.final_priors["poisson"]),
        "deletion": sorted(get_means.final_priors["poisson"])
    }
}


def fast_zipf(a_param, truncation=150):
    harmonic_series = np.arange(1,truncation+1)
    harmonic_series = np.power(harmonic_series, -a_param)
    harmonic_sum = np.sum(harmonic_series)
    return harmonic_series / harmonic_sum

length_dist_mapper = {
    "zipf": sf.CustomDistribution,
    "poisson": sf.PoissonDistribution,
    "geometric": sf.GeometricDistribution
}

def protocol_updater(protocol: sf.SimProtocol, params: list) -> None:
    protocol.set_sequence_size(params[0])
    protocol.set_insertion_rates(insertion_rate=params[1])
    protocol.set_deletion_rates(deletion_rate=params[2])
    protocol.set_insertion_length_distributions(insertion_dist=params[3])
    protocol.set_deletion_length_distributions(deletion_dist=params[4])

class PriorSampler:
    def __init__(self, conf_file=None,
                       len_dist="zipf",
                       rate_priors=[[-4,-1],[-1,1]], # log
                       seq_lengths=[100,500],
                       indel_model="sim",
                       seed = 1):
        self.seed = seed
        random.seed(seed)
        self.indel_model = indel_model

        self.length_distribution = length_dist_mapper[len_dist]
        self.len_dist = len_dist

        self.len_prior_dict = length_distribution_priors[len_dist]

        self.rate_prior_dict = {
            "sum_rates": rate_priors[0],
            "ratio_rates": rate_priors[1]
        }

        self.sequence_length_prior = [int(seq_lengths[0]*0.8), int(seq_lengths[1]*1.1)]

    def sample_root_length(self):
        while True:
            root_length = random.randint(*self.sequence_length_prior)
            yield root_length

    def sample_length_distributions(self):
        while True:
            x = random.uniform(*self.len_prior_dict["insertion"])
            if self.indel_model == "sim":
                indel_length_dist = self.length_distribution(fast_zipf(x, TRUNCATION))
                indel_length_dist.p = x
                yield self.len_dist, indel_length_dist, indel_length_dist
            else:
                y = random.uniform(*self.len_prior_dict["deletion"])
                indel_length_dist_insertion = self.length_distribution(fast_zipf(x, TRUNCATION))
                indel_length_dist_insertion.p = x
                indel_length_dist_deletion = self.length_distribution(fast_zipf(y, TRUNCATION))
                indel_length_dist_deletion.p = y
                yield self.len_dist, indel_length_dist_insertion, indel_length_dist_deletion



    def sample_rates(self):
        while True:
            sum_of_rates = 10**random.uniform(*self.rate_prior_dict["sum_rates"])
            ratio_of_rates = 10**random.uniform(*self.rate_prior_dict["ratio_rates"])
            if self.indel_model == "sim":
                yield (sum_of_rates, sum_of_rates)
            else:
                deletion_rate = sum_of_rates/(ratio_of_rates+1)
                insertion_rate = sum_of_rates - deletion_rate
                yield (insertion_rate, deletion_rate)


    def sample(self, n=1):
        root_length = self.sample_root_length()
        indel_rates = self.sample_rates()
        length_dists = self.sample_length_distributions()
        params_sample = []
        for params in zip(root_length, indel_rates, length_dists):
            if n==0:
                break
            params_sample.append(params)
            n = n - 1
        return params_sample