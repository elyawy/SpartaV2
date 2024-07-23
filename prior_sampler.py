import pickle, pathlib, os
import numpy as np

from msasim import sailfish as sf

from getting_priors import get_means


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

length_dist_mapper = {
    "zipf": sf.ZipfDistribution,
    "poisson": sf.PoissonDistribution,
    "geometric": sf.GeometricDistribution
}

def protocol_updater(protocol: sf.SimProtocol, params: list) -> None:
    protocol.set_sequence_size(params[0])
    protocol.set_insertion_rates(params[1])
    protocol.set_deletion_rates(params[2])
    protocol.set_insertion_length_distributions(params[4])
    protocol.set_deletion_length_distributions(params[5])

class SimConfig:
    def __init__(self, conf_file=None,
                       len_dist="zipf",
                       rate_priors=[[-4,-1],[-1,1]], # log
                       seq_lengths=[100,500],
                       indel_model="sim",
                       seed = 1):
        self.seed = seed
        np.random.seed(seed)
        self.indel_model = indel_model

        self.length_distribution = length_dist_mapper[len_dist]


        self.len_prior_dict = length_distribution_priors[len_dist]

        self.rate_prior_dict = {
            "sum_rates": rate_priors[0],
            "ratio_rates": rate_priors[1]
        }

        self.sequence_length_prior = [seq_lengths[0]*0.8, seq_lengths[1]*1.1]




    def get_random_sim(self, num_sims):
        
        insertion_lengths = [self.length_distribution(i) for i in np.random.uniform(*self.len_prior_dict["insertion"], num_sims)]
        sum_of_rates = 10**np.random.uniform(*self.rate_prior_dict["sum_rates"], num_sims)
        ratio_of_rates = 10**np.random.uniform(*self.rate_prior_dict["ratio_rates"], num_sims)


        if self.indel_model == "rim":
            deletion_lengths = [self.length_distribution(i) for i in np.random.uniform(*self.len_prior_dict["deletion"], num_sims)]
            deletion_rates = sum_of_rates/(ratio_of_rates+1)
            insertion_rates = sum_of_rates - deletion_rates

        elif self.indel_model == "sim":
            deletion_lengths = insertion_lengths
            deletion_rates = sum_of_rates
            insertion_rates = sum_of_rates


        root_lengths = np.random.randint(*self.sequence_length_prior, num_sims)

        length_distribution = np.repeat(self.length_distribution, num_sims)

        #     root_lengths: int
        #     insertion_rates: float
        #     deletion_rates: float
        #     length_distribution: str
        #     insertion_lengths: float
        #     deletion_lengths: float

        return np.array([
                root_lengths,
                insertion_rates,
                deletion_rates,
                length_distribution,
                insertion_lengths, deletion_lengths
                ], dtype=object).T
