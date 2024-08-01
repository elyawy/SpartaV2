import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import msastats
from msasim import sailfish as sf

from prior_sampler import protocol_updater
from utility import *


def parse_args(arg_list: list[str] | None):
    _parser = argparse.ArgumentParser(allow_abbrev=False)
    _parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
    _parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
    _parser.add_argument('-s','--seed', action='store',metavar="Simulation config" , type=int, required=True)
    _parser.add_argument('-l','--lengthdist', action='store',metavar="Simulation config" , type=str, required=True)
    _parser.add_argument('-m','--model', action='store',metavar="Simulation config" , type=str, required=True)
    _parser.add_argument('-v','--verbose', action='store_true')


    args = _parser.parse_args()
    return args


def simulate_data(prior_sampler: PriorSampler, num_sims: int, tree_path: str, seed: int):
    sim_protocol = sf.SimProtocol(tree=tree_path)
    print(sim_protocol.get_tree().get_num_nodes())
    sim_protocol.set_seed(seed)
    simulator = sf.Simulator(sim_protocol)

    # simulated_msas = []
    sum_stats = []
    root_sampler = prior_sampler.sample_root_length()
    indel_rate_sampler = prior_sampler.sample_rates()
    length_distribution_sampler = prior_sampler.sample_length_distributions()

    for idx in range(num_sims):
        root_length = next(root_sampler)
        insertion_rate, deletion_rate = next(indel_rate_sampler)
        lendist, insertion_length_dist, deletion_length_dist = next(length_distribution_sampler)

        # print(root_length, insertion_rate, deletion_rate, insertion_length_dist, deletion_length_dist)
        print(idx)
        numeric_params = [root_length ,insertion_rate, deletion_rate, insertion_length_dist.p, deletion_length_dist.p]
        protocol_updater(sim_protocol, [root_length, insertion_rate, deletion_rate,
                         insertion_length_dist, deletion_length_dist])

        sim_msa = simulator()
        sim_stats = msastats.calculate_msa_stats(sim_msa.get_msa().splitlines())

        sum_stats.append(numeric_params + sim_stats + [lendist])

    return np.array(sum_stats)

# TODO: simulate only indels lime in the old Sparta -> major speedups
def main(arg_list: list[str] | None = None):
    args = parse_args(arg_list)
    print(args)

    MAIN_PATH = Path(args.input).resolve()
    SEED = args.seed
    NUM_SIMS = args.numsim
    LENGTH_DISTRIBUTION = args.lengthdist
    INDEL_MODEL = args.model
    VERBOSE = args.verbose

    TREE_PATH = get_tree_path(MAIN_PATH)
    MSA_PATH = get_msa_path(MAIN_PATH)


    prior_sampler = prepare_prior_sampler(MSA_PATH, LENGTH_DISTRIBUTION, INDEL_MODEL, SEED)
    msa_stats = simulate_data(prior_sampler, num_sims=NUM_SIMS, tree_path=TREE_PATH, seed=SEED)

    data_full = msa_stats
    data_full = pd.DataFrame(data_full, columns=PARAMS_LIST + SUMSTATS_LIST + ["length_distribution"])
    data_full.to_parquet(MAIN_PATH / f"full_data_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}.parquet.gzip", compression="gzip")


    
if __name__ == "__main__":
    main()