import argparse
import logging
import time
from pathlib import Path
import pkg_resources

import numpy as np
import pandas as pd

import msastats
from msasim import sailfish as sf

from spartaabc.prior_sampler import protocol_updater
from spartaabc.utility import logger, setLogHandler
from spartaabc.utility import PriorSampler, prepare_prior_sampler
from spartaabc.utility import get_msa_path, get_tree_path
from spartaabc.utility import PARAMS_LIST, SUMSTATS_LIST
from spartaabc.utility import default_prior_config_path


default_config_path = pkg_resources.resource_filename("spartaabc", "default_prior.json")


def parse_args(arg_list: list[str] | None):
    _parser = argparse.ArgumentParser(allow_abbrev=False)
    _parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
    _parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
    _parser.add_argument('-s','--seed', action='store',metavar="Simulation config" , type=int, required=False)
    _parser.add_argument('-m','--model', action='store',metavar="Simulation config" , type=str, required=True)
    _parser.add_argument('-p','--prior', action='store',metavar="Prior config path" , type=str, required=False, default=default_prior_config_path)

    _parser.add_argument('-v','--verbose', action='store_true')


    args = _parser.parse_args(arg_list)
    return args


def simulate_data(prior_sampler: PriorSampler, num_sims: int, tree_path: str, seed: int):
    logger.debug(f"num_sims {num_sims} tree_path {tree_path} seed {seed}")
    sim_protocol = sf.SimProtocol(tree=tree_path)

    sim_protocol.set_seed(seed)
    simulator = sf.Simulator(sim_protocol,simulation_type=sf.SIMULATION_TYPE.NOSUBS)

    # simulated_msas = []
    sum_stats = []
    root_sampler = prior_sampler.sample_root_length()
    indel_rate_sampler = prior_sampler.sample_rates()
    length_distribution_sampler = prior_sampler.sample_length_distributions()

    logger.info("Starting msa simulation")
    for _ in range(num_sims):
        root_length = next(root_sampler)
        insertion_rate, deletion_rate = next(indel_rate_sampler)
        lendist, insertion_length_dist, deletion_length_dist = next(length_distribution_sampler)

        numeric_params = [root_length ,insertion_rate, deletion_rate, insertion_length_dist.p, deletion_length_dist.p]
        protocol_updater(sim_protocol, [root_length, insertion_rate, deletion_rate,
                         insertion_length_dist, deletion_length_dist])

        sim_msa = simulator()
        sim_stats = msastats.calculate_msa_stats(sim_msa.get_msa().splitlines())

        sum_stats.append(numeric_params + sim_stats)
    logger.info(f"Done with {num_sims} msa simulations")

    return np.array(sum_stats)

def main(arg_list: list[str] | None = None):
    logging.basicConfig()
    args = parse_args(arg_list)

    MAIN_PATH = Path(args.input).resolve()
    SEED = args.seed if args.seed else time.time_ns()

    NUM_SIMS = args.numsim
    INDEL_MODEL = args.model
    PRIOR_PATH = Path(args.prior).resolve()

    VERBOSE = args.verbose
    
    setLogHandler(MAIN_PATH)
    logger.info("\n\tMAIN_PATH: {},\n\tSEED: {}, NUM_SIMS: {},\n\tINDEL_MODEL {}".format(
        MAIN_PATH, SEED, NUM_SIMS, INDEL_MODEL
    ))


    TREE_PATH = get_tree_path(MAIN_PATH)
    MSA_PATH = get_msa_path(MAIN_PATH)


    prior_sampler = prepare_prior_sampler(MSA_PATH, INDEL_MODEL, SEED, PRIOR_PATH)
    logger.info("\nLoaded prior configuration from file {} is:\n\t{}".format(PRIOR_PATH, prior_sampler))
    LENGTH_DISTRIBUTION = prior_sampler.len_dist

    msa_stats = simulate_data(prior_sampler, num_sims=NUM_SIMS, tree_path=TREE_PATH, seed=SEED)

    data_full = msa_stats
    data_full = pd.DataFrame(data_full, columns=PARAMS_LIST + SUMSTATS_LIST)
    data_full.to_parquet(MAIN_PATH / f"full_data_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}.parquet.gzip", compression="gzip")


    
if __name__ == "__main__":
    main()