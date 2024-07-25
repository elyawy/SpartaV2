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


def simulate_data(sim_params, num_sims: int, tree_path: str, seed: int):
    sim_protocol = sf.SimProtocol(tree=tree_path)
    sim_protocol.set_seed(seed)
    simulator = sf.Simulator(sim_protocol,
                             simulation_type=sf.SIMULATION_TYPE.PROTEIN)

    # simulated_msas = []
    sum_stats = []

    simulator.set_replacement_model(model=sf.MODEL_CODES.WAG,
                                    gamma_parameters_alpha=1.0,
                                    gamma_parameters_catergories=1)


    for idx,params in enumerate(sim_params):
        numeric_params = [params[0],params[1], params[2], params[4].p, params[5].p]
        print(idx)
        protocol_updater(sim_protocol, params)

        sim_msa = simulator()
        sim_stats = msastats.calculate_msa_stats(sim_msa.get_msa().splitlines()[1::2])
        # print(sim_stats)
        # simulated_msas.append(sim_msa)
        sum_stats.append(numeric_params + sim_stats)    

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
    sim_params = prior_sampler.sample(NUM_SIMS)

    msa_stats = simulate_data(sim_params, num_sims=NUM_SIMS, tree_path=TREE_PATH, seed=SEED)

    data_full = np.concatenate([msa_stats, sim_params[:,[3]]], axis=1)
    data_full = pd.DataFrame(data_full, columns=PARAMS_LIST + SUMSTATS_LIST + ["length_distribution"])
    data_full.to_parquet(MAIN_PATH / f"full_data_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}.parquet.gzip", compression="gzip")

    
if __name__ == "__main__":
    main()