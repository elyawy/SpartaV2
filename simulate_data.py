import argparse
from pathlib import Path
import numpy as np

import msastats
from msasim import sailfish

from prior_sampler import protocol_updater
from utility import *


def parse_args(arg_list: list[str] | None):
    _parser = argparse.ArgumentParser(allow_abbrev=False)
    _parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
    _parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)
    _parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
    _parser.add_argument('-s','--seed', action='store',metavar="Simulation config" , type=int, required=False)
    _parser.add_argument('-l','--lengthdist', action='store',metavar="Simulation config" , type=str, required=True)
    _parser.add_argument('-m','--model', action='store',metavar="Simulation config" , type=str, required=True)
    _parser.add_argument('-v','--verbose', action='store_true')


    args = _parser.parse_args()
    return args

# TODO: simulate only indels lime in the old Sparta -> major speedups
def main(arg_list: list[str] | None = None):
    args = parse_args(arg_list)
    print(args)

    MAIN_PATH = Path(args.input).resolve()
    SEED = args.seed
    SEQUENCE_TYPE = args.type
    NUM_SIMS = args.numsim
    LENGTH_DISTRIBUTION = args.lengthdist
    INDEL_MODEL = args.model
    KEEP_STATS = args.keep_stats
    VERBOSE = args.verbose

    TREE_PATH = get_tree_path(MAIN_PATH)
    MSA_PATH = get_msa_path(MAIN_PATH)




    prior_sampler = prepare_prior_sampler(MSA_PATH, LENGTH_DISTRIBUTION, INDEL_MODEL)


    sim_params = sim_config.get_random_sim(NUM_SIMULATIONS)

sim_protocol = sailfish.SimProtocol(INPUT_TREE)
sim_protocol.set_seed(SEED)


sim = sailfish.Simulator(sim_protocol, simulation_type=sailfish.SIMULATION_TYPE.PROTEIN)
sim.set_replacement_model(model=sailfish.MODEL_CODES.WAG)

msa_stats: list[list[float]] = []
# protocol_updater(sim_protocol, sim_params[0])

# msa = sim()

for idx, param_set in enumerate(sim_params):
    protocol_updater(sim_protocol, param_set)
    msa = sim()
    sequeces = msa.get_msa().splitlines()[1::2]
    sim_stats = msastats.calculate_msa_stats(sequeces)
    msa_stats.append(sim_stats)
