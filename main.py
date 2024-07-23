import numpy as np

import msastats
from msasim import sailfish

from prior_sampler import SimConfig

INPUT_MSA = "ENOG5039UWA.fasta"
INPUT_TREE = "ENOG5039UWA.newick"
SEED = 42
NUM_SIMULATIONS = int(1e2)

min_length_index = msastats.stats_names().index("MSA_MIN_LEN")
max_length_index = msastats.stats_names().index("MSA_MAX_LEN")


empirical_stats = msastats.calculate_fasta_stats(INPUT_MSA)
smallest_sequence_size = empirical_stats[min_length_index]
largest_sequence_size = empirical_stats[max_length_index]

sim_config = SimConfig(seq_lengths=[smallest_sequence_size, largest_sequence_size], seed=SEED)

sim_params = sim_config.get_random_sim(NUM_SIMULATIONS)

sim_protocol = sailfish.SimProtocol(INPUT_TREE)
sim_protocol.set_seed(SEED)


def protocol_updater(protocol: sailfish.SimProtocol, params: list) -> None:
    protocol.set_sequence_size(params[0])
    protocol.set_insertion_rates(params[1])
    protocol.set_deletion_rates(params[2])
    protocol.set_insertion_length_distributions(params[4])
    protocol.set_deletion_length_distributions(params[5])


sim = sailfish.Simulator(sim_protocol, simulation_type=sailfish.SIMULATION_TYPE.PROTEIN)
sim.set_replacement_model(model=sailfish.MODEL_CODES.WAG)

msa_stats: list[float] = []
# protocol_updater(sim_protocol, sim_params[0])

# msa = sim()

for idx, param_set in enumerate(sim_params):
    protocol_updater(sim_protocol, param_set)
    print(idx)
    # print(idx,"\nparams=", param_set)
    msa = sim()
    # print("num_sequences=", msa.get_num_sequences())
    # print("length=", msa.get_length())

    sequeces = msa.get_msa().splitlines()[1::2]
    sim_stats = msastats.calculate_msa_stats(sequeces)

    msa_stats.append(sim_stats)
    # break



# print(msa_stats[20])
# print(sim_params[20])