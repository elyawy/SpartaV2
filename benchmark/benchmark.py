from pathlib import Path
import random
import tempfile

from msasim import sailfish as sf

from spartaabc.utility import get_tree_path
from spartaabc.prior_sampler import PriorSampler
from spartaabc.aligner_interface import Aligner
from spartaabc.abc_inference import IndelParams
from spartaabc.spabc import parallelized_inference

def create_fake_data_path(data_path: Path) -> Path:
    """
    Prepares folder for ABC pipeline:
    1. Simulates MSA based on tree
    2. Realigns with MAFFT
    
    Args:
        data_path: Path to tree folder
    
    """
    tree_path = get_tree_path(data_path)
    seed = random.randint(1,1e6)
    selected_model = "sim" if random.random() < 0.5 else "rim"

    prior_sampler = PriorSampler(indel_model=selected_model, seed=seed)
    sampled_params = prior_sampler.sample()

    root_length = sampled_params[0][0]
    insertion_rate, deletion_rate = sampled_params[0][1]
    lendist, length_dist_insertion, length_dist_deletion = sampled_params[0][2]

    indel_params = IndelParams(insertion_rate, deletion_rate, 
                               length_dist_insertion.p, length_dist_deletion.p,
                               lendist, selected_model)
    print(indel_params)
    (data_path/ "true_params.txt").write_text(str(indel_params) + f"\nSEED: {seed}" )

    sim_protocol = sf.SimProtocol(tree=str(tree_path),
                                  root_seq_size=root_length,
                                  deletion_rate=deletion_rate,
                                  insertion_rate=insertion_rate,
                                  deletion_dist=length_dist_deletion,
                                  insertion_dist=length_dist_insertion,
                                  seed=seed)

    sim = sf.Simulator(simProtocol=sim_protocol,
                       simulation_type=sf.SIMULATION_TYPE.PROTEIN)
    
    sim.set_replacement_model(sf.MODEL_CODES.WAG)

    # Simulate MSA
    simulated_msa = sim()
    
    # Run MAFFT realignment
    sequence_aligner = Aligner("MAFFT")
    msa_string = simulated_msa.get_msa()

    sim_fasta_unaligned = msa_string.replace("-","").encode()
    with tempfile.NamedTemporaryFile(suffix='.fasta') as tempf:
        tempf.write(sim_fasta_unaligned)
        tempf.seek(0)
        sequence_aligner.set_input_file(tempf.name, tree_file=tree_path)
        realigned_msa = sequence_aligner.get_realigned_msa()

    (data_path / "realigned.fasta").write_text(realigned_msa)


    return data_path        


# for idx,dir in enumerate(list(Path("benchmark/data").iterdir())):
#     print(idx)
#     print(dir.stem)
    # create_fake_data_path(dir) # generate simulated data if missing

test_path = Path("benchmark/data/BBS11035_5")

parallelized_inference(test_path, "AA", 100_000, 500, ["sim", "rim"], "mafft", True)
