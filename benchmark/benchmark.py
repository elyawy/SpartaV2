from pathlib import Path
import random
import tempfile
import re

from msasim import sailfish as sf

from spartaabc.utility import get_tree_path
from spartaabc.prior_sampler import PriorSampler
from spartaabc.aligner_interface import Aligner
from spartaabc.abc_inference import IndelParams
from spartaabc.spabc import parallelized_inference
from spartaabc.raxml_parser import parse_raxml_bestModel

def scale_tree(tree_path: str, scale_factor: float, overwrite: bool=False):
    # Read the tree
    tree_path = Path(tree_path)
    print(tree_path)

    with open(tree_path, 'r') as f:
        newick_str = f.read().strip()
    
    # Function to scale branch lengths in Newick string
    def scale_length(match):
        length = float(match.group(1))
        scaled = length / scale_factor
        return f":{scaled:.10f}"
    
    # Scale all branch lengths using regex
    scaled_newick = re.sub(r':([0-9.]+)', scale_length, newick_str)

    scaled_tree_path = (tree_path.parent / "scaled.tree") if not overwrite else tree_path
    # Write the scaled tree directly
    with open(scaled_tree_path, 'w') as f:
        f.write(scaled_newick)
    print(f"Scaled tree written to: {scaled_tree_path}")

def create_fake_data_path(data_path: Path) -> Path:
    """
    Prepares folder for ABC pipeline:
    1. Simulates MSA based on tree
    2. Realigns with MAFFT
    
    Args:
        data_path: Path to tree folder
    
    """
    tree_path = get_tree_path(data_path)
    substitution_model = parse_raxml_bestModel(data_path)
    print(substitution_model)
    # scale_tree(tree_path, scale_factor=10.0, overwrite=True)

    seed = random.randint(1,1e6)
    selected_model = "sim" if random.random() < 0.5 else "rim"

    prior_sampler = PriorSampler(indel_model=selected_model, seed=seed)
    sampled_params = prior_sampler.sample()

    root_length = sampled_params[0][0]
    insertion_rate, deletion_rate = sampled_params[0][1]
    lendist, length_dist_insertion, length_dist_deletion = sampled_params[0][2]

    indel_params = IndelParams(root_length,
                               insertion_rate, deletion_rate, 
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
    

    
    sim.set_replacement_model(model=substitution_model["submodel"],
                              gamma_parameters_categories=substitution_model["gamma_cats"],
                              gamma_parameters_alpha=substitution_model["gamma_shape"])
    

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


for idx,dir in enumerate(list(Path("benchmark/data").iterdir())):
    print(idx)
    print(dir.stem)
    # create_fake_data_path(dir) # generate simulated data if missing

# test_path = Path("benchmark/data/BBS20002_20")

# parallelized_inference(test_path, "AA", 100_000, 500, ["sim", "rim"], "mafft", True)
