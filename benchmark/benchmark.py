from pathlib import Path
import random
import tempfile
import re

from msasim import sailfish as sf

from spartaabc.utility import get_tree_path
from spartaabc.prior_sampler import PriorSampler
from spartaabc.aligner_interface import Aligner
from spartaabc.abc_inference import IndelParams
from spartaabc.raxml_parser import get_substitution_model

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

def create_fake_data_path(data_path: Path, indel_model: str) -> Path:
    """
    Prepares folder for ABC pipeline:
    1. Simulates MSA based on tree
    2. Realigns with MAFFT
    
    Args:
        data_path: Path to tree folder
    
    """
    tree_path = get_tree_path(data_path)
    substitution_model = get_substitution_model(data_path)
    print(substitution_model)
    # scale_tree(tree_path, scale_factor=2.0, overwrite=True)

    seed = random.randint(1,1e6)
    selected_model = indel_model

    prior_sampler = PriorSampler(indel_model=selected_model,
                                 rate_priors=[[0.0,0.05],[-1.0,1.0]],
                                 seq_lengths=[100,1000],
                                 seed=seed)
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
    sequence_aligner = Aligner("MAFFTFAST")
    msa_string = simulated_msa.get_msa()
    (data_path / "original_msa.txt").write_text(msa_string)

    sim_fasta_unaligned = msa_string.replace("-","").encode()
    with tempfile.NamedTemporaryFile(suffix='.fasta') as tempf:
        tempf.write(sim_fasta_unaligned)
        tempf.seek(0)
        sequence_aligner.set_input_file(tempf.name, tree_file=tree_path)
        realigned_msa = sequence_aligner.get_realigned_msa()

    (data_path / "realigned.fasta").write_text(realigned_msa)


    return data_path        


random.seed(42)
for indel_model_ in ["sim", "rim"]:
    for idx,dir in enumerate(list(Path(f"benchmark/data_fast_{indel_model_}").iterdir())):
        print(idx)
        print(dir.stem)
        create_fake_data_path(dir, indel_model_) # generate simulated data if missing
