import tempfile
import random
from pathlib import Path
from spartaabc.aligner_interface import Aligner
import ngesh
from ete3 import Tree
from msasim import sailfish as sf

from spartaabc import correction

def create_fake_data_path(tree: Tree) -> Path:
    """
    Runs a complete phylogenetic analysis pipeline:
    1. Writes tree to temp directory
    2. Simulates MSA based on tree
    3. Realigns with MAFFT
    
    Args:
        tree: Phylogenetic tree object
        simulator: MSA simulation object
        inference_pipeline: Function to run inference
    
    Returns:
        pd.DataFrame: Results from the analysis
    """
    # Create temp directory
    data_path = Path("tests/correction_test/fake_data") / str(len(tree))
    data_path.mkdir(exist_ok=True, parents=True)
    
    # Write tree file
    tree_path = data_path / "tree.newick"
    tree.write(format=1, outfile=str(tree_path))

    seed = random.randint(1,1e6)
    
    sim_protocol = sf.SimProtocol(tree=str(tree_path),
                                  root_seq_size=100, deletion_rate=0.03, insertion_rate=0.03,
                                  deletion_dist=sf.ZipfDistribution(2.2, 150),
                                  insertion_dist=sf.ZipfDistribution(2.2, 150),
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

for num_leaves in range(5,16):

    tree = ngesh.gen_tree(2.0, 2.0, num_leaves=num_leaves, labels="human")
    print(sum([node.dist for node in tree.iter_descendants()]))
    data_path = create_fake_data_path(tree)

    correction.compute_alignment_correction(
        MAIN_PATH=data_path, SEED=5, SEQUENCE_TYPE="AA",
        NUM_SIMS=500, ALIGNER="MAFFT", LENGTH_DISTRIBUTION="zipf",
        INDEL_MODEL="SIM", KEEP_STATS=False
    )



    print(tree)