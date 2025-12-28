import logging
import shutil
import sys
from pathlib import Path

# For Python 3.9+
try:
    from importlib.resources import files
except ImportError:
    # For Python < 3.9
    from importlib_resources import files

default_prior_config_path = files("spartaabc").joinpath("default_prior.json")


def validate_input_directory(main_path: Path) -> dict:
    """Validate input directory and provide helpful error messages"""
    issues = []
    
    # Check for required files
    fasta_files = list(main_path.glob("*.fasta"))
    tree_files = list(main_path.glob("*.tree")) + list(main_path.glob("*.newick"))
    
    if not fasta_files:
        issues.append("❌ No FASTA file found. Please provide a .fasta file with your MSA.")
    elif len(fasta_files) > 1:
        issues.append(f"⚠️  Multiple FASTA files found: {[f.name for f in fasta_files]}. Using {fasta_files[0].name}")
    
    if not tree_files:
        issues.append("❌ No tree file found. Please provide a .tree or .newick file.")
    elif len(tree_files) > 1:
        issues.append(f"⚠️  Multiple tree files found: {[f.name for f in tree_files]}. Using {tree_files[0].name}")
    
    # Check for optional .bestModel file
    model_files = list(main_path.glob("*.bestModel"))
    correction_recommended = len(model_files) > 0
    
    if issues:
        print("\n".join(issues))
        if any("❌" in issue for issue in issues):
            sys.exit(1)
    
    return {
        "fasta_file": fasta_files[0] if fasta_files else None,
        "tree_file": tree_files[0] if tree_files else None,
        "model_file": model_files[0] if model_files else None,
        "correction_recommended": correction_recommended
    }


def check_dependencies():
    """Check if all required external tools are available"""
    required_tools = ["mafftpy"]
    
    print("🔍 Checking dependencies...")
    
    for tool in required_tools:
        if shutil.which(tool):
            print(f"✅ {tool}: Found")
        else:
            print(f"❌ {tool}: Not found")
            print(f"   Please make sure that {tool} is installed on your system.")
            sys.exit(1)

def get_tree_path(main_path: Path) -> str:
    tree_path = None
    if len( n := list(main_path.glob("*.tree")) + list(main_path.glob("*.newick"))) == 1:
        tree_path = str(n[0])

    if tree_path is None:
        print("no tree file")
        exit()

    return tree_path

def get_msa_path(main_path: Path) -> str:
    msa_path = None
    if len( n := list(main_path.glob("*.fasta"))) == 1:
        msa_path = str(n[0])

    if msa_path is None:
        print("no fasta file")
        exit()

    return msa_path


def prepare_prior_sampler(empirical_msa_path: str, indel_model:str,
                          seed: int, prior_conf_path: Path):
    import msastats
    from spartaabc.prior_sampler import PriorSampler

    MIN_LENGTH_STAT_INDEX = msastats.stats_names().index("MSA_MIN_LEN")
    MAX_LENGTH_STAT_INDEX = msastats.stats_names().index("MSA_MAX_LEN")

    empirical_stats = msastats.calculate_fasta_stats(empirical_msa_path)
    smallest_sequence_size = empirical_stats[MIN_LENGTH_STAT_INDEX]
    largest_sequence_size = empirical_stats[MAX_LENGTH_STAT_INDEX]

    seq_lengths_in_msa = [smallest_sequence_size, largest_sequence_size]

    prior_sampler = PriorSampler(conf_file=prior_conf_path,
                        seq_lengths=seq_lengths_in_msa,
                        indel_model=indel_model,
                        seed=seed)
    return prior_sampler

def parse_model_params(file_path: Path):
    if not file_path.exists():
        print("No inferred indel model file found!\n Try rerunning spartaabc.")
        sys.exit(1)
    model = {}

    model_ = file_path.read_text().splitlines()
    model["indel_model"] = (model_[0].split(": ")[1])

    model_ = [line.split(": ") for line in model_[1:]]
    model_ = {key: float(val) for (key,val) in model_}

    model["root_length"] = int(model_.get("Root_length"))
    model["insertion_rate"] = model_.get("R_I") or model_.get("R_ID")
    model["deletion_rate"] = model_.get("R_D") or model_.get("R_ID")
    model["length_param_insertion"] = model_.get("A_I") or model_.get("A_ID")
    model["length_param_deletion"] = model_.get("A_D") or model_.get("A_ID")

    return model

PARAMS_LIST = [
    "root_length",
    "insertion_rate",
    "deletion_rate",
    "length_param_insertion",
    "length_param_deletion"
]

SUMSTATS_LIST = [f'SS_{i}' for i in range(0,27)]
SUMSTATS_DEFINITION = {
    'SS_0': "AVG_GAP_SIZE",
    'SS_1': "MSA_LEN",
    'SS_2': "MSA_MAX_LEN",
    'SS_3': "MSA_MIN_LEN",
    'SS_4': "TOT_NUM_GAPS",
    'SS_5': "NUM_GAPS_LEN_ONE",
    'SS_6': "NUM_GAPS_LEN_TWO",
    'SS_7': "NUM_GAPS_LEN_THREE",
    'SS_8': "NUM_GAPS_LEN_AT_LEAST_FOUR",
    'SS_9': "AVG_UNIQUE_GAP_SIZE",
    'SS_10': "TOT_NUM_UNIQUE_GAPS",
    'SS_11': "NUM_GAPS_LEN_ONE\nPOS_1_GAPS",
    'SS_12': "NUM_GAPS_LEN_ONE\nPOS_2_GAPS",
    'SS_13': "NUM_GAPS_LEN_ONE\nPOS_N_MINUS_1_GAPS",
    'SS_14': "NUM_GAPS_LEN_TWO\nPOS_1_GAPS",
    'SS_15': "NUM_GAPS_LEN_TWO\nPOS_2_GAPS",
    'SS_16': "NUM_GAPS_LEN_TWO\nPOS_N_MINUS_1_GAPS",
    'SS_17': "NUM_GAPS_LEN_THREE\nPOS_1_GAPS",
    'SS_18': "NUM_GAPS_LEN_THREE\nPOS_2_GAPS",
    'SS_19': "NUM_GAPS_LEN_THREE\nPOS_N_MINUS_1_GAPS",
    'SS_20': "NUM_GAPS_LEN_AT_LEAST_FOUR\nPOS_1_GAPS",
    'SS_21': "NUM_GAPS_LEN_AT_LEAST_FOUR\nPOS_2_GAPS",
    'SS_22': "NUM_GAPS_LEN_AT_LEAST_FOUR\nPOS_N_MINUS_1_GAPS",
    'SS_23': "MSA_POSITION_WITH_0_GAPS",
    'SS_24': "MSA_POSITION_WITH_1_GAPS",
    'SS_25': "MSA_POSITION_WITH_2_GAPS",
    'SS_26': "MSA_POSITION_WITH_N_MINUS_1_GAPS"
}




logger = logging.getLogger('main')
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s[%(levelname)s][%(filename)s][%(funcName)s]: %(message)s')



def setLogHandler(path: Path, mode: str="a"):
    handler = logging.FileHandler(path / 'info.log', mode=mode)  # Adjust the path
    handler.setFormatter(formatter)
    handler.setLevel(logging.INFO)
    logger.addHandler(handler)

    # handler = logging.FileHandler(path / 'error.log')  # Adjust the path
    # handler.setFormatter(formatter)
    # handler.setLevel(logging.ERROR)
    # logger.addHandler(handler)


