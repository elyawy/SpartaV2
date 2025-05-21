import logging
from pathlib import Path
from sklearn.base import BaseEstimator, TransformerMixin# define the transformer
import msastats

from spartaabc.prior_sampler import PriorSampler

# For Python 3.9+
try:
    from importlib.resources import files
except ImportError:
    # For Python < 3.9
    from importlib_resources import files

default_prior_config_path = files("spartaabc").joinpath("default_prior.json")




MIN_LENGTH_STAT_INDEX = msastats.stats_names().index("MSA_MIN_LEN")
MAX_LENGTH_STAT_INDEX = msastats.stats_names().index("MSA_MAX_LEN")


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
    
    empirical_stats = msastats.calculate_fasta_stats(empirical_msa_path)
    smallest_sequence_size = empirical_stats[MIN_LENGTH_STAT_INDEX]
    largest_sequence_size = empirical_stats[MAX_LENGTH_STAT_INDEX]

    seq_lengths_in_msa = [smallest_sequence_size, largest_sequence_size]

    prior_sampler = PriorSampler(conf_file=prior_conf_path,
                        seq_lengths=seq_lengths_in_msa,
                        indel_model=indel_model,
                        seed=seed)
    return prior_sampler

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


class StandardMemoryScaler(BaseEstimator, TransformerMixin):

    def __init__(self, epsilon=1e-4):
        self._epsilon = epsilon
        
    def fit(self, X, y = None):
        self._mean = X.mean()
        self._std = X.std()

        return self

    def transform(self, X):
        X = (X-self._mean)/(self._std+self._epsilon)
       
        return X



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


