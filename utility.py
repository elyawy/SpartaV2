from pathlib import Path
import msastats

from prior_sampler import PriorSampler

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


def prepare_prior_sampler(empirical_msa_path: str, length_distribution: str, indel_model:str):
    
    empirical_stats = msastats.calculate_fasta_stats(empirical_msa_path)
    smallest_sequence_size = empirical_stats[MIN_LENGTH_STAT_INDEX]
    largest_sequence_size = empirical_stats[MAX_LENGTH_STAT_INDEX]

    seq_lengths_in_msa = [smallest_sequence_size, largest_sequence_size]

    prior_sampler = PriorSampler(seq_lengths=seq_lengths_in_msa,
                        len_dist=length_distribution,
                        indel_model=indel_model)
    return prior_sampler