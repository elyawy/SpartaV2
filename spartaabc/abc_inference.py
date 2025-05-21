import argparse
import pickle
import logging
from pathlib import Path
import numpy as np
import pandas as pd
from dataclasses import dataclass

import msastats
from spartaabc.utility import get_msa_path
from spartaabc.utility import PARAMS_LIST, SUMSTATS_LIST
from spartaabc.utility import logger, setLogHandler



@dataclass
class IndelParams:
    root_length: int
    insertion_rate: float
    deletion_rate: float
    insertion_length_parameter: float
    deletion_length_parameter: float
    length_distribution: str
    indel_model: str

    def __repr__(self):
        model_str = f"Model: {self.indel_model}\n"
        model_str += f"Root_length: {self.root_length}\n"
        if self.indel_model in ["SIM", "sim"]:
            model_str += f"R_ID: {self.insertion_rate}\n"
            model_str += f"A_ID: {self.insertion_length_parameter}"
        elif self.indel_model in ["RIM", "rim"]:
            model_str += f"R_I: {self.insertion_rate}\n"
            model_str += f"R_D: {self.deletion_rate}\n"
            model_str += f"A_I: {self.insertion_length_parameter}\n"
            model_str += f"A_D: {self.deletion_length_parameter}"

        return model_str

def parse_args(arg_list: list[str] | None):
    _parser = argparse.ArgumentParser(allow_abbrev=False)
    _parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
    _parser.add_argument('-a','--aligner', action='store',metavar="Aligner", type=str,default="mafft" , required=False)
    _parser.add_argument('-d','--distance', action='store',metavar="Distance metric", type=str, default="mahal", required=False)
    _parser.add_argument('-noc','--no-correction', action='store_false')
    _parser.add_argument('-s', '--exclude_stats', nargs= '+', type=str, default=[], help='List of stats to exclude from original')

    args = _parser.parse_args()
    return args


def load_data(main_path: Path):
    full_data = {}
    for data_path in main_path.glob("*.parquet.gzip"):
        model = data_path.stem.split('.')[0].split("_", maxsplit=2)[2]
        temp_df = pd.read_parquet(data_path)
        full_data[model] = temp_df

    return full_data

def load_correction_regressors(main_path: Path, aligner: str):
    regressors = {}
    for regressor_path in (main_path / f"{aligner}_correction").glob("*.pickle"):
        model = regressor_path.stem.split("_", maxsplit=1)[1]
        with open(regressor_path, 'rb') as f:
            regressors[model] = pickle.load(f)
    return regressors

def load_correction_regressor_scores(main_path: Path, aligner: str):
    scores = pd.DataFrame(len(SUMSTATS_LIST)*[1.0], columns=["pearsonr"])
    for score_path in (main_path / f"{aligner}_correction").glob("*.csv"):
        score_df = pd.read_csv(score_path, index_col=0)[["pearsonr"]]
        scores[scores["pearsonr"] > score_df["pearsonr"]] = score_df

    return scores["pearsonr"].to_list()

def bias_correction(regressors, data: pd.DataFrame, regressor_scores: list[float], kept_statistics: list[int], r_threshold=0.8):
    data = data.to_numpy()

    kept_stats = []
    infered_data = []
    for idx, regressor in enumerate(regressors):
        if idx not in kept_statistics:
            continue
        if regressor_scores[idx] > r_threshold:
            kept_stats.append(idx)
            infered_data.append(regressor.predict(data).T)

    temp_data = np.array(infered_data)
    temp_data = pd.DataFrame(temp_data.T, columns=[SUMSTATS_LIST[i] for i in kept_stats])

    inferred_realigned_stats = temp_data.iloc[:, [kept_stats.index(i) for i in kept_statistics if i in kept_stats]]

    return inferred_realigned_stats, kept_stats

def correct_and_merge_models_data(main_path: Path, aligner: str, kept_statistics, parameters_and_stats) -> tuple[pd.DataFrame, list[int]]:
    regressors = load_correction_regressors(main_path, aligner)
    regressor_scores = load_correction_regressor_scores(main_path, aligner)
    
    stats_data = []
    for model in  parameters_and_stats.keys():
        current_regressors = regressors.get(model, None)

        realigned_statistics, kept_statistics = bias_correction(current_regressors,
                                                                parameters_and_stats[model],
                                                                regressor_scores, kept_statistics)
        stats_data.append(realigned_statistics)

    return pd.concat(stats_data), kept_statistics




def run(main_path: Path, aligner: str, distance_metric: str="mahal", correction=True, top_cutoff: int=100, exclude_stats: list[str] = []) -> IndelParams:

    invalid = [stat for stat in exclude_stats if stat not in SUMSTATS_LIST]
    if invalid:
        raise ValueError(f"Invalid summary stats to exclude: {invalid}")

    MSA_PATH = get_msa_path(main_path)

    kept_statistics_indices = [i for i, stat in enumerate(SUMSTATS_LIST) if stat not in exclude_stats]
    SUMSTATS_LIST_SUBSET = [f"SS_{i}" for i in kept_statistics_indices]

    empirical_stats = msastats.calculate_fasta_stats(MSA_PATH)

    parameters_and_stats = load_data(main_path)
    params_df = pd.concat([parameters_and_stats[model][PARAMS_LIST] for model in  parameters_and_stats.keys()])


    if correction:
        stats_df, kept_statistics_indices = correct_and_merge_models_data(main_path, aligner,
                                                                          kept_statistics_indices,
                                                                          parameters_and_stats)
        stats_df.drop(columns=exclude_stats, inplace=True, errors="ignore")
    else:
        stats_df = pd.concat([parameters_and_stats[model][SUMSTATS_LIST_SUBSET] for model in  parameters_and_stats.keys()])

    empirical_stats = [empirical_stats[i] for i in kept_statistics_indices]

    calculated_distances = None
    if distance_metric == "mahal":
        cov = np.cov(stats_df.T)
        cov = cov + np.eye(len(cov))*1e-4
        inv_covmat = np.linalg.inv(cov)
        u_minus_v = empirical_stats-stats_df
        left = np.dot(u_minus_v, inv_covmat)
        calculated_distances = np.sqrt(np.sum(u_minus_v*left, axis=1))
    if distance_metric == "euclid":
        weights = 1/(stats_df.std(axis=0) + 0.001)
        calculated_distances = np.sum(weights*(stats_df - empirical_stats)**2, axis=1)
    
    stats_df["distances"] = calculated_distances
    stats_df[PARAMS_LIST] = params_df

    top_stats = stats_df.nsmallest(top_cutoff, "distances")

    top_stats[["distances"] + PARAMS_LIST].to_csv(main_path / "top_params.csv", index=False)

    full_sim_data = stats_df[stats_df["insertion_rate"] == stats_df["deletion_rate"]]
    top_sim_data = full_sim_data.nsmallest(top_cutoff, "distances")
    top_sim_data[["distances"] + PARAMS_LIST].to_csv(main_path / "top_params_sim.csv", index=False)

    full_rim_data = stats_df[stats_df["insertion_rate"] != stats_df["deletion_rate"]]
    top_rim_data = full_rim_data.nsmallest(top_cutoff, "distances")
    top_rim_data[["distances"] + PARAMS_LIST].to_csv(main_path / "top_params_rim.csv", index=False)


    abc_indel_params = None
    if len(top_stats[top_stats["insertion_rate"] == top_stats["deletion_rate"]]) > (top_cutoff // 2):
        root_length = int(top_sim_data["root_length"].mean())
        R_ID = float(top_sim_data["insertion_rate"].mean())
        A_ID = float(top_sim_data["length_param_insertion"].mean())
        abc_indel_params = IndelParams(root_length,
                                       R_ID, R_ID,
                                       A_ID, A_ID,
                                       length_distribution="zipf",
                                       indel_model="SIM")
    else:
        root_length = int(top_rim_data["root_length"].mean())
        R_I = float(top_rim_data["insertion_rate"].mean())
        R_D = float(top_rim_data["deletion_rate"].mean())
        A_I = float(top_rim_data["length_param_insertion"].mean())
        A_D = float(top_rim_data["length_param_deletion"].mean())
        abc_indel_params = IndelParams(root_length,
                                       R_I, R_D,
                                       A_I, A_D,
                                       length_distribution="zipf",
                                       indel_model="RIM")
    (main_path / "model_params.txt").write_text(str(abc_indel_params))
    return abc_indel_params

def main(arg_list: list[str] | None = None):
    logging.basicConfig()

    args = parse_args(arg_list)


    MAIN_PATH = Path(args.input).resolve()
    ALIGNER = args.aligner
    DISTANCE_METRIC = args.distance
    CORRECTION = args.no_correction
    EXCLUDE_STATS= args.exclude_stats
    print(MAIN_PATH)
    
    setLogHandler(MAIN_PATH)
    logger.info("\n\tMAIN_PATH: {}".format(
        MAIN_PATH
    ))

    run(MAIN_PATH, ALIGNER, DISTANCE_METRIC, correction=CORRECTION, exclude_stats= EXCLUDE_STATS)


if __name__ == "__main__":
    main()