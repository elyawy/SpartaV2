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

def bias_correction(regressors, data: pd.DataFrame, regressor_scores: list[float], r_threshold=0.8):
    data = data.to_numpy()

    kept_stats = []
    infered_data = []
    for idx, regressor in enumerate(regressors):
        if regressor_scores[idx] > r_threshold:
            kept_stats.append(idx)
            infered_data.append(regressor.predict(data).T)

    temp_data = np.array(infered_data)
    temp_data = pd.DataFrame(temp_data.T, columns=[SUMSTATS_LIST[i] for i in kept_stats])
    return temp_data, kept_stats

def run(main_path: Path, aligner: str, distance_metric: str="mahal", correction=True, top_cutoff: int=100) -> IndelParams:


    MSA_PATH = get_msa_path(main_path)

    empirical_stats = msastats.calculate_fasta_stats(MSA_PATH)

    stats_data = load_data(main_path)
    if correction:
        regressors = load_correction_regressors(main_path, aligner)
        regressor_scores = load_correction_regressor_scores(main_path, aligner)
    else:
        regressors = {}

    params_data = []
    full_stats_data = []
    kept_statistics = range(len(SUMSTATS_LIST))
    for model in  stats_data.keys():
        current_regressors = regressors.get(model, None)
        params_data.append(stats_data[model][PARAMS_LIST])

        if current_regressors is not None:
            temp_df, kept_statistics = bias_correction(current_regressors, stats_data[model], regressor_scores)
            full_stats_data.append(temp_df)

    empirical_stats = [empirical_stats[i] for i in kept_statistics]

    params_data = pd.concat(params_data)
    if correction:
        full_stats_data = pd.concat(full_stats_data)
    else:
        full_stats_data = [val[SUMSTATS_LIST] for key, val in stats_data.items()]
        full_stats_data = pd.concat(full_stats_data)
    calculated_distances = None

    if distance_metric == "mahal":
        cov = np.cov(full_stats_data.T)
        cov = cov + np.eye(len(cov))*1e-4
        inv_covmat = np.linalg.inv(cov)
        u_minus_v = empirical_stats-full_stats_data
        left = np.dot(u_minus_v, inv_covmat)
        calculated_distances = np.sqrt(np.sum(u_minus_v*left, axis=1))
    if distance_metric == "euclid":
        weights = 1/(full_stats_data.std(axis=0) + 0.001)
        calculated_distances = np.sum(weights*(full_stats_data - empirical_stats)**2, axis=1)
    
    full_stats_data["distances"] = calculated_distances
    full_stats_data[PARAMS_LIST] = params_data

    top_stats = full_stats_data.nsmallest(top_cutoff, "distances")

    top_stats[["distances"] + PARAMS_LIST].to_csv(main_path / "top_params.csv", index=False)

    full_sim_data = full_stats_data[full_stats_data["insertion_rate"] == full_stats_data["deletion_rate"]]
    top_sim_data = full_sim_data.nsmallest(top_cutoff, "distances")
    top_sim_data[["distances"] + PARAMS_LIST].to_csv(main_path / "top_params_sim.csv", index=False)

    full_rim_data = full_stats_data[full_stats_data["insertion_rate"] != full_stats_data["deletion_rate"]]
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

    setLogHandler(MAIN_PATH)
    logger.info("\n\tMAIN_PATH: {}".format(
        MAIN_PATH
    ))

    run(MAIN_PATH, ALIGNER, DISTANCE_METRIC, correction=CORRECTION)


if __name__ == "__main__":
    main()