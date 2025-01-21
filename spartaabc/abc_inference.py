import argparse, pickle
import numpy as np
import pandas as pd
from dataclasses import dataclass

from utility import *

@dataclass
class IndelParams:
    insertion_rate: float
    deletion_rate: float
    insertion_length_parameter: float
    deletion_length_parameter: float
    length_distribution: str
    indel_model: str


def parse_args(arg_list: list[str] | None):
    _parser = argparse.ArgumentParser(allow_abbrev=False)
    _parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
    _parser.add_argument('-a','--aligner', action='store',metavar="Aligner", type=str,default="mafft" , required=False)
    _parser.add_argument('-d','--distance', action='store',metavar="Distance metric", type=str, default="mahal", required=False)


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

def bias_correction(regressors, data: pd.DataFrame):
    data = data.to_numpy()
    temp_data = np.array([regressor.predict(data).T for regressor in regressors])
    temp_data = pd.DataFrame(temp_data.T, columns=SUMSTATS_LIST)
    return temp_data

def run(main_path: Path, aligner: str, distance_metric: str="mahal", top_cutoff: int=100) -> IndelParams:


    MSA_PATH = get_msa_path(main_path)

    empirical_stats = msastats.calculate_fasta_stats(MSA_PATH)

    stats_data = load_data(main_path)
    regressors = load_correction_regressors(main_path, aligner)

    params_data = []
    full_stats_data = []
    for model in  stats_data.keys():
        current_regressors = regressors.get(model, None)
        params_data.append(stats_data[model][PARAMS_LIST])

        if current_regressors is not None:
            temp_df = bias_correction(current_regressors, stats_data[model])
            full_stats_data.append(temp_df)
    
    params_data = pd.concat(params_data)
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

    if len(top_stats[top_stats["insertion_rate"] == top_stats["deletion_rate"]]) > top_cutoff // 2:
        print("SIM")
        full_sim_data = full_stats_data[full_stats_data["insertion_rate"] == full_stats_data["deletion_rate"]]
        top_sim_data = full_sim_data.nsmallest(top_cutoff, "distances")
        R_ID = top_sim_data["insertion_rate"].mean()
        A_ID = top_sim_data["length_param_insertion"].mean()
        return IndelParams(R_ID, R_ID,
                           A_ID, A_ID,
                           length_distribution="zipf",
                           indel_model="SIM")

    else:
        full_rim_data = full_stats_data[full_stats_data["insertion_rate"] != full_stats_data["deletion_rate"]]
        top_rim_data = full_rim_data.nsmallest(top_cutoff, "distances")
        R_I = top_rim_data["insertion_rate"].mean()
        R_D = top_rim_data["deletion_rate"].mean()
        A_I = top_rim_data["length_param_insertion"].mean()
        A_D = top_rim_data["length_param_deletion"].mean()

        return IndelParams(R_I, R_D,
                           A_I, A_D,
                           length_distribution="zipf",
                           indel_model="RIM")



def main(arg_list: list[str] | None = None):
    logging.basicConfig()

    args = parse_args(arg_list)


    MAIN_PATH = Path(args.input).resolve()
    ALIGNER = args.aligner
    DISTANCE_METRIC = args.distance

    setLogHandler(MAIN_PATH)
    logger.info("\n\tMAIN_PATH: {}".format(
        MAIN_PATH
    ))

    run(MAIN_PATH, ALIGNER, DISTANCE_METRIC)


if __name__ == "__main__":
    main()