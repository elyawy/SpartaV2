import os
import tempfile
import pickle
import argparse
import warnings
import copy
import logging
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn import linear_model, model_selection, exceptions
from sklearn.pipeline import Pipeline
from scipy.stats import pearsonr

from msasim import sailfish as sf
import msastats

from spartaabc.prior_sampler import PriorSampler, protocol_updater
from spartaabc.aligner_interface import Aligner
from spartaabc.raxml_parser import get_substitution_model
from spartaabc.utility import get_msa_path, get_tree_path, prepare_prior_sampler, setLogHandler
from spartaabc.utility import StandardMemoryScaler, logger
from spartaabc.utility import default_prior_config_path

warnings.simplefilter("ignore", category=exceptions.ConvergenceWarning)

def parse_args(arg_list: list[str] | None):
    _parser = argparse.ArgumentParser(allow_abbrev=False)
    _parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
    # _parser.add_argument('-c','--config', action='store',metavar="Simulation config" , type=str, required=True)
    _parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)
    _parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
    _parser.add_argument('-s','--seed', action='store',metavar="Simulator seed" , type=int, required=True)
    _parser.add_argument('-a','--aligner', action='store',metavar="Alignment program to use" , type=str, default="mafft", required=False)
    _parser.add_argument('-p','--prior', action='store',metavar="Prior config path" , type=str, required=False, default=default_prior_config_path)

    _parser.add_argument('-m','--model', action='store',metavar="Simulation config" , type=str, required=True)
    _parser.add_argument('-k','--keep-stats', action='store_true')
    _parser.add_argument('-v','--verbose', action='store_true')


    args = _parser.parse_args()
    return args



# prepare indelible control file for subsitutions:
def prepare_substitution_model(main_path: Path, sequence_type: str):

    substitution_model = get_substitution_model(main_path)# if sequence_type == "NT" else {}
    substitution_model["mode"] = "DNA" if sequence_type == "NT" else "PROTEIN"
    
    return substitution_model


def simulate_data(prior_sampler: PriorSampler, num_sims: int, tree_path: str, substitution_model: dict, seed: int):
    print(seed)
    sim_protocol = sf.SimProtocol(tree_path, seed=seed)
    simulator = sf.Simulator(sim_protocol,
                             simulation_type=sf.SIMULATION_TYPE[substitution_model["mode"]])

    simulated_msas = []
    sum_stats = []

    sim_params_correction = prior_sampler.sample(num_sims)

    simulator.set_replacement_model(model=substitution_model["submodel"],
                                    model_parameters=substitution_model.get("params", None),
                                    gamma_parameters_alpha=substitution_model.get("gamma_shape", 1.0),
                                    gamma_parameters_categories=substitution_model.get("gamma_cats", 1),
                                    invariant_sites_proportion=substitution_model.get("invariant_sites", 0.0))

    logger.info("Starting msa simulation")
    for idx,params in enumerate(sim_params_correction):
        root_length = params[0]
        insertion_rate, deletion_rate = params[1]
        lendist, insertion_length_dist, deletion_length_dist = params[2]

        numeric_params = [root_length, insertion_rate, deletion_rate, insertion_length_dist.p, deletion_length_dist.p]
        protocol_updater(sim_protocol, [root_length, insertion_rate, deletion_rate,
                         insertion_length_dist, deletion_length_dist])

        sim_msa = simulator()
        sim_stats = msastats.calculate_msa_stats(sim_msa.get_msa().splitlines()[1::2])
        # print(sim_stats)
        simulated_msas.append(sim_msa)
        sum_stats.append(numeric_params + sim_stats)
    logger.info(f"Done with {num_sims} msa simulations")

    return simulated_msas, sum_stats


def compute_realigned_stats(msas: list[sf.Msa], sum_stats: list[list[float]],
                            sequence_aligner: Aligner, tree_path: str):
    logger.info("Realigning MSAs and recomputing stats")
    realigned_sum_stats = []

    for msa in msas:
        sim_fasta_unaligned = msa.get_msa().replace("-","").encode()
        with tempfile.NamedTemporaryFile(suffix='.fasta') as tempf:
            tempf.write(sim_fasta_unaligned)
            tempf.seek(0)
            sequence_aligner.set_input_file(tempf.name, tree_file=tree_path)
            realigned_msa = sequence_aligner.get_realigned_msa()
            
        realigned_msa = [s[s.index("\n"):].replace("\n","") for s in realigned_msa.split(">")[1:]]
        realigned_stats = msastats.calculate_msa_stats(realigned_msa)
        realigned_sum_stats.append(realigned_stats)
    
    logger.info("Done recomputing stats")
    return realigned_sum_stats

def compute_regressors(true_stats: list[list[float]], corrected_stats: list[list[float]]):
    logger.info("Performing regression on all realigned stats")

    X = np.array(true_stats, dtype=float)
    Y = np.array(corrected_stats, dtype=float).T

    reg = linear_model.Lasso()
    parameters = {'alpha':np.logspace(-7,4,20)}
    clf_lassocv = model_selection.GridSearchCV(estimator = reg,
                                param_grid = parameters, cv=3,
                                scoring = 'neg_mean_squared_error')
    regression_pipline = Pipeline([("scaler", StandardMemoryScaler()),('regression', clf_lassocv)])
    regressors = []
    performance_metrics = []
    for y in Y:
        regression_pipline.fit(X, y)
        saved_estimator = copy.deepcopy(regression_pipline)
        regressors.append(saved_estimator)
        
        Y_pred = regression_pipline.predict(X)
        r_val, p_val = pearsonr(Y_pred,y)
        performance_metrics.append({
            'pearsonr': r_val,
            'p_val': p_val,
            'mean_test_score': np.min(np.sqrt(-clf_lassocv.cv_results_['mean_test_score']))
        })
    logger.info("Done with regression")
    return regressors, performance_metrics


def main(arg_list: list[str] | None = None):
    logging.basicConfig()
    args = parse_args(arg_list)

    MAIN_PATH = Path(args.input).resolve()
    SEED = args.seed
    SEQUENCE_TYPE = args.type
    NUM_SIMS = args.numsim
    ALIGNER = Aligner(args.aligner.upper())
    INDEL_MODEL = args.model
    KEEP_STATS = args.keep_stats
    PRIOR_PATH = Path(args.prior).resolve()

    VERBOSE = args.verbose

    setLogHandler(MAIN_PATH)
    logger.info("\n\tMAIN_PATH: {},\n\tSEED: {}, NUM_SIMS: {}, SEQUENCE_TYPE: {},\n\tINDEL_MODEL: {},\n\tALIGNER: {}".format(
        MAIN_PATH, SEED, NUM_SIMS, SEQUENCE_TYPE, INDEL_MODEL, ALIGNER._aligner_name
    ))


    TREE_PATH = get_tree_path(MAIN_PATH)
    MSA_PATH = get_msa_path(MAIN_PATH)

    prior_sampler = prepare_prior_sampler(MSA_PATH, INDEL_MODEL, SEED, PRIOR_PATH)
    logger.info("\nLoaded prior configuration from file {} is:\n\t{}".format(PRIOR_PATH, prior_sampler))
    LENGTH_DISTRIBUTION = prior_sampler.len_dist

    substitution_model = prepare_substitution_model(MAIN_PATH, SEQUENCE_TYPE)
    msas, sum_stats = simulate_data(prior_sampler=prior_sampler, num_sims=NUM_SIMS,
                                    tree_path=TREE_PATH, substitution_model=substitution_model,
                                    seed=SEED)
    realigned_sum_stats = compute_realigned_stats(msas, sum_stats, ALIGNER, TREE_PATH)
    regressors, regressors_performence = compute_regressors(true_stats=sum_stats, corrected_stats=realigned_sum_stats)

    full_correction_path = MAIN_PATH / f"{args.aligner}_correction"
    try:
        os.mkdir(full_correction_path)
    except:
        print("correction folder exists already")
    
    with open(full_correction_path / f'regressors_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}.pickle', 'wb') as f:
        pickle.dump(regressors, f)
    pd.DataFrame(regressors_performence).to_csv(
        full_correction_path / f'regression_performance_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}.csv')

    if KEEP_STATS:
        print("saving stats...")
        true_stats = pd.DataFrame(sum_stats)
        true_stats.columns = map(str, range(len(true_stats.columns)))
        true_stats.to_parquet(full_correction_path / "true_stats.parquet.gzip", compression='gzip', index=False)

        realigned_stats = pd.DataFrame(realigned_sum_stats)
        realigned_stats.columns = map(str, realigned_stats.columns)
        realigned_stats.to_parquet(full_correction_path / "realigned_stats.parquet.gzip", compression='gzip', index=False)

        infered_stats = np.array([regressor.predict(true_stats.values).T for regressor in regressors])
        infered_stats = pd.DataFrame(infered_stats.T, columns=map(str, range(27)))
        infered_stats.to_parquet(full_correction_path / "infered_stats.parquet.gzip", compression='gzip', index=False)



if __name__ == '__main__':
    main()