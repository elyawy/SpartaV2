import os, pathlib, tempfile, pickle, argparse, warnings, copy
from io import StringIO
import numpy as np
import pandas as pd
# from Bio import Phylo
from sklearn.base import BaseEstimator, TransformerMixin# define the transformer
from sklearn import linear_model, model_selection, exceptions
from sklearn.pipeline import Pipeline
warnings.simplefilter("ignore", category=exceptions.ConvergenceWarning)
from scipy.stats import pearsonr

from msasim import sailfish as sf
import msastats

from prior_sampler import SimConfig, protocol_updater
from aligner_interface import Aligner
from raxml_parser import get_substitution_model

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

_parser = argparse.ArgumentParser(allow_abbrev=False)
_parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
# _parser.add_argument('-c','--config', action='store',metavar="Simulation config" , type=str, required=True)
_parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)
_parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
# _parser.add_argument('-s','--seed', action='store',metavar="Simulation config" , type=int, required=False)
_parser.add_argument('-a','--aligner', action='store',metavar="Alignment program to use" , type=str, required=True)

_parser.add_argument('-l','--lengthdist', action='store',metavar="Simulation config" , type=str, required=True)
_parser.add_argument('-m','--model', action='store',metavar="Simulation config" , type=str, required=True)
_parser.add_argument('-k','--keep-stats', action='store_true')
_parser.add_argument('-v','--verbose', action='store_true')


args = _parser.parse_args()

MAIN_PATH = args.input
MODE = args.type
NUM_SIMS = args.numsim
ALIGNER = Aligner(args.aligner.upper())
LENGTH_DISTRIBUTION = args.lengthdist
INDEL_MODEL = args.model
KEEP_STATS = args.keep_stats
VERBOSE = args.verbose


MAIN_PATH = pathlib.Path(args.input).resolve()
TREE_PATH = None
MSA_PATH = None
if len( n := list(MAIN_PATH.glob("*.tree")) + list(MAIN_PATH.glob("*.newick"))) == 1:
    TREE_PATH = str(n[0])

if len( n := list(MAIN_PATH.glob("*.fasta"))) == 1:
    MSA_PATH = str(n[0])


if TREE_PATH is None or MSA_PATH is None:
    print("no fasta or tree file")
    exit()
if TREE_PATH is None or MSA_PATH is None:
    print("no fasta or tree file")
    exit()

min_length_index = msastats.stats_names().index("MSA_MIN_LEN")
max_length_index = msastats.stats_names().index("MSA_MAX_LEN")

    
empirical_stats = msastats.calculate_fasta_stats(MSA_PATH)
smallest_sequence_size = empirical_stats[min_length_index]
largest_sequence_size = empirical_stats[max_length_index]

seq_lengths_in_msa = [smallest_sequence_size, largest_sequence_size]

sim_config = SimConfig(seq_lengths=seq_lengths_in_msa,
                       len_dist=LENGTH_DISTRIBUTION,
                       indel_model=INDEL_MODEL)
# np.random.seed(int(sim_config.seed))


# prepare indelible control file for subsitutions:
substitution_model = get_substitution_model(MAIN_PATH) if MODE == "NT" else {}
substitution_model["mode"] = "DNA" if MODE == "NT" else "PROTEIN"

if substitution_model["mode"] == "PROTEIN":
    substitution_model["submodel"] = sf.MODEL_CODES.WAG,
    substitution_model["gamma_shape"] = 1.0
    substitution_model["gamma_cats"] = 4
else:
    substitution_model["params"] = substitution_model["freq"] + substitution_model["rates"]


def init_correction(sim_config, num_sims):
    sim_protocol = sf.SimProtocol(TREE_PATH)
    simulator = sf.Simulator(sim_protocol,
                             simulation_type=sf.SIMULATION_TYPE[substitution_model["mode"]])

    correction_list = []
    correction_list_sum_stats = []

    sim_params_correction = sim_config.get_random_sim(num_sims)
    simulator.set_replacement_model(model=substitution_model["submodel"][0],
                                    model_parameters=substitution_model.get("params", None),
                                    gamma_parameters_alpha=substitution_model.get("gamma_shape", 1.0),
                                    gamma_parameters_catergories=substitution_model.get("gamma_cats", 1))

    for idx,params in enumerate(sim_params_correction):
        numeric_params = [params[0],params[1], params[2], params[4].p, params[5].p]
        protocol_updater(sim_protocol, params)

        sim_msa = simulator()
        sim_stats = msastats.calculate_msa_stats(sim_msa.get_msa().splitlines()[1::2])
        # print(sim_stats)
        correction_list.append(sim_msa)
        correction_list_sum_stats.append(numeric_params + sim_stats)    

    return correction_list, correction_list_sum_stats
msas_list, sum_stats_list = init_correction(sim_config, NUM_SIMS)


def get_regressors(msa_list: list[sf.Msa], sum_stats_list):
    realigned_sum_stats = []
    for msa in msa_list:
        sim_fasta_unaligned = msa.get_msa().replace("-","").encode()

        with tempfile.NamedTemporaryFile(suffix='.fasta') as tempf:
            tempf.write(sim_fasta_unaligned)
            tempf.seek(0)
            ALIGNER.set_input_file(tempf.name, tree_file=TREE_PATH)
            realigned_msa = ALIGNER.get_realigned_msa()
        
        realigned_msa = [s[s.index("\n"):].replace("\n","") for s in realigned_msa.split(">")[1:]]
        realigned_stats = msastats.calculate_msa_stats(realigned_msa)
        realigned_sum_stats.append(realigned_stats)

    def compute_regressors(true_stats, corrected_stats):
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
        return regressors, performance_metrics

    regressors, performance = compute_regressors(sum_stats_list, realigned_sum_stats)        

    return regressors, performance, realigned_sum_stats
regressors, performance, realigned_sum_stats = get_regressors(msas_list, sum_stats_list)


if MAIN_PATH is None:
    for weight in [reg.best_model_.coef_ for reg in regressors]:
        print(",".join([str(w) for w in weight]))
else:
    full_correction_path = os.path.join(MAIN_PATH, f"{args.aligner}_correction")
    path_joiner = lambda x: os.path.join(full_correction_path, x)
    try:
        os.mkdir(full_correction_path)
    except:
        print("correction folder exists already")
    pickle.dump(regressors, open(path_joiner(f'regressors_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}'), 'wb'))
    pd.DataFrame(performance).to_csv(path_joiner(f'regression_performance_{LENGTH_DISTRIBUTION}_{INDEL_MODEL}.csv'))

if KEEP_STATS:
    print("saving stats...")
    true_stats = pd.DataFrame(sum_stats_list)
    true_stats.columns = map(str, range(len(true_stats.columns)))
    true_stats.to_parquet("true_stats.parquet.gzip", compression='gzip', index=False)

    realigned_stats = pd.DataFrame(realigned_sum_stats)
    realigned_stats.columns = map(str, realigned_stats.columns)
    realigned_stats.to_parquet("realigned_stats.parquet.gzip", compression='gzip', index=False)

    infered_stats = np.array([regressor.predict(true_stats.values).T for regressor in regressors])
    infered_stats = pd.DataFrame(infered_stats.T, columns=map(str, range(27)))
    infered_stats.to_parquet("infered_stats.parquet.gzip", compression='gzip', index=False)



