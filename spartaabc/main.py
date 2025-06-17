import logging
import argparse
import subprocess
import sys
import time
from pathlib import Path

from spartaabc.utility import logger, setLogHandler, default_prior_config_path, validate_input_directory, parse_model_params



interpreter=sys.executable

def parse_args(arg_list: list[str] | None):
    _parser = argparse.ArgumentParser(allow_abbrev=False)
    _parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
    # _parser.add_argument('-c','--config', action='store',metavar="Simulation config" , type=str, required=True)
    _parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)
    _parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
    _parser.add_argument('-nc','--numsim-correction', action='store',metavar="Number of correction simulations", default=0 , type=int, required=False)

    _parser.add_argument('-s','--seed', action='store',metavar="Simulator seed" , type=int, required=False)
    _parser.add_argument('-p','--prior', action='store',metavar="Prior config path" , type=str, required=False, default=default_prior_config_path)

    _parser.add_argument('-a','--aligner', action='store',metavar="Alignment program to use" , type=str, default="mafft", required=False)

    _parser.add_argument('-k','--keep-stats', action='store_true')
    # _parser.add_argument('-v','--verbose', action='store_true')


    args = _parser.parse_args()
    return args

def format_results_summary(results_path: Path, sequence_type: str) -> str:
    """Create a human-readable results summary"""
    summary = []
    summary.append("🎯 SpartaABC Analysis Results")
    summary.append("=" * 40)
    
    # Read model parameters
    model_params = parse_model_params(results_path / "model_params.txt")
    
    summary.append(f"📊 Selected Model: {model_params['indel_model']}")
    summary.append(f"🧬 Root Sequence Length: {model_params['root_length']}" + ("bp" if sequence_type == "NT" else "aa"))
    
    if model_params['indel_model'] == "SIM":
        summary.append(f"🔄 Indel Rate: {model_params['insertion_rate']:.4f}")
        summary.append(f"📏 Indel Length Parameter: {model_params['length_param_insertion']:.4f}")
    else:
        summary.append(f"➕ Insertion Rate: {model_params['insertion_rate']:.4f}")
        summary.append(f"📏 Insertion Length Parameter: {model_params['length_param_insertion']:.4f}")

        summary.append(f"➖ Deletion Rate: {model_params['deletion_rate']:.4f}")
        summary.append(f"📏 Deletion Length Parameter: {model_params['length_param_deletion']:.4f}")

    
    # Add confidence intervals if available
    # summary.append("\n📈 Parameter Confidence:")
    # ... add confidence interval calculations ...
    
    return "\n".join(summary)

def main(arg_list: list[str] | None = None):
    logging.basicConfig()

    CURRENT_SCRIPT_DIR = Path(__file__).parent
    args = parse_args(arg_list)

    MAIN_PATH = Path(args.input).resolve()
    validation_status = validate_input_directory(MAIN_PATH)


    SEED = args.seed if args.seed else time.time_ns()
    SEQUENCE_TYPE = args.type
    NUM_SIMS = args.numsim

    NUM_SIMS_CORRECTION = args.numsim_correction
    CORRECTION = False
    if NUM_SIMS_CORRECTION > 0:
        CORRECTION = True
    
    PRIOR_PATH = args.prior

    ALIGNER = args.aligner.upper()
    KEEP_STATS = args.keep_stats
    # VERBOSE = args.verbose


    setLogHandler(MAIN_PATH, "w")
    logger.info("\n\tMAIN_PATH: {},\n\tSEED: {}, NUM_SIMS: {}, NUM_SIMS_CORRECTION: {}, SEQUENCE_TYPE: {}, PRIOR: {}".format(
        MAIN_PATH, SEED, NUM_SIMS, NUM_SIMS_CORRECTION, SEQUENCE_TYPE, PRIOR_PATH
    ))
    logger.info(validation_status)
    INDEL_MODELS = ["sim", "rim"]

    processes = []
    for model in INDEL_MODELS:
        simulate_cmd = [interpreter, CURRENT_SCRIPT_DIR / "simulate_data.py",
                        "-i", str(MAIN_PATH), "-n", str(NUM_SIMS),
                        "-s", str(SEED), "-m", f"{model}", "-p", PRIOR_PATH]
    
        if not CORRECTION:
            SEED += 1
            processes.append(subprocess.Popen(simulate_cmd))
            continue
        
        correction_cmd_sim = [interpreter, CURRENT_SCRIPT_DIR / "correction.py",
                              "-i", str(MAIN_PATH), "-n", str(NUM_SIMS_CORRECTION),
                              "-s", str(SEED+1), "-m", f"{model}",
                              "-t", SEQUENCE_TYPE, "-a", ALIGNER, "-p", PRIOR_PATH]
        if KEEP_STATS:
            correction_cmd_sim += ["-k"]
        
        SEED += 2
        processes.append(subprocess.Popen(simulate_cmd))
        processes.append(subprocess.Popen(correction_cmd_sim))

    exit_codes = [p.wait() for p in processes]
    
    abc_cmd = [interpreter, CURRENT_SCRIPT_DIR / "abc_inference.py",
               "-i", str(MAIN_PATH),"-a", ALIGNER, 
               ]
    if not CORRECTION:
        abc_cmd.append("-noc")

    subprocess.run(abc_cmd)

    print(format_results_summary(MAIN_PATH, SEQUENCE_TYPE))




if __name__=="__main__":
    main()