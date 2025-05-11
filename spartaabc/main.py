import logging
import argparse
import subprocess
import sys
import time

from pathlib import Path

from spartaabc.utility import logger, setLogHandler

interpreter=sys.executable

def parse_args(arg_list: list[str] | None):
    _parser = argparse.ArgumentParser(allow_abbrev=False)
    _parser.add_argument('-i','--input', action='store',metavar="Input folder", type=str, required=True)
    # _parser.add_argument('-c','--config', action='store',metavar="Simulation config" , type=str, required=True)
    _parser.add_argument('-t','--type', action='store',metavar="Type of MSA NT/AA" , type=str, required=True)
    _parser.add_argument('-n','--numsim', action='store',metavar="Number of simulations" , type=int, required=True)
    _parser.add_argument('-nc','--numsim-correction', action='store',metavar="Number of correction simulations" , type=int, required=True)
    _parser.add_argument('-noc','--no-correction', action='store_false')

    _parser.add_argument('-s','--seed', action='store',metavar="Simulator seed" , type=int, required=False)
    _parser.add_argument('-a','--aligner', action='store',metavar="Alignment program to use" , type=str, default="mafft", required=False)

    _parser.add_argument('-k','--keep-stats', action='store_true')
    _parser.add_argument('-v','--verbose', action='store_true')


    args = _parser.parse_args()
    return args

def run_pipeline_parallel(MAIN_PATH: Path, SEED: int, SEQUENCE_TYPE: str,
                          NUM_SIMS: int, NUM_SIMS_CORRECTION: int, INDEL_MODELS: list[str],
                          ALIGNER: str, KEEP_STATS: bool):
    logging.basicConfig()

    CURRENT_SCRIPT_DIR = Path(__file__).parent
    print(CURRENT_SCRIPT_DIR)


    print()

    setLogHandler(MAIN_PATH, "w")
    logger.info("\n\tMAIN_PATH: {},\n\tSEED: {}, NUM_SIMS: {}, NUM_SIMS_CORRECTION: {}, SEQUENCE_TYPE: {}".format(
        MAIN_PATH, SEED, NUM_SIMS, NUM_SIMS_CORRECTION, SEQUENCE_TYPE
    ))


    processes = []
    for model in INDEL_MODELS:
        simulate_cmd = [interpreter, CURRENT_SCRIPT_DIR / "simulate_data.py",
                        "-i", str(MAIN_PATH), "-n", str(NUM_SIMS),
                        "-s", str(SEED), "-l", "zipf", "-m", f"{model}"]
        
    
        correction_cmd_sim = [interpreter, CURRENT_SCRIPT_DIR / "correction.py",
                              "-i", str(MAIN_PATH), "-n", str(NUM_SIMS_CORRECTION),
                              "-s", str(SEED+1), "-l", "zipf", "-m", f"{model}",
                              "-t", SEQUENCE_TYPE]
        SEED += 2

        processes.append(subprocess.Popen(simulate_cmd))
        processes.append(subprocess.Popen(correction_cmd_sim))

    exit_codes = [p.wait() for p in processes]
    
    abc_cmd = [interpreter, CURRENT_SCRIPT_DIR / "abc_inference.py", "-i", str(MAIN_PATH)]
    subprocess.run(abc_cmd)



def main(arg_list: list[str] | None = None):
    logging.basicConfig()

    CURRENT_SCRIPT_DIR = Path(__file__).parent
    print(CURRENT_SCRIPT_DIR)
    args = parse_args(arg_list)

    MAIN_PATH = Path(args.input).resolve()
    SEED = args.seed if args.seed else time.time_ns()
    SEQUENCE_TYPE = args.type
    NUM_SIMS = args.numsim
    NUM_SIMS_CORRECTION = args.numsim_correction
    CORRECTION = args.no_correction
    print(SEED)

    ALIGNER = args.aligner.upper()
    KEEP_STATS = args.keep_stats
    VERBOSE = args.verbose


    setLogHandler(MAIN_PATH, "w")
    logger.info("\n\tMAIN_PATH: {},\n\tSEED: {}, NUM_SIMS: {}, NUM_SIMS_CORRECTION: {}, SEQUENCE_TYPE: {}".format(
        MAIN_PATH, SEED, NUM_SIMS, NUM_SIMS_CORRECTION, SEQUENCE_TYPE
    ))

    INDEL_MODELS = ["sim", "rim"]

    processes = []
    for model in INDEL_MODELS:
        simulate_cmd = [interpreter, CURRENT_SCRIPT_DIR / "simulate_data.py",
                        "-i", str(MAIN_PATH), "-n", str(NUM_SIMS),
                        "-s", str(SEED), "-l", "zipf", "-m", f"{model}"]
    
        if not CORRECTION:
            SEED += 1
            processes.append(subprocess.Popen(simulate_cmd))
            continue
        
        correction_cmd_sim = [interpreter, CURRENT_SCRIPT_DIR / "correction.py",
                              "-i", str(MAIN_PATH), "-n", str(NUM_SIMS_CORRECTION),
                              "-s", str(SEED+1), "-l", "zipf", "-m", f"{model}",
                              "-t", SEQUENCE_TYPE, "-a", ALIGNER]
        SEED += 2

        processes.append(subprocess.Popen(simulate_cmd))
        processes.append(subprocess.Popen(correction_cmd_sim))

    exit_codes = [p.wait() for p in processes]
    
    abc_cmd = [interpreter, CURRENT_SCRIPT_DIR / "abc_inference.py",
               "-i", str(MAIN_PATH),"-a", ALIGNER
               ]
    if not CORRECTION:
        abc_cmd.append("-noc")

    subprocess.run(abc_cmd)




if __name__=="__main__":
    main()