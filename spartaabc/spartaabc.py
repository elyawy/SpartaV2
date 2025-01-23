
import random
from pathlib import Path
from dataclasses import dataclass

import simulate_data
import correction
import abc_inference
from abc_inference import IndelParams


@dataclass
class IndelModel:
    model: str
    length_distribution: str

def prepare_data_for_inference(data_path: Path, sequnce_type: str="AA", 
                               number_of_simulations: int=100000, number_of_correction_simulations: int=500,
                               seed: int=42, aligner: str="mafft", keep_correction_data: bool=False,
                               indel_model: IndelModel=IndelModel("sim", "zipf"),
                               verbose: int=0) -> None:

    # simulate_data.generate_summary_statistics(data_path, seed, number_of_simulations, 
    #                                           indel_model.length_distribution, indel_model.model)
    
    correction.compute_alignment_correction(data_path, seed, sequnce_type, number_of_correction_simulations,
                                            aligner, indel_model.length_distribution,
                                            indel_model.model, keep_correction_data)
    

def infer_indel_params(data_path: Path, aligner: str="MAFFT", 
                      distance_metric: str="mahal", top_cutoff: int=100) -> IndelParams:
    indel_params =  abc_inference.run(data_path, aligner, distance_metric, top_cutoff)

    print(indel_params)



def standard_inference(data_path: Path, seed=random.randint(1, 1e6)):
    prepare_data_for_inference(data_path=test_path, number_of_simulations=100_000,
                            number_of_correction_simulations=1000, seed=seed, keep_correction_data=True)
    
    prepare_data_for_inference(data_path=test_path, 
                            number_of_simulations=100_000, 
                            number_of_correction_simulations=1000, 
                            indel_model=IndelModel(model="rim", length_distribution="zipf"),
                            seed=seed+2, keep_correction_data=True)

    infer_indel_params(data_path=test_path)

test_path = Path("tests/pipeline_test")
# standard_inference(test_path)
prepare_data_for_inference(data_path=test_path, number_of_simulations=100,
                        number_of_correction_simulations=150, keep_correction_data=True)
