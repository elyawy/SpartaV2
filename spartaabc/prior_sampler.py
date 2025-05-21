import numpy as np
import random
import json
from pathlib import Path
from typing import Dict

from msasim import sailfish as sf


def fast_zipf(a_param, truncation):
    harmonic_series = np.arange(1, truncation+1)
    harmonic_series = np.power(harmonic_series, -a_param)
    harmonic_sum = np.sum(harmonic_series)
    return harmonic_series / harmonic_sum

def fast_geo(p_param, truncation):
    geometric_series = np.arange(1, truncation+1)
    geometric_series = np.power((1-p_param), (geometric_series-1))* p_param
    geometric_sum = np.sum(geometric_series)
    return geometric_series / geometric_sum

len_dist_mapper = {
    "zipf": fast_zipf,
    "geometric": fast_geo
}




def protocol_updater(protocol: sf.SimProtocol, params: list) -> None:
    protocol.set_sequence_size(params[0])
    protocol.set_insertion_rates(insertion_rate=params[1])
    protocol.set_deletion_rates(deletion_rate=params[2])
    protocol.set_insertion_length_distributions(insertion_dist=params[3])
    protocol.set_deletion_length_distributions(deletion_dist=params[4])


class SamplingMethod:
    """Class to handle different sampling methods for parameters"""
    @staticmethod
    def uniform(range_min: float, range_max: float) -> float:
        """Sample uniformly from [range_min, range_max]"""
        return random.uniform(range_min, range_max)
    
    @staticmethod
    def log_uniform(range_min: float, range_max: float) -> float:
        """Sample log-uniformly, returns 10^uniform([range_min, range_max])"""
        return 10 ** random.uniform(range_min, range_max)
    
    @staticmethod
    def integer_uniform(range_min: int, range_max: int) -> int:
        """Sample integer uniformly from [range_min, range_max]"""
        return random.randint(range_min, range_max)
    
    @staticmethod
    def get_sampler(method: str):
        """Return the appropriate sampling method function"""
        samplers = {
            "uniform": SamplingMethod.uniform,
            "log_uniform": SamplingMethod.log_uniform,
            "integer_uniform": SamplingMethod.integer_uniform
        }
        if method not in samplers:
            raise ValueError(f"Unknown sampling method: {method}. Available methods: {list(samplers.keys())}")
        return samplers[method]


class PriorSampler:
    def __init__(self, conf_file: Path=None,
                 len_dist="zipf",
                 rate_priors=[[0.0, 0.05], [-1, 1]],
                 length_distribution_priors=[1.001,2.0],
                 truncation=50,
                 seq_lengths=[100, 500],
                 indel_model="sim",
                 seed=1):
        self.seed = seed
        random.seed(seed)

        # Set length distribution and indel model directly, not from config
        self.length_distribution = sf.CustomDistribution
        self.indel_model = indel_model
        # Default configuration - exclude indel_model (sim/rim)
        self.config = {
            "sequence_length": {
                "method": "integer_uniform",
                "range": seq_lengths,
                "scale_factor": [0.8, 1.1]  # Used to adjust the range as in original code
            },
            "indel_rates": {
                "sum_rates": {
                    "method": "uniform",
                    "range": rate_priors[0]
                },
                "ratio_rates": {
                    "method": "log_uniform",
                    "range": rate_priors[1]
                }
            },
            "length_distribution_params": {
                    "distribution": len_dist,
                    "method": "uniform",
                    "range": length_distribution_priors,
                    "truncation": truncation
            },
        }

        # Load configuration from file if provided
        if conf_file:
            self._load_config(conf_file)
        
        # Initialize based on configuration
        self._initialize_from_config()

    def _load_config(self, conf_file: Path) -> None:
        """Load configuration from JSON file"""
        try:
            with open(conf_file, 'r') as f:
                if conf_file.suffix == '.json':
                    loaded_config = json.load(f)
                else:
                    raise ValueError("Config file must be JSON (.json)")
                
                # Update config with loaded values, keeping defaults for missing fields
                self._update_config_recursive(self.config, loaded_config)
        except Exception as e:
            print(f"Error loading configuration file: {e}")
            print("Using default configuration")

    def _update_config_recursive(self, default_dict: Dict, update_dict: Dict) -> None:
        """Recursively update default dictionary with values from update dictionary"""
        for key, value in update_dict.items():
            if key in default_dict and isinstance(default_dict[key], dict) and isinstance(value, dict):
                self._update_config_recursive(default_dict[key], value)
            else:
                default_dict[key] = value

    def _initialize_from_config(self) -> None:
        """Initialize sampler properties from configuration"""
        # Note: length_distribution and indel_model are already set in __init__
        # and are not read from the config file
        
        
        # Set sequence length prior - using constructor-provided seq_lengths
        # Do not override with config values
        scale = self.config["sequence_length"]["scale_factor"]
        sequence_length_range = self.config["sequence_length"]["range"]
        self.sequence_length_prior = [
            int(sequence_length_range[0] * scale[0]), 
            int(sequence_length_range[1] * scale[1])
        ]
        
        # Get samplers
        self.seq_length_sampler = SamplingMethod.get_sampler(
            self.config["sequence_length"]["method"]
        )
        
        self.sum_rates_sampler = SamplingMethod.get_sampler(
            self.config["indel_rates"]["sum_rates"]["method"]
        )
        self.sum_rates_range = self.config["indel_rates"]["sum_rates"]["range"]
        
        self.ratio_rates_sampler = SamplingMethod.get_sampler(
            self.config["indel_rates"]["ratio_rates"]["method"]
        )
        self.ratio_rates_range = self.config["indel_rates"]["ratio_rates"]["range"]
        
        self.len_dist = self.config["length_distribution_params"]["distribution"]

        self.length_param_sampler = SamplingMethod.get_sampler(self.config["length_distribution_params"]["method"])
        self.length_param_prior = self.config["length_distribution_params"]["range"]
        self.truncation = self.config["length_distribution_params"]["truncation"]


    def sample_root_length(self):
        while True:
            root_length = self.seq_length_sampler(*self.sequence_length_prior)
            yield root_length

    def sample_length_distributions(self):
        while True:
            x = self.length_param_sampler(*self.length_param_prior)
            
            if self.indel_model == "sim":
                probabilities = len_dist_mapper[self.len_dist](x, self.truncation)
                indel_length_dist = self.length_distribution(probabilities)
                indel_length_dist.p = x
                yield self.len_dist, indel_length_dist, indel_length_dist
            else:
                y = self.length_param_sampler(*self.length_param_prior)
                probabilities = len_dist_mapper[self.len_dist](x, self.truncation)
                indel_length_dist_insertion = self.length_distribution(probabilities)
                indel_length_dist_insertion.p = x
                probabilities = len_dist_mapper[self.len_dist](y, self.truncation)
                indel_length_dist_deletion = self.length_distribution(probabilities)
                indel_length_dist_deletion.p = y
                yield self.len_dist, indel_length_dist_insertion, indel_length_dist_deletion

    def sample_rates(self):
        while True:
            sum_of_rates = self.sum_rates_sampler(*self.sum_rates_range)                
            ratio_of_rates = self.ratio_rates_sampler(*self.ratio_rates_range)
            
            if self.indel_model == "sim":
                ratio_of_rates = 1.0
            # The ratio_of_rates is insertion_rate / deletion_rate
            # So to get the individual rates from the sum:
            # insertion_rate = (ratio_of_rates * deletion_rate)
            # insertion_rate + deletion_rate = sum_of_rates
            # Substituting: 
            # (ratio_of_rates * deletion_rate) + deletion_rate = sum_of_rates
            # deletion_rate * (ratio_of_rates + 1) = sum_of_rates
            deletion_rate = sum_of_rates / (ratio_of_rates + 1)
            insertion_rate = ratio_of_rates * deletion_rate
            
            # Verify the sum constraint is maintained
            assert abs((insertion_rate + deletion_rate) - sum_of_rates) < 1e-10, "Sum constraint violated"
            
            yield (insertion_rate, deletion_rate)

    def sample(self, n=1):
        root_length = self.sample_root_length()
        indel_rates = self.sample_rates()
        length_dists = self.sample_length_distributions()
        params_sample = []
        for params in zip(root_length, indel_rates, length_dists):
            if n == 0:
                break
            params_sample.append(params)
            n = n - 1
        return params_sample
        
    def __repr__(self):
        """Provide a string representation of the PriorSampler object."""
        representation = [
            f"PriorSampler(seed={self.seed})",
            f"  Indel Model: {self.indel_model}",
            f"  Sequence Length: {self.sequence_length_prior}",
            "  Indel Rates:",
            f"    Sum Rates: method={self.config['indel_rates']['sum_rates']['method']}, range={self.sum_rates_range}",
            f"    Ratio Rates: method={self.config['indel_rates']['ratio_rates']['method']}, range={self.ratio_rates_range}",
            "  Length Distribution Parameters:",
            f"    Length Distribution: {self.len_dist}",
            f"    Truncation: {self.truncation}",
            f"    Length Parameter: method={self.config['length_distribution_params']['method']}, range={self.length_param_prior}",
        ]
        return "\n".join(representation)