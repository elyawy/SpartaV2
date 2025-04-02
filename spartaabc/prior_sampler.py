import numpy as np
import random
import json
from typing import Dict

from msasim import sailfish as sf

from spartaabc.getting_priors import get_means

# Default truncation value, can be overridden in config
DEFAULT_TRUNCATION = 150

length_distribution_priors = {
    "zipf": {
        "insertion": sorted(get_means.final_priors["zipf"]),
        "deletion": sorted(get_means.final_priors["zipf"])
    },
    "geometric": {
        "insertion": sorted(get_means.final_priors["geometric"]),
        "deletion": sorted(get_means.final_priors["geometric"])
    },
    "poisson": {
        "insertion": sorted(get_means.final_priors["poisson"]),
        "deletion": sorted(get_means.final_priors["poisson"])
    }
}


def fast_zipf(a_param, truncation):
    harmonic_series = np.arange(1, truncation+1)
    harmonic_series = np.power(harmonic_series, -a_param)
    harmonic_sum = np.sum(harmonic_series)
    return harmonic_series / harmonic_sum


length_dist_mapper = {
    "zipf": sf.CustomDistribution,
    "poisson": sf.PoissonDistribution,
    "geometric": sf.GeometricDistribution
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
    def __init__(self, conf_file=None,
                 len_dist="zipf",
                 rate_priors=[[-4, -1], [-1, 1]],  # log
                 seq_lengths=[100, 500],
                 indel_model="sim",
                 seed=1,
                 truncation=DEFAULT_TRUNCATION):
        self.seed = seed
        random.seed(seed)
        self.indel_model = indel_model
        self.truncation = truncation

        # Set length distribution and indel model directly, not from config
        self.length_distribution = length_dist_mapper[len_dist]
        self.len_dist = len_dist
        self.indel_model = indel_model
        self.seq_lengths = seq_lengths
        # Default configuration - exclude length_distribution and indel_model
        self.config = {
            "sequence_length": {
                "method": "integer_uniform",
                "range": seq_lengths,
                "scale_factor": [0.8, 1.1]  # Used to adjust the range as in original code
            },
            "indel_rates": {
                "sum_rates": {
                    "method": "log_uniform",
                    "range": rate_priors[0]
                },
                "ratio_rates": {
                    "method": "log_uniform",
                    "range": rate_priors[1]
                }
            },
            "length_distribution_params": {
                "insertion": {
                    "method": "uniform",
                    "range": length_distribution_priors[len_dist]["insertion"]
                },
                "deletion": {
                    "method": "uniform",
                    "range": length_distribution_priors[len_dist]["deletion"]
                }
            },
            "truncation": truncation  # Add truncation parameter to config
        }

        # Load configuration from file if provided
        if conf_file:
            self._load_config(conf_file)
        
        # Initialize based on configuration
        self._initialize_from_config()

    def _load_config(self, conf_file: str) -> None:
        """Load configuration from JSON file"""
        try:
            with open(conf_file, 'r') as f:
                if conf_file.endswith('.json'):
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
        
        # Get truncation value
        self.truncation = self.config.get("truncation", DEFAULT_TRUNCATION)
        
        # Set sequence length prior - using constructor-provided seq_lengths
        # Do not override with config values
        scale = self.config["sequence_length"]["scale_factor"]
        self.sequence_length_prior = [
            int(self.seq_lengths[0] * scale[0]), 
            int(self.seq_lengths[1] * scale[1])
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
        
        self.ins_dist_sampler = SamplingMethod.get_sampler(
            self.config["length_distribution_params"]["insertion"]["method"]
        )
        self.ins_dist_range = self.config["length_distribution_params"]["insertion"]["range"]
        
        self.del_dist_sampler = SamplingMethod.get_sampler(
            self.config["length_distribution_params"]["deletion"]["method"]
        )
        self.del_dist_range = self.config["length_distribution_params"]["deletion"]["range"]

    def sample_root_length(self):
        while True:
            root_length = self.seq_length_sampler(*self.sequence_length_prior)
            yield root_length

    def sample_length_distributions(self):
        while True:
            x = self.ins_dist_sampler(*self.ins_dist_range)
            
            if self.indel_model == "sim":
                indel_length_dist = self.length_distribution(fast_zipf(x, self.truncation))
                indel_length_dist.p = x
                yield self.len_dist, indel_length_dist, indel_length_dist
            else:
                y = self.del_dist_sampler(*self.del_dist_range)
                indel_length_dist_insertion = self.length_distribution(fast_zipf(x, self.truncation))
                indel_length_dist_insertion.p = x
                indel_length_dist_deletion = self.length_distribution(fast_zipf(y, self.truncation))
                indel_length_dist_deletion.p = y
                yield self.len_dist, indel_length_dist_insertion, indel_length_dist_deletion

    def sample_rates(self):
        while True:
            sum_of_rates = self.sum_rates_sampler(*self.sum_rates_range)                
            ratio_of_rates = self.ratio_rates_sampler(*self.ratio_rates_range)
            
            
            if self.indel_model == "sim":
                # In the "sim" case, we set both insertion and deletion rates equal to sum_of_rates
                # This means the total mutation rate is actually 2*sum_of_rates
                # This maintains consistency with the original code's behavior
                yield (sum_of_rates, sum_of_rates)
            else:
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
            f"  Length Distribution: {self.len_dist}",
            f"  Indel Model: {self.indel_model}",
            f"  Truncation: {self.truncation}",
            f"  Sequence Length: {self.sequence_length_prior}",
            "  Indel Rates:",
            f"    Sum Rates: method={self.config['indel_rates']['sum_rates']['method']}, range={self.sum_rates_range}",
            f"    Ratio Rates: method={self.config['indel_rates']['ratio_rates']['method']}, range={self.ratio_rates_range}",
            "  Length Distribution Parameters:",
            f"    Insertion: method={self.config['length_distribution_params']['insertion']['method']}, range={self.ins_dist_range}",
            f"    Deletion: method={self.config['length_distribution_params']['deletion']['method']}, range={self.del_dist_range}"
        ]
        return "\n".join(representation)