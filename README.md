# SpartaABC

A Python package for inferring indel (insertion-deletion) model parameters using Approximate Bayesian Computation (ABC).

## Description

`spartaabc` implements ABC methods to estimate parameters of indel models from observed data. This approach is particularly useful when dealing with complex evolutionary models where traditional likelihood-based methods may be computationally intractable.

## Installation

You can install `spartaabc` directly from PyPI:

```bash
pip install spartaabc
```

## Quick Start

Here's a basic example of how to use `spartaabc`:

```bash
sparta -i /path/to/data/dir -t AA -n 10000 -nc 100 -s 42
```

Arguments:
- `-i`: Input directory containing phylogenetic tree in newick format, and raxml ".bestModel" file.
- `-t`: Sequence type (AA for amino acid or NT for nucleotide)
- `-n`: Number of simulations (10000)
- `-nc`: Number of corrective simulations (100)
- `-s`: Random number seed (42)



## License

This project is licensed under the Academic Free License v3.0.

## Contact

- GitHub Issues: [Create an issue](https://github.com/elyawy/SpartaV2/issues)
- Email: elya.wygoda@gmail.com
