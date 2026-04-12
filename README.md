# SpartaABC

![PyPI version](https://img.shields.io/pypi/v/spartaabc)
![License](https://img.shields.io/badge/license-AFL--3.0-blue)

SpartaABC is a tool for inferring indel model parameters from multiple sequence alignments using Approximate Bayesian Computation (ABC).

If you encounter any problems, please open an [issue](https://github.com/elyawy/SpartaV2/issues).

---

## Table of Contents

- [What is SpartaABC?](#what-is-spartaabc)
- [Installation](#installation)
- [Input Files](#input-files)
- [Available Commands](#available-commands)
- [Simple Example (No Correction)](#simple-example-no-correction)
- [Advanced Inference with Correction](#advanced-inference-with-correction)
- [Using the sparta Shortcut](#using-the-sparta-shortcut)
- [Citations](#citations)
- [License](#license)

---

## What is SpartaABC?

When sequences evolve, they don't only accumulate substitutions — they also undergo **insertions** and **deletions** of sequence segments, collectively called **indels**. An **indel model** describes how these events occur: how frequently they happen, and what lengths the inserted or deleted segments tend to be.

Estimating indel model parameters from a real multiple sequence alignment (MSA) is not straightforward. The likelihood of an indel model given an alignment is generally intractable to compute directly, which rules out standard maximum-likelihood approaches.

SpartaABC addresses this using **Approximate Bayesian Computation (ABC)**. Instead of computing the likelihood directly, it:

1. Simulates a large number of alignments under different parameter values sampled from a prior distribution.
2. Computes summary statistics on both the simulated and the observed (empirical) alignments.
3. Selects the simulated parameter values whose summary statistics are closest to the empirical ones — these form the approximate posterior, from which point estimates are derived.

---

## Installation

SpartaABC requires Python 3.9 or later. Install it from PyPI:

```bash
pip install spartaabc
```

---

## Input Files

Every SpartaABC analysis requires the following files, placed in a single directory:

**1. A multiple sequence alignment (`.fasta`)**
Your empirical MSA in FASTA format. Sequences should already be aligned.

**2. A phylogenetic tree (`.tree` or `.newick`)**
A tree in Newick format with branch lengths. Branch lengths determine the evolutionary distances used during simulation.

**Optionally**, for the correction step (described later):

**3. A substitution model file (`.bestModel`)**
Describes the nucleotide or amino acid substitution model. See the [Advanced Inference](#advanced-inference-with-correction) section for details.

> SpartaABC expects exactly one file of each type in the input directory. If multiple files of the same type are found, it will use the first one and issue a warning.

---

## Available Commands

| Command | Description |
|---------|-------------|
| `sparta-simulate` | Simulate alignments under a specified prior |
| `sparta-correction` | Compute the correction regressors (requires a `.bestModel` file) |
| `sparta-abc` | Run the final ABC inference |
| `sparta` | Shortcut that runs all steps in sequence (see [Using the sparta Shortcut](#using-the-sparta-shortcut)) |

---

## Simple Example (No Correction)

For a first analysis, no substitution model is needed. We will use a small protein alignment of five species.

**`alignment.fasta`**
```
>SpeciesA
MCNCLLQLRELYKLCNERNITRDDVLELIDPLIKTRCFAYSLVVLANANPIALSILPRKILINGEPLLLEYGNIYGKDFLYRPSLQVILEEEELN
>SpeciesB
MCNCLYQVKALVEYSKTH--GKTDVLELLDPLVKTRCFAYTLVVCINANPVAFSILPRKLLINGEPLLIEYGNVYGKDFIIRPSLQVILEEEC--
>SpeciesC
MCNCIKQVAALVQHCKATNIHPSDVLELNDPLVAVKCLAYTLVLVTNADPVAFSILPRKILINGEPLLIEHGNVYGKDFLVRPSLQVILEEEVTD
>SpeciesD
MCNCHLQLRDLYRLCNKRHIRREDVPELIDPLVKTRCFAYSLVVLANANPIAFSILPRKILINGEPLLLEYGSIYGKDFIIRPSLQVILEDELN-
>SpeciesE
MCKCAQYIKVFNTN--IHGHSNTSLLTIDDHMLKYKCFAFALANLLVDNPIAGALLPRKMLINGKPILIEYGKIKADKFLINPPSKVEFYDD---
```

**`tree.newick`**
```
(SpeciesE:1.12098,(SpeciesB:2.1089e-06,(SpeciesC:0.328448,(SpeciesA:0.107594,SpeciesD:0.137885):0.333936):0.12119):1.12098);
```

Place both files in a directory, for example `my_analysis/`.

### Step 1: Simulate

SpartaABC supports two indel models — **SIM** and **RIM** (described below). To enable model selection, run a simulation for each:

```bash
sparta-simulate -i my_analysis/ -n 10000 -m sim
sparta-simulate -i my_analysis/ -n 10000 -m rim
```

| Flag | Description |
|------|-------------|
| `-i` | Path to the input directory |
| `-n` | Number of simulations to run |
| `-m` | Indel model to simulate under: `sim` or `rim` |

### Step 2: ABC Inference

Once both simulations are complete, run the inference:

```bash
sparta-abc -i my_analysis/ -noc
```

The `-noc` flag tells SpartaABC to skip the correction step.

### Output

Results are written to the input directory. The key output file is `model_params.txt`, which contains the inferred indel model. SpartaABC selects between two models:

#### SIM — Simple Indel Model

Insertions and deletions occur at the same rate and follow the same length distribution.

```
Model: SIM
Root_length: 95
R_ID: 0.0182
A_ID: 1.74
```

| Parameter | Meaning |
|-----------|---------|
| `Root_length` | Estimated length of the ancestral (root) sequence, in residues |
| `R_ID` | Shared indel rate per site per unit branch length |
| `A_ID` | Length distribution parameter — higher values mean shorter indels are more probable |

#### RIM — Rich Indel Model

Insertions and deletions have independent rates and length distributions.

```
Model: RIM
Root_length: 95
R_I: 0.0140
R_D: 0.0230
A_I: 1.61
A_D: 1.88
```

| Parameter | Meaning |
|-----------|---------|
| `Root_length` | Estimated length of the ancestral (root) sequence, in residues |
| `R_I` | Insertion rate per site per unit branch length |
| `R_D` | Deletion rate per site per unit branch length |
| `A_I` | Length distribution parameter for insertions |
| `A_D` | Length distribution parameter for deletions |

SpartaABC selects SIM if the majority of the top-scoring simulations have equal insertion and deletion rates, and RIM otherwise.

---

## Advanced Inference with Correction

### Why Correction is Needed

The empirical MSA was produced by an alignment program, while simulated alignments come directly from the simulator with perfect gap information. This mismatch introduces a systematic **aligner bias**: statistics from simulated data are not directly comparable to those from the empirical data.

The correction step re-aligns the simulated sequences with MAFFT, then trains a regression model to map uncorrected simulated statistics to corrected ones, leading to better-calibrated inference.

For best results, **re-align your empirical MSA using MAFFT's G-INSi mode** before running SpartaABC with correction, so that the empirical and simulated statistics are computed under the same alignment strategy.

### The `.bestModel` File

The correction step simulates full sequence evolution including substitutions, so SpartaABC needs to know the substitution model. This is provided as a `.bestModel` file in the input directory.

The format follows the substitution model notation used by **RAxML-NG** and **IQ-TREE**. For example:

**Nucleotide:**
```
GTR{1.2/3.4/0.8/1.1/2.5}+F{0.31/0.20/0.25/0.24}+G4{0.85}
```

**Protein:**
```
LG+G4{0.92}
```

> **Note:** SpartaABC supports a subset of the model syntax used by these tools. Since there is no single accepted format standard across phylogenetics software, some model specifications may use syntax variants that are not currently supported. If you encounter a parsing error, please open an [issue](https://github.com/elyawy/SpartaV2/issues).

### Running with Correction

**`alignment.fasta`**
```
>SpeciesA
GTATGTTAGACCATTTAGAACGAGTCAGGATCCACTGATTTATTCCGTAAGCG----ACATTATATTAACCAACCCTTTATA---G
>SpeciesB
GTATGTTGGATTACGCAAAGCAATTCAGGATCTACTGATCCAATC---AAGCTTGTGATTTTTTATTAACAAATTCTCTATAATAG
>SpeciesC
GTATGTTAAATTATGCAAAACAACTCAGGATCTTCTGGTCCAGTC---AAGCGTGTAATTT-TTATTAACAAGTTCTTTATAATAG
>SpeciesD
GTATGTTAGATTGTGCAAGAAAATTCAGATTCTACTTTCTCAATT---AAGCGTGTAATATTTTATTAACAAATCCTTTTTA---G
```

**`tree.newick`**
```
(SpeciesA:0.255377,(SpeciesC:0.081007,SpeciesB:0.084970):0.030395,SpeciesD:0.178788);
```

**`model.bestModel`**
```
GTR{1.0/4.0/1.0/1.0/4.0}+F{0.25/0.25/0.25/0.25}+G4{1.0}
```

**Step 1: Simulate**

```bash
sparta-simulate -i my_analysis/ -n 10000 -m sim
sparta-simulate -i my_analysis/ -n 10000 -m rim
```

**Step 2: Compute Correction Regressors**

The `-t` flag is required here, as the correction step needs to know whether to use a DNA or protein substitution model:

```bash
sparta-correction -i my_analysis/ -t NT -n 1000 -m sim
sparta-correction -i my_analysis/ -t NT -n 1000 -m rim
```

| Flag | Description |
|------|-------------|
| `-t` | Sequence type: `NT` for nucleotides, `AA` for amino acids |
| `-n` | Number of simulations for correction (200 simulations is typical) |
| `-m` | Indel model: `sim` or `rim` |

**Step 3: ABC Inference**

```bash
sparta-abc -i my_analysis/
```

---

## Using the `sparta` Shortcut

The `sparta` command runs all of the above steps in one go, spawning parallel processes. It requires at least **2 cores** without correction, or **4 cores** with correction:

```bash
# Without correction (2 cores recommended):
sparta -i my_analysis/ -t NT -n 10000

# With correction (4 cores recommended):
sparta -i my_analysis/ -t NT -n 10000 -nc 1000
```

| Flag | Description |
|------|-------------|
| `-i` | Path to the input directory |
| `-t` | Sequence type: `NT` or `AA` |
| `-n` | Number of simulations |
| `-nc` | Number of correction simulations; triggers the correction pipeline when provided |

---

## Citations

---

## License

This project is licensed under the Academic Free License v3.0.