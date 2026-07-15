# Probing Synthetic Caroli-de Gennes-Matricon States Through Critical Current in Full-Shell Nanowire Josephson Junctions
[![DOI](https://zenodo.org/badge/1016588588.svg)](https://doi.org/10.5281/zenodo.21381725)
![arXiv](https://img.shields.io/badge/arXiv-soon-b31b1b.svg)
[![Julia v1.12+](https://img.shields.io/badge/Julia-v1.12+-blue.svg)](https://julialang.org/)
[![Quantica badge](https://raw.githubusercontent.com/pablosanjose/Quantica.jl/master/docs/src/assets/badge.svg)](https://github.com/pablosanjose/Quantica.jl)
[![License](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](LICENSE.md)

This repository contains the code, experimental data and numerical results used to generate all figures for the preprint "Probing Synthetic Caroli-de Gennes-Matricon States Through Critical Current in Full-Shell Nanowire Josephson Junctions". A preprint link and citation will be added upon submission.

## Abstract
Full-shell hybrid nanowires consisting of a semiconductor core fully enveloped by a superconducting shell have emerged as a platform to study Caroli-de Gennes-Matricon (CdGM) analogs. These subgap states can be considered a synthetic version of CdGM states in Abrikosov vortices. Unlike conventional CdGM states, these analogs exhibit a level spacing comparable to the superconducting gap, making them readily observable via tunneling spectroscopy techniques. The spectral density of CdGM analogs follows a characteristic skewed pattern as a function of applied axial magnetic field, an effect that is superimposed on the Little-Parks oscillations of the shell's gap induced by fluxoid quantization. Here, we provide experimental evidence for CdGM analogs through a distinctive skewness fingerprint in the critical current and zero-bias resistance of overdamped Josephson junctions based on full-shell nanowires. 

## About

This work compares full-shell Josephson junction theory against transport measurements, fitting the skewed critical current and core critical temperature observed in the experiments (data by Edu Lee et al.) to a microscopic model built with [Quantica.jl](https://github.com/pablosanjose/Quantica.jl) and [FullShell.jl](https://github.com/CarlosP24/FullShell.jl).

### Prerequisites

- Julia v1.12 or higher
- Recommended: 8+ CPU cores for parallel calculations

### Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/CarlosP24/SP45_Skweness_Experiment.git
   cd SP45_Skweness_Experiment
   ```

2. Install dependencies:

   ```bash
   julia --project=. -e "using Pkg; Pkg.instantiate()"
   ```

### Running Calculations

#### Local Execution

For local calculations with parallel processing (spawns 8 worker processes):

```bash
julia bin/launch_local.jl <name>
```

where `<name>` selects a wire (LDOS calculation) or junction (Josephson current calculation) registered in [src/models/wires.jl](src/models/wires.jl) and [src/models/junctions.jl](src/models/junctions.jl), e.g. `fitted_C3_left` or `fitted_C3`.

#### HPC Cluster Execution

Calculations can be submitted to a SLURM cluster via `make`, using the cluster definitions in [config/clusters.yaml](config/clusters.yaml):

```bash
make CLUSTER=esbirro ARG=fitted_C3 run
```

This generates a job launcher, deploys the code to the cluster, submits the job array, and syncs the results back to `data/` once finished. Available targets:

```bash
make                 # list available targets
make deploy           # sync local code to CLUSTER
make ARG=<arg> run     # run code on CLUSTER with ARG <arg>
make sync             # sync data from CLUSTER to local
```

### Generating Figures

1. Activate the plots environment and run the figure-generation script:

   ```bash
   cd plots
   julia --project=. plots.jl
   ```

   This includes and runs each `fig_*.jl` script in turn, saving the output.

All figures are saved as PDF files in `plots/figures/`.

## Repository Structure

```text
├── src/                     # Main source code
│   ├── main.jl              # Entry point, dispatches on wire/junction name
│   ├── main_local.jl        # Local driver used for interactive exploration
│   ├── builders/            # Tight-binding Hamiltonian construction
│   │   └── JosephsonJunction.jl
│   ├── operators/           # Green's function utilities
│   │   └── greens.jl
│   ├── calculations/        # Numerical calculation modules
│   │   ├── LDOS.jl          # Local density of states
│   │   └── Josephson.jl     # Josephson current
│   ├── models/               # Physical system definitions and fitted parameter sets
│   │   ├── params.jl
│   │   ├── wires.jl
│   │   └── junctions.jl
│   ├── parallelizers/        # Parallel computing utilities
│   └── utilities.jl
├── plots/                    # Figure generation (separate Julia environment)
│   ├── plots.jl               # Runs all fig_*.jl scripts
│   ├── fig_exp.jl              # Experimental data figure
│   ├── fig_theo_exp.jl          # Theory vs. experiment comparison figure
│   ├── plotters/                # Shared plotting panels
│   ├── kitchen/                  # Experimental data readers and critical-current extraction
│   └── figures/                   # Generated PDF figures
├── bin/                       # Execution scripts
│   ├── launch_local.jl         # Local parallel execution
│   └── launcher.jl              # SLURM job entry point
├── config/                    # Cluster and job configuration
├── data/                      # Computational results (JLD2 format)
└── data_exp/                  # Experimental measurements used in the figures
```

## Data

- `data/`: numerical results (LDOS, Josephson currents) produced by the calculations above, stored in JLD2 format.
- `data_exp/`: raw experimental transport measurements (lock-in and DC signals) used to produce the experimental and theory-vs-experiment figures.

## Dependencies

### Core Packages

- **[Quantica.jl](https://github.com/pablosanjose/Quantica.jl)** (v1.3.0): Quantum tight-binding calculations
- **[FullShell.jl](https://github.com/CarlosP24/FullShell.jl)** (v1.0.0): Full-shell hybrid nanowire utilities
- **[JLD2.jl](https://github.com/JuliaIO/JLD2.jl)**: Data serialization
- **[Parameters.jl](https://github.com/mauro3/Parameters.jl)**: Parameter management
- **[Roots.jl](https://github.com/JuliaMath/Roots.jl)** / **[QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl)** / **[Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)**: Root finding, quadrature and interpolation

### Computational Environment

- **[SlurmClusterManager.jl](https://github.com/JuliaParallel/SlurmClusterManager.jl)**: HPC integration
- **[ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl)**: Progress tracking

### Plotting (`plots/` environment)

- **[CairoMakie.jl](https://github.com/JuliaPlots/Makie.jl)**: Scientific plotting
- **[SavitzkyGolay.jl](https://github.com/lnacquaroli/SavitzkyGolay.jl)** / **[Peaks.jl](https://github.com/halleysfifthinc/Peaks.jl)**: Experimental data smoothing and peak extraction
- **[Bessels.jl](https://github.com/JuliaMath/Bessels.jl)**: Ambegaokar-Halperin critical current fitting
- **[CSV.jl](https://github.com/JuliaData/CSV.jl)** / **[DataFrames.jl](https://github.com/JuliaData/DataFrames.jl)** / **DelimitedFiles**: Data loading

## Citation

If you use this code or data in your research, please cite the associated preprint (citation details to be added upon submission).

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details.

## Support

For questions about the code or manuscript, please open a new issue with a detailed description or contact the corresponding authors listed in the manuscript.

## Acknowledgments
We acknowledge funding by the EU through the European Research Council (ERC) Starting Grant agreement 716559 (TOPOQDot), the FET-Open contract AndQC, the Danish National Research Foundation, Innovation Fund Denmark, the Carlsberg Foundation, and the Spanish AEI through Grants No.~PID2023-150224NB-I00, PRE2022-101362, PID2024-161156NB-I00, PID2024-161665NB-I00, and through the "María de Maeztu" Program for Units of Excellence in R\&D (CEX2018-000805-M), and Severo Ochoa Centers of Excellence program (CEX2024-001445-S); funded by MICIU/AEI/10.13039/501100011033, "ERDF A way of making Europe'' and "ESF+"".
