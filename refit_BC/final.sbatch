#!/bin/bash -l
# Final curves (fine-grid I_c and LDOS) of a selected candidate, refit_BC/final.jl, on a cluster.
# Needs refit_BC/out/<grid>_<mode>_scores.jld2 (analyze.jl) next to it on the cluster.
# Submit from the SP45 root:  sbatch [-p … -N … --ntasks-per-node …] refit_BC/final.sbatch <grid> [mode] [rank]
#SBATCH --partition=esbirro
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1900M
#SBATCH --time=03:00:00
#SBATCH --output=refit_BC/logs/%j.out
#SBATCH --job-name=sp45_BC_final

cd "$SLURM_SUBMIT_DIR"
mkdir -p refit_BC/logs refit_BC/out
JULIA=${JULIA_BIN:-julia +1.12.4}
$JULIA --project=. refit_BC/final.jl "${1:-B}" "${2:-fixed}" "${3:-1}"
