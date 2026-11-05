#!/bin/bash -l
# Devices-B/C fit scan (cluster use approved by Carlos, 2026-11-05).
# Submit from the SP45 root:  sbatch [-p … -N … --ntasks-per-node …] refit_BC/scan.sbatch <grid> <part> <tag>
# e.g.  sbatch -p esbirro -N 8 --ntasks-per-node=32 refit_BC/scan.sbatch B 1/1 esbirro
# Defaults below are atto's; command-line options override them.
#SBATCH --partition=long
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1800M
#SBATCH --time=12:00:00
#SBATCH --output=refit_BC/logs/%j.out
#SBATCH --job-name=sp45_BC
#SBATCH --mail-user=carlos.paya@csic.es
#SBATCH --mail-type=END,FAIL

cd "$SLURM_SUBMIT_DIR"
mkdir -p refit_BC/logs refit_BC/out
# Julia 1.12.4 explicitly (see wip/refit_zeeman/scan.sbatch): JULIA_BIN overrides
JULIA=${JULIA_BIN:-julia +1.12.4}
$JULIA --project=. refit_BC/scan.jl "${1:-smokeB}" "${2:-1/1}" "${3:-all}"
