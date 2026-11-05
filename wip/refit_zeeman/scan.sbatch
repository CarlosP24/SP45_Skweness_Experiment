#!/bin/bash -l
# Zeeman-corrected refit scan (Carlos approved cluster use, all clusters, 2026-11-05).
# Submit from the SP45 root:  sbatch [-p … -N … --ntasks-per-node …] wip/refit_zeeman/scan.sbatch <grid> <part> <tag>
# e.g.  sbatch -p long -N 4 --ntasks-per-node=96 wip/refit_zeeman/scan.sbatch scan1 1,2,3,4,5,6/19 atto_long
# Defaults below are atto's; command-line options override them.
#SBATCH --partition=long
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1800M
#SBATCH --time=12:00:00
#SBATCH --output=wip/refit_zeeman/logs/%j.out
#SBATCH --job-name=sp45_refit
#SBATCH --mail-user=carlos.paya@csic.es
#SBATCH --mail-type=END,FAIL

cd "$SLURM_SUBMIT_DIR"
# Julia 1.12.4 explicitly: juliaup directory overrides are not seen on every cluster's compute
# nodes, and on CESGA a system Julia 1.11 precedes juliaup in PATH there, so pass the absolute
# binary with --export=ALL,JULIA_BIN=~/.julia/juliaup/julia-1.12.4+0.x64.linux.gnu/bin/julia
JULIA=${JULIA_BIN:-julia +1.12.4}
$JULIA --project=. wip/refit_zeeman/scan.jl "${1:-scan1}" "${2:-1/1}" "${3:-all}"
