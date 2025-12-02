using Pkg
Pkg.activate("plots")
Pkg.instantiate()
using Parameters, JLD2, Quantica, FullShell, CSV, DataFrames
using CairoMakie, DelimitedFiles, Statistics, SavitzkyGolay, Roots, Bessels
include("plotters/LDOS.jl")
include("plotters/current.jl")
##