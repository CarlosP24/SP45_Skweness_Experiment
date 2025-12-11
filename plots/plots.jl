using Pkg
Pkg.activate("plots")
Pkg.instantiate()
using Parameters, JLD2, Quantica, FullShell, CSV, DataFrames, DelimitedFiles, Parameters
using CairoMakie, Statistics, SavitzkyGolay, Roots, Bessels, Peaks
include("kitchen/data_reader.jl")
include("kitchen/Ic_reader.jl")
include("plotters/LDOS.jl")
include("plotters/current.jl")
##