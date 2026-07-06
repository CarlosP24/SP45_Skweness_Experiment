using Pkg
Pkg.activate("plots")
Pkg.instantiate()
using Parameters, JLD2, Quantica, FullShell, CSV, DataFrames, DelimitedFiles
using CairoMakie, Statistics, SavitzkyGolay, Roots, Bessels, Peaks
include("kitchen/data_reader.jl")
include("kitchen/Ic_reader.jl")
##