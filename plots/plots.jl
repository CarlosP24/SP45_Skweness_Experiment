using Pkg
Pkg.activate("plots")
Pkg.instantiate()
using Parameters, JLD2, Quantica, FullShell, CSV, DataFrames, DelimitedFiles
using CairoMakie, Statistics, SavitzkyGolay, Roots, Bessels, Peaks
include("kitchen/data_reader.jl")
include("kitchen/Ic_reader.jl")
##

## Fig 2: Experiments
include("fig_exp.jl")
fig = fig_exp(Exp_Data(;
    keysight = "data_exp/250mK_DEV13_277_18_dmm_keysight_signal.dat",
    lockin = "data_exp/250mK_DEV13_277_18_mf1_X.dat"
    ),
    Exp_Data(;
    keysight = "data_exp/250mK_DEV13_294_35_V_drop.dat",
    lockin = "data_exp/250mK_DEV13_294_35_R_demod.dat"
    )
)
save("plots/figures/fig_exp.pdf", fig)
fig


## Fig 3: Experiments vs Theory
include("fig_theo_exp.jl")
fig = fig_theo_exp(Exp_Data(), "C3")
save("plots/figures/fig_theo_exp.pdf", fig)
fig
