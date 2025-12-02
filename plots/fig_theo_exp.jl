using DelimitedFiles
using Parameters

const πoΦ0 = 1.519267e-3  

@with_kw struct Exp_Params
    Rs::Float64 = 5000
    β::Float64 = 1
    dI::Float64 = 5e-9 * β
    G_0::Float64 = 7.74809173e-5
end

@with_kw struct Exp_Data
    keysight = "data_exp/250mK_DEV13_275_16_dmm_keysight_signal.dat"
    lockin = "data_exp/250mK_DEV13_275_16_mf1_X.dat"
end

@with_kw struct Filter_Params
    window_size::Int = 201
    poly_order::Int = 3
end

theo_params= Params_mm(;
    R = 60,
    ξd = 55,
    θ = 0.1,    
    Δ0 = 0.204
)

@with_kw struct Data
    Ts = Float64[]
    Is = Float64[]
    Vmat = Array{Float64,2}(undef, 0, 0)
    dVmat = Array{Float64,2}(undef, 0, 0)
    dIdVmat = Array{Float64,2}(undef, 0, 0)
    dVdImat = Array{Float64,2}(undef, 0, 0)
end

function read_data(exp_data::Exp_Data; exp_params::Exp_Params = Exp_Params())
    @unpack keysight, lockin = exp_data
    @unpack Rs, dI, G_0 = exp_params

    # load files (DelimitedFiles.readdlm returns an Array)


    keysight = DelimitedFiles.readdlm(keysight)
    lockin = DelimitedFiles.readdlm(lockin)

    # unique sorted indices for rows/columns
    Ts = sort(unique(keysight[:, 1]))
    Is = sort(unique(keysight[:, 2]))

    # mapping from value -> matrix index
    t_idx = Dict(val => i for (i, val) in enumerate(Ts))
    i_idx = Dict(val => j for (j, val) in enumerate(Is))

    # create matrix with missing values for pairs not present
    Vmat = fill(NaN, length(Ts), length(Is))
    dVmat = fill(NaN, length(Ts), length(Is))
    dIdVmat = fill(NaN, length(Ts), length(Is))
    dVdImat = fill(NaN, length(Ts), length(Is))
    # populate matrix: column 1 -> row index, column 2 -> col index, column 4 -> value
    for r in eachindex(keysight[:, 1])
        t = keysight[r, 1]
        i = keysight[r, 2]
        v = keysight[r, 4]
        dv = lockin[r, 4]
        # Real voltage and lock-in correction
        Vmat[t_idx[t], i_idx[i]] = v - Is[i_idx[i]] * Rs
        dVmat[t_idx[t], i_idx[i]] = dv - dI * Rs
        dIdVmat[t_idx[t], i_idx[i]] = dI / dVmat[t_idx[t], i_idx[i]] / G_0
        dVdImat[t_idx[t], i_idx[i]] = dVmat[t_idx[t], i_idx[i]] / dI
    end


    return Data(; Ts = Ts, Is = Is, Vmat = Vmat, dVmat = dVmat, dIdVmat = dIdVmat, dVdImat = dVdImat)
end

function find_lobes(Bs, theo_params)
    BtoΦ = get_Φ(theo_params)
    Φs = BtoΦ.(Bs)
    is = [findmin(abs.(Φs .- nh))[2] for nh in [0.5, 1.5, 2.5]]
    return is, Bs[is]
end

function find_Isw(dVdImat, filter_params)
    @unpack window_size, poly_order = filter_params
    sgfilter = SGolay(window_size, poly_order)
    iIsw = mapslices(dVdImat, dims=2) do dVdI_row
        filtered = dVdI_row |> sgfilter 
        imax = findmax(filtered.y)[2]
        dVdImean = mean(dVdI_row[imax:end])
        iIsw = findmin(abs.(filtered.y[1:imax] .- dVdImean))[2]
        return imax
    end
    return [iIsw...]
end

function integrate_Is(dVdImat)
    it = sum.(dVdImat)
    return [it...]
end

function plot_conductance(ax, name::String; dir = "data/Conductance")
    data = load(joinpath(dir, name * ".jld2"))["res"]
    @unpack Gs, params = data
    @unpack Brng = params
    lines!(ax, Brng, 1 ./ abs.(Gs) .|> log10; color = :blue)
    xlims!(ax, minimum(Brng), maximum(Brng))
    return ax
end

function Ic_AH(Rp, kBT, RN)
    α = Rp / RN
    if α < 1
        f(γ) = Rp/RN - 1 / besseli0(γ/2)^2
        γ = find_zero(f, (0, 100);)
        return 240 * kBT * γ # in nA
    end
    return 0
end

function fig_theo_exp(exp_data::Exp_Data, name::String; Imax = 130, B = 0, theo_params = theo_params, filter_params = Filter_Params())
    fig = Figure()
    data = read_data(exp_data)
    @unpack Ts, Is, dVdImat, dVmat = data 
    @unpack Rs = Exp_Params()
    Is = Is .* 1e9
    idx = findall(Is .<= Imax .&& Is .>= 0) 
    Is = Is[idx]
    dVdImat = dVdImat[:, idx]
    dVmat = dVmat[:, idx]
    
    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I$ (nA)",
        )
    heatmap!(ax, Ts, Is, log10.(abs.(dVdImat)); colormap = :thermal, colorrange = (2.9, 3.2))
    plot_Ic(ax, name; conv = 74.4, color = :white)

    vlines!(ax, [B]; color = :white, linestyle = :dash)
    is, Bis = find_lobes(Ts, theo_params)
    vlines!(ax, Bis; color = :red, linestyle = :dash)
    hidexdecorations!(ax)
    iIsw = find_Isw(dVdImat, filter_params)
    Isw = Is[iIsw]

    scatter!(ax, Ts, Isw; color = :lightblue,)
    Bx = findmin(abs.(Ts .- B))[2]
    ax = Axis(fig[1, 2];
        xlabel = L"$dV/dI$ ($\Omega$)",
        ylabel =  L"$I$ (nA)",
        )


    l = dVdImat[Bx, :]

    lines!(ax, l./maximum(l), Is; color = :blue)
    ylims!(ax, minimum(Is), maximum(Is))
    hideydecorations!(ax, ticks = false)

    # ax = Axis(fig[2, 1];
    #     xlabel = L"$B$ (T)",
    #     ylabel = L"$\int dV/dI$ (a.u.)",
    #     )
    
    # IdV = integrate_Is([dVdImat[i, 1:iIsw[i]] for i in eachindex(Ts)])
    # lines!(ax, Ts, IdV./ maximum(IdV); color = :green,)
    # vlines!(ax, Bis; color = :red, linestyle = :dash)
    # xlims!(ax, minimum(Ts), maximum(Ts))
    # hidexdecorations!(ax)

    ax = Axis(fig[2, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I_c$ (nA)",
        )
    lines!(ax, Ts, Ic_AH.(dVdImat[:, 1], 0.022, dVdImat[:, end]))
    vlines!(ax, Bis; color = :red, linestyle = :dash)
    #text!(ax, 0.35, 5.5; text = L"$I = 0$", color = :black, align = (:center, :center))
    xlims!(ax, minimum(Ts), maximum(Ts))

    ax = Axis(fig[2, 2];
        xlabel = L"$B$ (T)",
        ylabel = L"$I_c$ (nA)",
        )
    plot_Ic(ax, name; conv = 74.4,)
   # hideydecorations!(ax, ticks = false)
    vlines!(ax, Bis;  color = :red, linestyle = :dash)

    #rowgap!(fig.layout, 1, -30)

    # ax = Axis(fig[2, 2];
    #     xlabel = L"$B$ (T)",
    #     ylabel = L"$G$ (a. u.)",
    #     )
    # plot_conductance(ax, name; )
    return fig
end 

fig = fig_theo_exp(Exp_Data(), "fitted_C"; B = 0.32)
fig

##
j = findmin(abs.(Ts .- 0.05))[2]
dVdI_B = dVdImat[j, :]
filtered = savitzky_golay(dVdI_B, 101, 3;)

fig = Figure()
ax = Axis(fig[1, 1];
    xlabel = L"$I$ (nA)",
    ylabel = L"$dV/dI$ (a.u.)",
    )
lines!(ax, Is, dVdI_B; color = :blue, linestyle = :dash)
lines!(ax, Is, filtered.y; color = :red, linewidth = 4)
fig
