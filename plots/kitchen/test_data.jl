include("data_reader.jl")
data = read_data(Exp_Data(); exp_params = Exp_Params(; β = 0.95, Rs = 5490), window_size = 50, order = 3)

@unpack Is, dVdImat, Ts = data 
Bs = Ts

i0 = findmin(abs.(Is))[2]

R0 = dVdImat[:, i0] |> vec

fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"$B$ (T)", ylabel = L"$dV/dI (I=0)$ (Ω)",)

lines!(ax, Bs, R0, color = :navyblue, linewidth = 2)
ylims!(ax, (0, 1.2 * maximum(R0)))
fig

##
fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"$dV/dI$", ylabel = L"$$ Counts")
hist!(ax, vec(dVdImat), bins = 100, color = :navyblue)
xlims!(ax, 0, 2000)
fig

##
fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"$I$ (nA)", ylabel = L"$dV/dI$")
lines!(ax, Is, dVdImat[1, :], color = :navyblue, linewidth = 2)
fig


## Device D
exp_data = Exp_Data(;
    keysight = "data_exp/250mK_DEV13_277_18_dmm_keysight_signal.dat",
    lockin = "data_exp/250mK_DEV13_277_18_mf1_X.dat"
)
data = read_data(exp_data; exp_params = Exp_Params(; β = 0.95, Rs = 5490), window_size = 50, order = 3)
Imax = 120
Blims = (0.08, 0.4)
T = 460
Ts, Is, dVdImat, _, _ = prepare_data(exp_data, Imax)

Isw, _, _ = extract_peak_currents(dVdImat, Is, Ts)
RNs = compute_normal_resistances(dVdImat, Is, Isw)

ilim = findmin(abs.(Ts .- Blims[1]))[2]
imax = findmin(abs.(Ts .- Blims[2]))[2]

i0 = findmin(abs.(Is))[2]

rnorm = dVdImat ./ RNs
r0 = rnorm[:, i0] |> vec
γs = γAH.(abs.(r0))
IcAH = 0.0209 * T .* γs  # T in mK

Ic = vcat(NaN .* ones(1:ilim), IcAH[ilim+1:imax])

Trng = Ts[1:imax]

Rs = r0[ilim:imax]
Trng_Rs = Ts[ilim:imax]

fig = Figure()
ax = Axis(fig[1, 1], xlabel = L"$B$ (T)", ylabel = L"$R_0/R_N$",)
lines!(ax, Trng,  r0[1:imax])

ax = Axis(fig[2, 1], xlabel = L"$B$ (T)", ylabel = L"$I$",)
lines!(ax, Trng,  r0[1:imax])
fig

