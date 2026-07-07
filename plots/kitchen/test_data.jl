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