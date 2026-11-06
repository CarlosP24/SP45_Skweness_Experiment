function include_all(dir::AbstractString)
    dir_path = endswith(pwd(), "src") ? joinpath(pwd(), dir) : joinpath(pwd(), "src", dir)
    for file in sort(readdir(dir_path))
        path = joinpath(dir_path, file)
        if endswith(file, ".jl") && isfile(path)
            include(path)
        end
    end
end

"""
    filter_πs(φrng)

Drop the phase points at 0, π and 2π, where the Josephson integrand is singular.

Lives here rather than beside a single calculation because both the multimode
(`calc_Josephson`) and the m_J-resolved (`calc_CPR`) paths use it, and so does
the diagnostic, which loads neither.
"""
function filter_πs(φrng)
    φrng1 = φrng[findall(x -> x > 0 && x < π, φrng)]
    φrng2 = φrng[findall(x -> x > π && x < 2π, φrng)]
    deleteat!(φrng1, findall(x -> isapprox(x, π), φrng1))
    deleteat!(φrng2, findall(x -> isapprox(x, 2π), φrng2))
    return vcat(φrng1, φrng2)
end

function ensure_pkgs(required_pkgs)
    for pkg in required_pkgs
        if !haskey(Pkg.project().dependencies, pkg)
            @info "Installing missing package: $pkg"
            Pkg.add(pkg)
        end
    end
end