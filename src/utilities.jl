function include_all(dir::AbstractString)
    dir_path = endswith(pwd(), "src") ? joinpath(pwd(), dir) : joinpath(pwd(), "src", dir)
    for file in sort(readdir(dir_path))
        path = joinpath(dir_path, file)
        if endswith(file, ".jl") && isfile(path)
            include(path)
        end
    end
end

function ensure_pkgs(required_pkgs)
    for pkg in required_pkgs
        if !haskey(Pkg.project().dependencies, pkg)
            @info "Installing missing package: $pkg"
            Pkg.add(pkg)
        end
    end
end