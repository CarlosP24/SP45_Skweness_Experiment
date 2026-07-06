function pfunction(func::Function, ranges; kw...)
    if length(ranges) != first(methods(func).ms).nargs - 1
        throw(ArgumentError("The number of arguments in func must match the number of ranges."))
    end
    pts = Iterators.product(ranges...)
    Fs = @showprogress pmap(pts) do args
        return func(args...; kw...)
    end
    return reshape(Fs, size(pts)...)
end