function γAH(R)
    """Compute the Ambegaokar-Helkin normalized current γ from normal resistance ratio R."""
    if R < 1
        f(γ) = R - besseli0(γ/2)^(-2)
        γ = find_zero(f, (0, 100);)
        return γ
    end
    return NaN
end


function select_top_peaks(all_indices::Vector, Is::Vector, n::Int=2)
    """Select the n highest peaks by Is value."""
    if length(all_indices) > n
        sorted_idx = sortperm(Is[all_indices], rev=true)
        return all_indices[sorted_idx[1:n]]
    else
        return all_indices
    end
end

function get_peak_indices_smoothed(dVdI_row::Vector; window=101, order=1, min_prom=6)
    """Find peak indices in a smoothed dVdI curve."""
    rf = savitzky_golay(dVdI_row, window, order).y
    pks = findmaxima(rf;)
    pks = peakproms(pks; min = min_prom)
    return pks.indices
end

function extract_peak_currents(dVdImat::Matrix, Is::Vector, Ts::Vector)
    """Extract critical currents from peak analysis across all B field values."""
    Ics = map(axes(dVdImat, 1)) do ir
        all_indices = get_peak_indices_smoothed(dVdImat[ir, :])
        selected = select_top_peaks(all_indices, Is, 2)
        
        if length(selected) == 2
            Ic_u, Ic_l = Is[selected[1]], Is[selected[2]]
            if Ic_u < Ic_l
                Ic_u, Ic_l = Ic_l, Ic_u
            end
            return Ic_u, Ic_l
        elseif length(selected) == 1
            if Is[selected[1]] > 0
                return Is[selected[1]], NaN
            else
                return NaN, Is[selected[1]]
            end
        else
            return NaN, NaN
        end
    end
    
    Icu = getindex.(Ics, 1)
    Icl = getindex.(Ics, 2)
    
    Ic = map(zip(Icu, Icl)) do (u, l)
        if abs(u + l) > 20
            return NaN
        end
        return (u - l) / 2
    end
    
    return Ic, Icu, Icl
end

function compute_normal_resistances(dVdImat::Matrix, Is::Vector, Ic::Vector)
    """Compute normal-state resistances."""
    RNs = map(axes(dVdImat, 1)) do ir
        if !isnan(Ic[ir]) && Ic[ir] != 0
            threshold = abs(Ic[ir]) + 20
        else
            threshold = 20
        end
        iIs = findall(x -> abs(x) >= threshold, Is)
        if !isempty(iIs)
            return mean(dVdImat[ir, iIs])
        else
            return NaN
        end
    end
    return RNs
end

function find_minima_indices(Rf::Vector; min_prom=2)
    """Find minima indices in Rf at I=0."""
    pks = findminima(Rf;)
    pks = peakproms(pks; min = min_prom)
    indeces = pks.indices[2:end] 
    pushfirst!(indeces, 1)
    return indeces
end