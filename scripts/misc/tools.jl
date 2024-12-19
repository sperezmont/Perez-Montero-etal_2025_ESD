using Pkg
Pkg.activate(".")
using NCDatasets, JLD2, Interpolations, DSP, LaTeXStrings, CairoMakie

function calc_Ieff(df, par)
    albedo_eff = df["albedo"][:]
    albedo_eff[df["m"][:].>0] .= par.albedo_newice

    absortion = 1 .- albedo_eff
    sw = absortion .* (df["I"] .- params.insol_ref)
    Ieff = copy(sw)
    Ieff[df["I"][:].<=par.insol_ref] .= 0.0
    Ieff[df["H"][:].<=10.0] .= 0.0

    return Ieff
end

function find_points_in_cycle(time, ice, interval, threshold; offset=1)
    dt = time[2] .- time[1]
    ice = ice[interval]

    maxice = findmax(ice)[2] + interval[1] - offset

    res = abs.(ice .- threshold)
    dicedt = (ice[2:end] .- ice[1:end-1]) ./ dt

    starts, endings = [], []
    for i in range(2, length(dicedt))
        if (dicedt[i] > 0) && (res[i - 1] <= threshold)
            push!(starts, i - 1 + interval[1] - offset)
        elseif  (dicedt[i] < 0) && (res[i + 1] <= threshold)
            push!(endings, i + interval[1])
        end
    end

    if starts == []
        push!(starts, interval[1])
    end
    if endings == []
        push!(endings, interval[end])
    end
    
    return starts[1], maxice, endings[end]
end

function find_closest_cycle(tivec::Vector, cycles::Vector, period::Real)
    periods = [abs(tivec[idx[1]] - tivec[idx[2]]) for idx in cycles]
    closest_cycle = findmin(abs.(periods .- period))[2]
    return cycles[closest_cycle]
end

function find_glacial_cycle(data, indexes, period, threshold, minval)

    cycles_found = [(indexes[1], indexes[end])]
    for i in range(indexes[1], indexes[end] - period - threshold, step=1)
        value1 = abs(data[i])

        if value1 <= abs(minval)

            i2 = i + period - 1
            w2look = range(i2 - threshold, i2 + threshold, step=1)
            res = []
            for j in w2look
                push!(res, abs(data[j] - minval))
            end

            minres = findmin(res)[2]

            if abs(data[w2look[minres]]) <= abs(minval)
                idx1, idx2 = i, w2look[minres]
                push!(cycles_found, (idx1, idx2))
                #display(("Found! " * string(period * 10), idx1, idx2, data[idx1], data[idx2], minval))
            end

        end
    end

    return cycles_found
end

function convert_strings_to_latex(ticks::Vector)
    return [latexstring("$(i)") for i in ticks]
end

function butterworth_filt(x, fmax, fs; order=4)
    newx = vcat(reverse(x), vcat(x, reverse(x)))    # minimum derivative condition
    responsetype = Lowpass(fmax; fs=fs)
    designmethod = Butterworth(order)
    return filt(digitalfilter(responsetype, designmethod), newx)[length(x)+1:2*length(x)]
end

function interp_time_series(t, x, t0, tend, dt)
    # First, linear intepolation to dt time step
    nodes = (t,)
    x_interpolant = interpolate(nodes, x, Gridded(Linear()))#linear_interpolation(grid, x)
    t_interpolant = interpolate(nodes, t, Gridded(Linear()))#linear_interpolation(grid, t)

    new_nodes = t[1]:dt:t[end]
    x_interp = x_interpolant.(new_nodes)
    t_interp = t_interpolant.(new_nodes)

    # Second, cut time series to interval [t0, tend]
    start_index = findmin(abs.(t_interp .- t0))[2]
    end_index = findmin(abs.(t_interp .- tend))[2]

    return t_interp[start_index:end_index], x_interp[start_index:end_index]
end

function eliminate_missing_from_series(t, x)
    x_nomiss = collect(skipmissing(x))
    t_nomiss = t[broadcast(!, ismissing.(x))]
    new_x = Vector{Float64}(undef, length(x_nomiss))
    new_x[:] = x_nomiss
    return t_nomiss, new_x
end

function convert_nc_to_jld2(path2data::String, t0, tend, dt)
    df = NCDataset(path2data)
    filename = split(path2data, "/")[end][1:end-2] * "jld2"
    variables = keys(df)
    mat = [df[variables[1]][:] df[variables[2]][:]]' # creates a matrix whose rows are time, variable   

    t_nomiss, x_nomiss = eliminate_missing_from_series(mat[1, :], mat[2, :])    # eliminates missing for frequency analysis purposes
    mindt = findmin(abs.(t_nomiss[2:end] .- t_nomiss[1:end-1]))
    new_t_to_find_dt = t_nomiss
    while mindt[1] == 0.0
        new_t_to_find_dt = vcat(new_t_to_find_dt[1:mindt[2]], new_t_to_find_dt[mindt[2]+2:end])   # we eliminate the posibility of having same timesteps as in Lüthi et al., 2008 record
        mindt = findmin(abs.(new_t_to_find_dt[2:end] .- new_t_to_find_dt[1:end-1]))
    end

    if t_nomiss[1] < t_nomiss[end]
        tdtmin, xdtmin = interp_time_series(t_nomiss, x_nomiss, t_nomiss[1], t_nomiss[end], mindt[1]) # interpolates to minimum dt yr to normalize all the timesteps without losing frequencies
    else
        tdtmin, xdtmin = interp_time_series(reverse(t_nomiss), reverse(x_nomiss), t_nomiss[end], t_nomiss[1], mindt[1]) # interpolates to minimum dt yr to normalize all the timesteps without losing frequencies
    end
    
    xfiltered = butterworth_filt(xdtmin, 1 / 10e3, 1 / mindt[1])  # filters the data with a lowpass filter from 10 kyr of frequency

    # newt, newx = interp_time_series(tdtmin, xfiltered, t0, tend, dt) # interpolates to dt yr to normalize all the timesteps
    newt, newx = interp_time_series(tdtmin, xfiltered, t0, tend, mindt[1]) # cut time series
    mat = [newt newx]'

    JLD2.save_object("data/paleo_records/$(filename)", mat)
end

function create_PSD(t::Vector, x::Vector, newt1::Real, newt2::Real, newdt::Real)
    t, x = interp_time_series(t, x, newt1, newt2, newdt)    # cut time series
    x_trend = sum(x) ./ length(x) .+ (x[end] .- x[1]) ./ (t[end] .- t[1]) .* t
    fs = 1 / (t[2] - t[1])
    G, fr = calc_spectrum(x .- x_trend, fs)
    G, fr = G[fr.>1/200e3], fr[fr.>1/200e3]   # filtering
    G = G ./ sum(G)
    periods = 1 ./ fr
    return G, periods
end

function create_WV(t::Vector, x::Vector)
    x_trend = sum(x) ./ length(x) .+ (x[end] .- x[1]) ./ (t[end] .- t[1]) .* t
    sigma2use = 3.1
    scale2use = 10
    fs = 1 / (t[2] - t[1])
    Wnorm, freqs = calc_wavelet(x .- x_trend, fs; sigma=sigma2use, s=scale2use)
    periods = 1 ./ freqs
    return Wnorm, periods
end