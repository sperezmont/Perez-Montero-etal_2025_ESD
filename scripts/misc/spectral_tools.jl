using Pkg
Pkg.activate(".")
using JLD2, DSP, Wavelets, ContinuousWavelets, Statistics, StatsBase, RecurrenceAnalysis

"""
    calc_spectrum(d, fs)
computes the power spectrum density of vector `d` using sample frequency `fs` and a Blackman window 
"""
function calc_spectrum(d::Vector, fs::Real)
    # -- spectrum blackman tuckey
    N = length(d)
    P = periodogram(d, fs=fs, window=blackman(N))
    G, freq = P.power, P.freq
    return G, freq
end

"""
    calc_wavelet(d, fs; sigma)
generate an array with the values of the wavelet applied to `d`
"""
function calc_wavelet(d::Vector, fs::Real; sigma::Real=π, s::Real=8)
    wvt = ContinuousWavelets.wavelet(Morlet(sigma), s=s, boundary=SymBoundary(), averagingType=NoAve(), β=1)
    S = ContinuousWavelets.cwt(d, wvt)
    freq = getMeanFreq(ContinuousWavelets.computeWavelets(length(d), wvt)[1], fs)
    S = abs.(S) .^ 2
    S = S ./ sum(S, dims=2)
    return S, freq
end

"""
    calc_coi(t, f, cf)
calculates the cone of influence from a wavelet analysis
adapted from Torrence and Compo, 1998, BAMS
https://es.mathworks.com/help/wavelet/ug/boundary-effects-and-the-cone-of-influence.html
"""
function calc_coi(t::Vector, f::Vector, cf::Real)
    predtimes = sqrt(2) .* cf ./ f
    tmax, tmin = maximum(t), minimum(t)
    t1 = tmin .+ predtimes
    t2 = tmax .- predtimes
    t_samples = vcat(t1, t2)
    p_samples = vcat(1 ./ f, 1 ./ f)
    return t_samples, p_samples
end
