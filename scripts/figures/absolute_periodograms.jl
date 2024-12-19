using Pkg
Pkg.activate(".")
using JLD2, Interpolations, CairoMakie, LaTeXStrings, Colors

include(pwd() * "/scripts/misc/figure_definitions.jl")
include(pwd() * "/scripts/misc/spectral_tools.jl")
include(pwd() * "/scripts/misc/tools.jl")

Css = vcat(1e-10:1e-10:9e-10, 1e-9:1e-9:9e-9, 1e-8:1e-8:9e-8, 1e-7:1e-7:9e-7, 1e-6:1e-6:9e-6, 1e-5:1e-5:9e-5, 1e-4)
Cs_levels = [[1, 27, length(Css)], ["10^{-10}", "10^{-7}", "10^{-4}"]]
taukins = [1e3, 5e3, 1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5, 1e6]
taukin_levels = [[1, 3, 7, 11], ["1", "10^{1}", "10^{2}", "10^{3}"]]  # [10, 50, 100, 500, 1000, 2500, 5000, 7500, 10000]
taualphas = vcat(10e3:5e3:120e3)
taualpha_levels = [[1, 2, 12, 24, 25, 26, 27], ["1", "10", "60", "120", "400", "10^{3}", "10^{4}"]]   # tau_alpha = vcat(1e3, 10e3:5e3:120e3, 400e3, 1e6, 1e7)
cmap = [:purple, :royalblue, :olive, :darkorange, :darkred]
cmap13 = [:plum3, :purple, :slategray2, :lightsteelblue3, :cornflowerblue, :royalblue3, :royalblue4, :orange, :darkorange, :darkorange2, :darkorange4]
cmap15 = [:slategray2, :lightsteelblue3, :cornflowerblue, :royalblue3, :royalblue4, :darkseagreen1, :darkolivegreen2, :darkolivegreen3, :olive, :darkolivegreen]

function calc_absolute_periodogram(t::Vector, x::Vector, newt1::Real, newt2::Real, newdt::Real)
    t, x = interp_time_series(t, x, newt1, newt2, newdt)    # cut time series
    x_trend = sum(x) ./ length(x) .+ (x[end] .- x[1]) ./ (t[end] .- t[1]) .* t
    fs = 1 / (t[2] - t[1])
    G, fr = calc_spectrum(x .- x_trend, fs)
    G, fr = G[fr.>1/200e3], fr[fr.>1/200e3]   # filtering
    #G = G ./ sum(G)
    periods = 1 ./ fr
    return G, periods
end

function get_nc_from_exp(path2experiment::String)
    elements = readdir(path2experiment)
    exps = []
    for e in eachindex(elements)
        if elements[e][end-2:end] == ".nc"
            push!(exps, elements[e])
        end
    end
    return exps
end

function add_plot_to_axes!(ax::Axis, path2experiment::String; cmap=[:purple, :royalblue, :olive, :darkorange, :darkred], colors=[])
    # Get experiments
    experiments = get_nc_from_exp(path2experiment)

    # Compute colors
    if colors == []
        clrmap = cgrad(cmap, length(experiments), categorical=true)
        colors = collect(clrmap)
    else
        colors = colors
    end

    # Add things to ax
    ax.aspect = 1
    ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
    xlims!(ax, ylims_periods)

    # Compute and plot
    for r in eachindex(experiments)
        dfr = NCDataset("$(path2experiment)/$(experiments[r])")
        G, periods = calc_absolute_periodogram(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
        lines!(ax, periods ./ 1e3, G ./ 1e11, color=colors[r], linewidth=2)
    end

    return nothing

end


fig = Figure(resolution=(1500, 1000), fontsize=28)
ax_LIN, ax_NONLIN, ax_ISOS, ax_RISOS = Axis(fig[1, 1], title=L"LIN$\,$"), Axis(fig[1, 2], title=L"NONLIN$\,$"), Axis(fig[1, 3], title=L"ISOS$\,$"), Axis(fig[1, 4], title=L"RISOS$\,$")
ax_BASE, ax_THERM, ax_AGING = Axis(fig[2, 1], title=L"BASE$\,$"), Axis(fig[2, 2], title=L"THERM$\,$"), Axis(fig[2, 3], title=L"AGING$\,$")
ax_THERM_taukin, ax_AGING_taualpha = Axis(fig[3, 2], title=L"THERM-$\tau_\mathrm{kin}$$\,$"), Axis(fig[3, 3], title=L"AGING-$\tau_\alpha$$\,$")

add_plot_to_axes!(ax_LIN, pwd() * "/data/runs/exp01/LIN-MV/", cmap=cmap)
add_plot_to_axes!(ax_NONLIN, pwd() * "/data/runs/exp02/NONLIN-MV/", cmap=cmap)
add_plot_to_axes!(ax_ISOS, pwd() * "/data/runs/exp03/ISOS-MV/", cmap=cmap)
add_plot_to_axes!(ax_RISOS, pwd() * "/data/runs/exp03/RISOS-MV/", cmap=cmap)
add_plot_to_axes!(ax_BASE, pwd() * "/data/runs/exp04/BASE-MV-Cs/", cmap=cmap)
add_plot_to_axes!(ax_THERM, pwd() * "/data/runs/exp05/THERM-MV-Cs/", cmap=cmap)
add_plot_to_axes!(ax_AGING, pwd() * "/data/runs/exp06/AGING-MV-Cs/", cmap=cmap)

add_plot_to_axes!(ax_THERM_taukin, pwd() * "/data/runs/exp05/THERM-MV-taukin/", cmap=cmap13)

# This is a special case, we have runs with τα = 1e3, 10e3:5e3:120e3, 400e3, 1000e3, 10000e3
colormap_glacial = cgrad(cmap15, length(10e3:5e3:120e3), categorical=true)
colorsalpha = vcat(:purple, collect(colormap_glacial), :darkorange, :firebrick, :darkred) # add colors of individual curves

add_plot_to_axes!(ax_AGING_taualpha, pwd() * "/data/runs/exp06/AGING-MV-taualpha/", colors=colorsalpha)

ax_LIN.ylabel = L"$H_\mathrm{PSD}$ ($10^{11}\mathrm{m^2 ~ yr}$) $\,$"
ax_BASE.ylabel = L"$H_\mathrm{PSD}$ ($10^{11}\mathrm{m^2 ~ yr}$) $\,$"
ax_THERM_taukin.ylabel = L"$H_\mathrm{PSD}$ ($10^{11}\mathrm{m^2 ~ yr}$) $\,$"

ax_RISOS.xlabel = L"Period (kyr) $\,$"
ax_BASE.xlabel = L"Period (kyr) $\,$"
ax_THERM.xlabel = L"Period (kyr) $\,$"
ax_AGING.xlabel = L"Period (kyr) $\,$"
ax_THERM_taukin.xlabel = L"Period (kyr) $\,$"
ax_AGING_taualpha.xlabel = L"Period (kyr) $\,$"

for i in 1:4
    colsize!(fig.layout, i, Relative(1 / 4))
end

clrmap = cgrad(cmap, length(Css), categorical=true)
Colorbar(fig[2, 4], width=Relative(2 / 3), height=Relative(0.1), colormap=clrmap,
    label=L"$C_s$ ($\mathrm{m\cdot yr^{-1}\cdot Pa^{-2}}$)",
    limits=(1, length(Css) + 1),
    ticks=(Cs_levels[1] .+ 0.5, convert_strings_to_latex(Cs_levels[2])),
    ticklabelsize=28,
    vertical=false,
)

clrmap = cgrad(cmap13, length(taukins), categorical=true)
Colorbar(fig[3, 1], width=Relative(0.1), height=Relative(2 / 3), colormap=clrmap, flipaxis=true,
    label=L"$\tau_{kin}$ (kyr)",
    limits=(1, length(taukins) + 1),
    ticks=(taukin_levels[1] .+ 0.5, convert_strings_to_latex(taukin_levels[2])),
    ticklabelsize=28,
    vertical=true, halign=0.5
)

# Special case again, colorbar only includes 10 kyr to 120 kyr, I add four lines that represent 1kyr, 400kyr, 1000kyr and 10000kyr
clrmap = cgrad(cmap15, length(taualphas), categorical=true)
Colorbar(fig[3, 4], width=Relative(0.1), height=Relative(2 / 3), colormap=clrmap,
    label=L"$\tau_\alpha$ ($\mathrm{kyr}$)",
    limits=(1, length(taualphas) + 1),
    ticks=([1, 11, 23] .+ 0.5, convert_strings_to_latex(["10", "60", "120"])),
    ticklabelsize=28,
    vertical=true, halign=0.2
)
elem_1, elem_400, elem_1e3, elem_1e4 = [LineElement(color=:purple, linestyle=:solid),
    LineElement(color=:darkorange, linestyle=:solid),
    LineElement(color=:firebrick, linestyle=:solid),
    LineElement(color=:darkred, linestyle=:solid)]

Legend(fig[3, 4],
    [elem_1, elem_400, elem_1e3, elem_1e4],
    [L" $\tau_\alpha$ = 1 kyr", L" $\tau_\alpha$ = 400 kyr", L" $\tau_\alpha$ = $10^3$ kyr", L" $\tau_\alpha$ = $10^4$ kyr"],
    patchsize=(20, 30), rowgap=0, linewidth=3, framevisible=false, halign=1.5)


colgap!(fig.layout, 0.0)
colgap!(fig.layout, 1, -180.0)
colgap!(fig.layout, 2, -80.0)
colgap!(fig.layout, 3, -100.0)
rowgap!(fig.layout, 0.0)
rowgap!(fig.layout, 1, -50.0)
rowgap!(fig.layout, 2, 50.0)

resize_to_layout!(fig)

save("figures/absolute_periodograms.png", fig)