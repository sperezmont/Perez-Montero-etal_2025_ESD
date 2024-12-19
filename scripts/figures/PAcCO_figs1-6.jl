using Pkg
Pkg.activate(".")
using JLD2, Interpolations, CairoMakie, LaTeXStrings, Colors

include(pwd() * "/scripts/misc/figure_definitions.jl")
include(pwd() * "/scripts/misc/spectral_tools.jl")
include(pwd() * "/scripts/misc/tools.jl")

Css = vcat(1e-10:1e-10:9e-10, 1e-9:1e-9:9e-9, 1e-8:1e-8:9e-8, 1e-7:1e-7:9e-7, 1e-6:1e-6:9e-6, 1e-5:1e-5:9e-5, 1e-4)
Cs_levels = [[1, 27, length(Css)], ["10^{-10}", "10^{-7}", "10^{-4}"]]
cmap = [:purple, :royalblue, :olive, :darkorange, :darkred]

function calc_zthr(Tsl::Vector, Tref::Vector, p)
    return (p.sref .+ p.ks .* (Tsl .- Tref) .- p.lambda .* (Tsl .- p.Tthreshold)) / ((p.ks .- p.lambda) .* p.Γ)
end


# ===========================================================================================================
# Figure LIN-MV
# ===========================================================================================================
figname = "exp01_LIN-MV_PSD"
df1 = NCDataset("data/runs/exp01/LIN-MV/LIN-MV01.nc")
params = JLD2.load_object("data/runs/exp01/LIN-MV/LIN-MV01_params.jld2")
runs2plot = range(1, length(Css))
colormap = cgrad(cmap, length(runs2plot), categorical=true)
colors = collect(colormap)

xticks = Int.(-8e2:1e2:0)
xlims = (-8e2, 0)
ylims = (-600, 5000)
yticks = Int.(0:1500:3000)

fig = Figure(resolution=(1500, 750), fonts=(; regular="TeX"), fontsize=28)

# Panel a
ax = Axis(fig[1, 1], ylabel=ins_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks_insol, convert_strings_to_latex(yticks_insol))
ylims!(ax, (400, 600))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b)
text!(ax, -7.95e2, 580, text=a_label, align=(:left, :center))
hlines!(ax, [params.insol_ref], color=background_color, linestyle=:solid)
lines!(ax, df1["time"] ./ 1e3, df1["I"], color=pacco_color, linewidth=3)

# Panel b
ax = Axis(fig[1, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.4))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.3, text=b_label, align=(:left, :center))
G, periods = create_PSD(df1["time"][:], df1["I"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=pacco_color, linewidth=2.5)

# Panel c
ax = Axis(fig[2, 1], ylabel=icethick_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks, convert_strings_to_latex(yticks))
ylims!(ax, ylims)
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b, :t)
text!(ax, -7.95e2, 4000, text=L"(c)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs/exp01/LIN-MV/LIN-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[r])
end

# Panel d
ax = Axis(fig[2, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.4))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.3, text=d_label, align=(:left, :center))
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp01/LIN-MV/LIN-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

# Panel e
ax = Axis(fig[3, 1], ylabel=v_label, xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 200, 400, 600], convert_strings_to_latex([0, 200, 400, 600]))
ylims!(ax, (-10, 700))
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -7.95e2, 650, text=L"(e)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp01/LIN-MV/LIN-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["v"], color=colors[r], linewidth=2)
end

# Panel f
ax = Axis(fig[3, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.4))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.3, text=f_label, align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp01/LIN-MV/LIN-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["v"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

Colorbar(fig[2:3, 3], width=Relative(0.3), height=Relative(1 / 3), colormap=colormap,
    label=L"$C_s$ ($\mathrm{m\cdot yr^{-1}\cdot Pa^{-2}}$)",
    limits=(1, length(runs2plot) + 1),
    ticks=(Cs_levels[1] .+ 0.5, convert_strings_to_latex(Cs_levels[2])),
    ticklabelsize=28,
    vertical=true, halign=0.0,
)

rowgap!(fig.layout, 0.0)
colsize!(fig.layout, 2, Relative(2 / 6))
colsize!(fig.layout, 3, Relative(0.04))
resize_to_layout!(fig)

save("figures/fig03.pdf", fig)

# ===========================================================================================================
# Figure NONLIN-MV
# ===========================================================================================================
figname = "exp02_NONLIN-MV_PSD"
df1 = NCDataset("data/runs//exp02/NONLIN-MV/NONLIN-MV01.nc")
params = JLD2.load_object("data/runs//exp02/NONLIN-MV/NONLIN-MV01_params.jld2")
runs2plot = range(1, length(Css))
colormap = cgrad(cmap, length(runs2plot), categorical=true)
colors = collect(colormap)

xticks = Int.(-8e2:1e2:0)
xlims = (-8e2, 0)
ylims = (-600, 5000)
yticks = Int.(0:1500:9000)

fig = Figure(resolution=(1500, 750), fonts=(; regular="TeX"), fontsize=28)

# Panel a
ax = Axis(fig[1, 1], ylabel=ins_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks_insol, convert_strings_to_latex(yticks_insol))
ylims!(ax, (400, 600))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b)
text!(ax, -7.95e2, 580, text=a_label, align=(:left, :center))
hlines!(ax, [params.insol_ref], color=background_color, linestyle=:solid)
lines!(ax, df1["time"] ./ 1e3, df1["I"], color=pacco_color, linewidth=3)

# Panel b
ax = Axis(fig[1, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.4))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.3, text=b_label, align=(:left, :center))
G, periods = create_PSD(df1["time"][:], df1["I"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=pacco_color, linewidth=2.5)

# Panel c
ax = Axis(fig[2, 1], ylabel=icethick_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks, convert_strings_to_latex(yticks))
ylims!(ax, (-600, 5500))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b, :t)
text!(ax, -7.95e2, 5200, text=L"(c)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp02/NONLIN-MV/NONLIN-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[r])
end

# Panel d
ax = Axis(fig[2, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.4))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.3, text=d_label, align=(:left, :center))
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp02/NONLIN-MV/NONLIN-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

# Panel e
ax = Axis(fig[3, 1], ylabel=v_label, xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 200, 400, 600], convert_strings_to_latex([0, 200, 400, 600]))
ylims!(ax, (-10, 700))
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -7.95e2, 650, text=L"(e)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp02/NONLIN-MV/NONLIN-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["v"], color=colors[r], linewidth=2)
end

# Panel f
ax = Axis(fig[3, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.4))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.3, text=f_label, align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp02/NONLIN-MV/NONLIN-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["v"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

Colorbar(fig[2:3, 3], width=Relative(0.3), height=Relative(1 / 3), colormap=colormap,
    label=L"$C_s$ ($\mathrm{m\cdot yr^{-1}\cdot Pa^{-2}}$)",
    limits=(1, length(runs2plot) + 1),
    ticks=(Cs_levels[1] .+ 0.5, convert_strings_to_latex(Cs_levels[2])),
    ticklabelsize=28,
    vertical=true, halign=0.0,
)

rowgap!(fig.layout, 0.0)
colsize!(fig.layout, 2, Relative(2 / 6))
colsize!(fig.layout, 3, Relative(0.04))
resize_to_layout!(fig)

save("figures/fig05.pdf", fig)

# ===========================================================================================================
# Figure ISOS-MV
# ===========================================================================================================
figname = "exp03_ISOS-MV_PSD"
df1 = NCDataset("data/runs//exp03/ISOS-MV/ISOS-MV01.nc")
params = JLD2.load_object("data/runs//exp03/ISOS-MV/ISOS-MV01_params.jld2")
runs2plot = range(1, length(Css))
colormap = cgrad(cmap, length(runs2plot), categorical=true)
colors = collect(colormap)

xticks = Int.(-8e2:1e2:0)
xlims = (-8e2, 0)
ylims = (-600, 5000)
yticks = Int.(0:1500:9000)

fig = Figure(resolution=(1500, 750), fonts=(; regular="TeX"), fontsize=28)

# Panel a
ax = Axis(fig[1, 1], ylabel=ins_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks_insol, convert_strings_to_latex(yticks_insol))
ylims!(ax, (400, 600))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b)
text!(ax, -7.95e2, 580, text=a_label, align=(:left, :center))
hlines!(ax, [params.insol_ref], color=background_color, linestyle=:solid)
lines!(ax, df1["time"] ./ 1e3, df1["I"], color=pacco_color, linewidth=3)

# Panel b
ax = Axis(fig[1, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.4))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.3, text=b_label, align=(:left, :center))
G, periods = create_PSD(df1["time"][:], df1["I"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=pacco_color, linewidth=2.5)

# Panel c
ax = Axis(fig[2, 1], ylabel=icethick_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks, convert_strings_to_latex(yticks))
ylims!(ax, (-600, 6500))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b, :t)
text!(ax, -7.95e2, 6200, text=L"(c)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp03/ISOS-MV/ISOS-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[r], linewidth=2)
end

# Panel d
ax = Axis(fig[2, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.4))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.3, text=d_label, align=(:left, :center))
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp03/ISOS-MV/ISOS-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

# Panel e
ax = Axis(fig[3, 1], ylabel=v_label, xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 100, 200, 300, 400], convert_strings_to_latex([0, 100, 200, 300, 400]))
ylims!(ax, (-10, 700))
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -7.95e2, 650, text=L"(e)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp03/ISOS-MV/ISOS-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["v"], color=colors[r], linewidth=2)
end

# Panel f
ax = Axis(fig[3, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.4))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.3, text=f_label, align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp03/ISOS-MV/ISOS-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["v"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

Colorbar(fig[2:3, 3], width=Relative(0.3), height=Relative(1 / 3), colormap=colormap,
    label=L"$C_s$ ($\mathrm{m\cdot yr^{-1}\cdot Pa^{-2}}$)",
    limits=(1, length(runs2plot) + 1),
    ticks=(Cs_levels[1] .+ 0.5, convert_strings_to_latex(Cs_levels[2])),
    ticklabelsize=28,
    vertical=true, halign=0.0,
)

rowgap!(fig.layout, 0.0)
colsize!(fig.layout, 2, Relative(2 / 6))
colsize!(fig.layout, 3, Relative(0.04))
resize_to_layout!(fig)

save("figures/fig07.pdf", fig)

# ===========================================================================================================
# Figure 5. NONLIN-MV vs ISOS-MV (analysis)
# ===========================================================================================================
figname = "exp03_NONLINvsISOS_zthr"
selected_index = 30
runs2plot = ["exp02/NONLIN-MV/NONLIN-MV$(string(selected_index))", "exp03/ISOS-MV/ISOS-MV$(string(selected_index))"]
df1 = NCDataset("data/runs//exp03/ISOS-MV/ISOS-MV$(string(selected_index)).nc")
params = JLD2.load_object("data/runs//exp03/ISOS-MV/ISOS-MV$(string(selected_index))_params.jld2")
#zthr = (df1["Tsl"] .- params.Tthreshold) ./ params.Γ    # only dependent on the forcing (the same for both runs)
zthr = calc_zthr(df1["Tsl"][:], df1["Tref"][:], params)
print(params.Cs)
colors = [:steelblue, :darkorange]
labels = [L"Isostasy OFF$\,$", L"Isostasy ON$\,$"]

xticks = Int.(-8e2:1e2:0)
xlims = (-8e2, 0)
ylims = (-1500, 5000)
yticks = Int.(0:1500:3000)

fig = Figure(resolution=(1300, 500), fonts=(; regular="TeX"), fontsize=28)

# Panel a
ax = Axis(fig[1, 1], ylabel=z_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks, convert_strings_to_latex(yticks))
ylims!(ax, ylims)
xlims!(ax, xlims)
hidespines!(ax, :b)
hidexdecorations!(ax)
text!(ax, -7.95e2, 4000, text=L"(a)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
lines!(ax, df1["time"] ./ 1e3, zthr, color=:darkred, linewidth=2, strokecolor=:darkred, strokewidth=0.0001, label=L"$z_\mathrm{thr}$")
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs/$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["z"], color=colors[r], linewidth=4, label=labels[r])
end
axislegend(ax, framevisible=false, position=:rt, labelsize=25, nbanks=3)

# Panel b
ax = Axis(fig[2, 1], xlabel=time_label, ylabel=L"$T_{surf} - T_\mathrm{ref}$ (K)", xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks_temp, convert_strings_to_latex(yticks_temp))
ylims!(ax, (-40, 30))
xlims!(ax, xlims)
# hidexdecorations!(ax)
hidespines!(ax, :t)
text!(ax, -7.95e2, 20, text=L"(b)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs/$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["Tsurf"] .- dfr["Tref"], color=colors[r], linewidth=4, label=labels[r])
end
hlines!(ax, params.Tthreshold - params.degK, color=:darkred, linestyle=:solid)

# Panel c
# ax = Axis(fig[3, 1], xlabel=time_label, ylabel=L"m$\cdot$yr$^{-1}$", xgridvisible=false, ygridvisible=false)
# ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
# #ax.yticks = (yticks_temp, convert_strings_to_latex(yticks_temp))
# ylims!(ax, (-5, 1))
# xlims!(ax, xlims)
# hidespines!(ax, :t)
# text!(ax, -7.95e2, 20, text=L"c)$\,$", align=(:left, :center))
# for r in eachindex(runs2plot)
#     dfr = NCDataset("data/runs//$(string(r, pad=ndigits(length(runs2plot)))).nc")
#     lines!(ax, dfr["time"] ./ 1e3, dfr["s"] .- dfr["a"], color=colors[r], linewidth=2)
#     lines!(ax, dfr["time"] ./ 1e3, -1 .* dfr["v"] .* dfr["H"] ./ params.L, color=colors[r], linestyle=:dash, linewidth=2)
# end
# dfr = NCDataset("data/runs//$(runs2plot[1]).nc")
# lines!(ax, dfr["time"] ./ 1e3, (dfr["s"] .- dfr["a"]) .* NaN, color=:black, linewidth=2, label=L"$\dot{s} \mathrm{-} \dot{a}$")
# lines!(ax, dfr["time"] ./ 1e3, -1 .* dfr["v"] .* dfr["H"] ./ params.L  .* NaN, color=:black, linestyle=:dash, linewidth=2, label=L"$\mathrm{-}v \cdot H / L$")
# axislegend(ax, framevisible=false, position=:rb, labelsize=25, nbanks=3)

# # Panel c
# ax = Axis(fig[1:2, 2], ylabel=z_label, xlabel=L"$T_{surf} - T_\mathrm{ref}$ (K)", xgridvisible=false, ygridvisible=false, yaxisposition=:right, aspect=1)
# ax.yticks = ([0, 1500, 3000], convert_strings_to_latex([0, 1500, 3000]))
# ax.xticks = ([-20, -10, 0], convert_strings_to_latex([-20, -10, 0]))
# text!(ax, -35, 3000, text=L"c)$\,$", align=(:left, :center))
# lines!(ax, df1["Tsurf"] .- df1["Tref"], zthr, zeros(length(zthr)), color=:darkred)
# for r in eachindex(runs2plot)
#     dfr = NCDataset("data/runs//$(string(r, pad=ndigits(length(runs2plot)))).nc")
#     lines!(ax, dfr["Tsurf"] .- dfr["Tref"], dfr["z"], color=colors[r], linewidth=2, label=labels[r])
# end
# vlines!(ax, params.Tthreshold - params.degK, color=:darkred, linestyle=:solid)

#colsize!(fig.layout, 2, Relative(2 / 6))
rowgap!(fig.layout, 0.0)

save("figures/fig08.pdf", fig)

# ===========================================================================================================
# Figure 6. RISOS-MV
# ===========================================================================================================
figname = "exp03_RISOS-MV_PSD"
df1 = NCDataset("data/runs/exp03/RISOS-MV/RISOS-MV01.nc")
runs2plot = range(1, length(Css))
colormap = cgrad(cmap, length(runs2plot), categorical=true)
colors = collect(colormap)

xticks = Int.(-8e2:1e2:0)
xlims = (-8e2, 0)
ylims = (-600, 5000)
yticks = Int.(0:1500:7000)

fig = Figure(resolution=(1500, 500), fonts=(; regular="TeX"), fontsize=28)

# Panel a
ax = Axis(fig[1, 1], ylabel=ins_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks_insol, convert_strings_to_latex(yticks_insol))
ylims!(ax, (400, 600))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b)
text!(ax, -7.95e2, 580, text=a_label, align=(:left, :center))
hlines!(ax, [params.insol_ref], color=background_color, linestyle=:solid)
lines!(ax, df1["time"] ./ 1e3, df1["I"], color=pacco_color, linewidth=3)

# Panel b
ax = Axis(fig[1, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=b_label, align=(:left, :center))
G, periods = create_PSD(df1["time"][:], df1["I"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=pacco_color, linewidth=2.5)

# Panel c
ax = Axis(fig[2, 1], xlabel=time_label, ylabel=icethick_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks, convert_strings_to_latex(yticks))
ylims!(ax, (-600, 8500))
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -7.95e2, 6500, text=L"(c)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp03/RISOS-MV/RISOS-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[r], linewidth=2)
end

# Panel d
ax = Axis(fig[2, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.23, text=d_label, align=(:left, :center))
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp03/RISOS-MV/RISOS-MV$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

Colorbar(fig[2, 3], width=Relative(0.3), height=Relative(2 / 3), colormap=colormap,
    label=L"$C_s$ ($\mathrm{m\cdot yr^{-1}\cdot Pa^{-2}}$)",
    limits=(1, length(runs2plot) + 1),
    ticks=(Cs_levels[1] .+ 0.5, convert_strings_to_latex(Cs_levels[2])),
    ticklabelsize=28,
    vertical=true, halign=0.0,
)

rowgap!(fig.layout, 0.0)
colsize!(fig.layout, 2, Relative(2 / 6))
colsize!(fig.layout, 3, Relative(0.04))
resize_to_layout!(fig)

save("figures/fig09.pdf", fig)
