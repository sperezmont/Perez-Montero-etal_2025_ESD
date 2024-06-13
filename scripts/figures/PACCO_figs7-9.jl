using Pkg
Pkg.activate(".")
using JLD2, NCDatasets, CairoMakie, Statistics, LaTeXStrings

function cut_run(t, x, time2cut)
    minidx = findmin(abs.(t .- time2cut))[2]
    return x[minidx:end]
end

include(pwd() * "/scripts/misc/figure_definitions.jl")
include(pwd() * "/scripts/misc/spectral_tools.jl")
include(pwd() * "/scripts/misc/tools.jl")

lisiecki2005_temp = JLD2.load_object("data/paleo_records/lisiecki-raymo_2005_T_800.jld2")
hodell2023_temp = JLD2.load_object("data/paleo_records/hodell-etal_2023_T_800.jld2")
barker2011_temp = JLD2.load_object("data/paleo_records/barker-etal_2011.jld2")
bintanja2008_T = JLD2.load_object("data/paleo_records/bintanja-vandewal_2008_T_800.jld2")
luthi2008_co2 = JLD2.load_object("data/paleo_records/luthi-etal_2008.jld2")
yamamoto2022_co2 = JLD2.load_object("data/paleo_records/yamamoto-etal_C_800.jld2")
berends2021_co2 = JLD2.load_object("data/paleo_records/berends-etal_2020_C_800.jld2")
spratt2016_vol = JLD2.load_object("data/paleo_records/spratt-lisiecki_2016.jld2")
bintanja2008_vol = JLD2.load_object("data/paleo_records/bintanja-vandewal_2008_800.jld2")
berends2021_vol = JLD2.load_object("data/paleo_records/berends-etal_2020_Vol_800.jld2")

Cs_levels = [[1, 13, 23], ["10^{-10}", "10^{-7}", "10^{-4}"]]
taukin_levels = [[1, 3, 5, 9], ["10^{1}", "10^{2}", "10^{3}", "10^{4}"]]  # [10, 50, 100, 500, 1000, 2500, 5000, 7500, 10000]
taualpha_levels = [[1, 3, 7, 11], ["10^{3}", "10^{4}", "10^{5}", "10^{6}"]]   # [1e3, 5e3, 1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5, 1e6]
hgeo_levels = [[1, 9, 19], ["10^{-3}", L"5\cdot10^{-3}", "10^{-2}"]] # collect(1:0.5:10) .* 1e-3

cmap = [:black, :purple, :royalblue, :olive, :darkorange, :darkred]
cmap2 = [:lightblue, :royalblue, :black, :goldenrod, :darkorange, :darkred]

# Ticks
xticks_time = (xticks_time_800, convert_strings_to_latex(-1 .* xticks_time_800))

function calc_zthr(Tsl::Vector, Tref::Vector, p)
    return (p.sref .+ p.ks .* (Tsl .- Tref) .- p.lambda .* (Tsl .- p.Tthreshold)) / ((p.ks .- p.lambda) .* p.Γ)
end

# ===========================================================================================================
# Figure BASE-Cs
# ===========================================================================================================

figname = "exp04_BASE-Cs_PSD"
runs2plot = range(1, 23)
df1 = NCDataset("data/runs//exp04/BASE-Cs/BASE-Cs1.nc")
params = JLD2.load_object("data/runs//exp04/BASE-Cs/BASE-Cs1_params.jld2")

colormap = cgrad(cmap, length(runs2plot), categorical=true)
colors = collect(colormap)

xticks = Int.(-8e2:1e2:0)
xlims = (-8e2, 0)
ylims = (-600, 3000)
yticks = Int.(0:1000:3000)

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
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=b_label, align=(:left, :center))
G, periods = create_PSD(df1["time"][:], df1["I"][:])
lines!(ax, periods ./ 1e3, G, color=pacco_color, linewidth=2.5)

# Panel c
ax = Axis(fig[2, 1], ylabel=icethick_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks, convert_strings_to_latex(yticks))
ylims!(ax, ylims)
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b, :t)
text!(ax, -7.95e2, 2500, text=L"(c)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp04/BASE-Cs/BASE-Cs$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[r])
end

# Panel d
ax = Axis(fig[2, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=d_label, align=(:left, :center))
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp04/BASE-Cs/BASE-Cs$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

# Panel e
ax = Axis(fig[3, 1], ylabel=temp_label, xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks_temp, convert_strings_to_latex(yticks_temp))
ylims!(ax, ylims_temp)
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -7.95e2, 10, text=L"(e)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp04/BASE-Cs/BASE-Cs$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[r], linewidth=2)
end

strt_idx = findmin(abs.(bintanja2008_T[1, :] .+ 125e3))[2]
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :], color=:black, linewidth=3.5, linestyle=:dash, label=L"Bintanja and van de Wal (2008) $\,$")
strt_idx = findmin(abs.(barker2011_temp[1, :] .+ 125e3))[2]
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :], color=:grey25, linewidth=3.5, linestyle=:dash, label=L"Barker et al. (2011) $\,$")
axislegend(ax, framevisible=false, position=:rt, labelsize=25, nbanks=2, patchsize=(40, 20))

# Panel f
ax = Axis(fig[3, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.23, text=f_label, align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp04/BASE-Cs/BASE-Cs$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["T"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :])
lines!(ax, periods ./ 1e3, G, color=:black, linestyle=:dash, linewidth=3.5)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :])
lines!(ax, periods ./ 1e3, G, color=:grey25, linestyle=:dash, linewidth=3.5)

Colorbar(fig[1, 2], width=Relative(1.2 / 3), height=Relative(1 / 10), colormap=colormap,
    label=L"$C_s$ ($\mathrm{m\cdot yr^{-1}\cdot Pa^{-2}}$)",
    limits=(0, length(runs2plot)),
    ticks=(Cs_levels[1], convert_strings_to_latex(Cs_levels[2])),
    ticklabelsize=28,
    vertical=false, halign=0.7, valign=0.4
)

rowgap!(fig.layout, 1, 0.0)
rowgap!(fig.layout, 2, 0.0)
colsize!(fig.layout, 2, Relative(2 / 6))

save("figures/fig11.pdf", fig)

# ===========================================================================================================
# Figure THERM
# ===========================================================================================================

figname = "exp05_THERM_PSD"
runs2plot = range(1, Int(length(readdir("data/runs//exp05/THERM-hgeo/"))/2))

colors = collect(cgrad(cmap2, length(runs2plot), categorical=true))

xticks = Int.(-8e2:1e2:0)
xlims = (-8e2, 0)
ylims = (-600, 3000)
yticks = Int.(0:1000:3000)

fig = Figure(resolution=(1500, 1500), fonts=(; regular="TeX"), fontsize=28)

# Panel a
ax = Axis(fig[1, 1], ylabel=temp_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks_temp, convert_strings_to_latex(yticks_temp))
ylims!(ax, (-30, 15))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b)
text!(ax, -7.95e2, 10, text=a_label, align=(:left, :center))
hlines!(ax, [params.insol_ref], color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-hgeo/THERM-hgeo$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[r])
end

strt_idx = findmin(abs.(bintanja2008_T[1, :] .+ 125e3))[2]
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :], color=:black, linewidth=2.5, linestyle=:dash, label=L"Bintanja and van de Wal (2008) $\,$")
strt_idx = findmin(abs.(barker2011_temp[1, :] .+ 125e3))[2]
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :], color=:grey25, linewidth=2.5, linestyle=:dash, label=L"Barker et al. (2011) $\,$")
axislegend(ax, framevisible=false, position=:rb, labelsize=25, nbanks=2)

# Panel b
ax = Axis(fig[1, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=b_label, align=(:left, :center))
G, periods = create_PSD(df1["time"][:], df1["I"][:])
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-hgeo/THERM-hgeo$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["T"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :])
lines!(ax, periods ./ 1e3, G, color=:black, linestyle=:dash, linewidth=2.5)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :])
lines!(ax, periods ./ 1e3, G, color=:grey25, linestyle=:dash, linewidth=2.5)

# Panel c
ax = Axis(fig[2, 1], ylabel=icethick_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks, convert_strings_to_latex(yticks))
ylims!(ax, ylims)
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b, :t)
text!(ax, -7.95e2, 2500, text=L"(c)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-hgeo/THERM-hgeo$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[r])
end

# Panel d
ax = Axis(fig[2, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=d_label, align=(:left, :center))
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-hgeo/THERM-hgeo$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

# Panel e
ax = Axis(fig[3, 1], ylabel=L"$T_{ice}$ ($^o$C)", xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([-20, -10, 0], convert_strings_to_latex([-20, -10, 0]))
ylims!(ax, (-25, 6))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :t, :b)
text!(ax, -7.95e2, 4, text=L"(e)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-hgeo/THERM-hgeo$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["Tice"] .- params.degK, color=colors[r], linewidth=2)
end

# Panel f
ax = Axis(fig[3, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hidexdecorations!(ax)
hideydecorations!(ax)
hidespines!(ax, :t, :b, :l, :r)
text!(ax, 5, 0.23, text=f_label, align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-hgeo/THERM-hgeo$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["Tice"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

# Panel g
ax = Axis(fig[4, 1], ylabel=L"$f_{str}$", xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0.0, 0.5, 1.0], convert_strings_to_latex([0.0, 0.5, 1.0]))
ylims!(ax, (-0.1, 1.2))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :t)
text!(ax, -7.95e2, 1, text=L"(g)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-hgeo/THERM-hgeo$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["fstr"], color=colors[r], linewidth=2)
end

# Panel h
ax = Axis(fig[4, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :t, :l, :b, :r)
text!(ax, 5, 0.23, text=h_label, align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-hgeo/THERM-hgeo$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["fstr"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

# Colorbar(fig[4, 3], colormap=colors,
#     label=L"$h_{geo}$ (W/m²)",
#     limits=(0, length(runs2plot)),
#     ticks=(hgeo_levels[1], convert_strings_to_latex(hgeo_levels[2])),
#     ticklabelsize=35,
#     vertical=true, halign=0.5, valign=0.7
# )

Colorbar(fig[4, 1], width=Relative(1.8 / 6), height=Relative(1 / 10), colormap=cgrad(cmap2, length(runs2plot), categorical=true),
    label=L"$h_{geo}$ ($\mathrm{W \cdot m^{-2}}$)",
    limits=(0, length(runs2plot)),
    ticks=(hgeo_levels[1], convert_strings_to_latex(hgeo_levels[2])),
    ticklabelsize=28,
    vertical=false, halign=0.5, valign=0.7
)

rowgap!(fig.layout, 1, 0.0)
rowgap!(fig.layout, 2, 0.0)
rowgap!(fig.layout, 3, -90.0)

# τkin ####
runs2plot = range(1, Int(length(readdir("data/runs//exp05/THERM-taukin/"))/2))
colors2 = collect(cgrad(cmap, length(runs2plot), categorical=true))

# Panel i
ax = Axis(fig[5, 1], ylabel=temp_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks_temp, convert_strings_to_latex(yticks_temp))
ylims!(ax, (-30, 15))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b)
text!(ax, -7.95e2, 10, text=i_label, align=(:left, :center))
hlines!(ax, [params.insol_ref], color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-taukin/THERM-taukin$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors2[r])
end

strt_idx = findmin(abs.(bintanja2008_T[1, :] .+ 125e3))[2]
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :], color=:black, linewidth=2.5, linestyle=:dash, label=L"Bintanja and van de Wal (2008) $\,$")
strt_idx = findmin(abs.(barker2011_temp[1, :] .+ 125e3))[2]
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :], color=:grey25, linewidth=2.5, linestyle=:dash, label=L"Barker et al. (2011) $\,$")
axislegend(ax, framevisible=false, position=:rb, labelsize=25, nbanks=2)

# Panel j
ax = Axis(fig[5, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=j_label, align=(:left, :center))
G, periods = create_PSD(df1["time"][:], df1["I"][:])
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-taukin/THERM-taukin$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["T"][:])
    lines!(ax, periods ./ 1e3, G, color=colors2[r], linewidth=2)
end

G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :])
lines!(ax, periods ./ 1e3, G, color=:black, linestyle=:dash, linewidth=2.5)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :])
lines!(ax, periods ./ 1e3, G, color=:grey25, linestyle=:dash, linewidth=2.5)

# Panel k
ax = Axis(fig[6, 1], ylabel=icethick_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks, convert_strings_to_latex(yticks))
ylims!(ax, ylims)
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b, :t)
text!(ax, -7.95e2, 2500, text=L"(k)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-taukin/THERM-taukin$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors2[r])
end

# Panel l
ax = Axis(fig[6, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=l_label, align=(:left, :center))
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-taukin/THERM-taukin$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:])
    lines!(ax, periods ./ 1e3, G, color=colors2[r], linewidth=2)
end

# Panel m
ax = Axis(fig[7, 1], ylabel=L"$T_{ice}$ ($^o$C)", xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([-20, -10, 0], convert_strings_to_latex([-20, -10, 0]))
ylims!(ax, (-25, 6))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :t, :b)
text!(ax, -7.95e2, 4, text=L"(m)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-taukin/THERM-taukin$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["Tice"] .- params.degK, color=colors2[r], linewidth=2)
end

# Panel n
ax = Axis(fig[7, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hidexdecorations!(ax)
hideydecorations!(ax)
hidespines!(ax, :t, :b, :l, :r)
text!(ax, 5, 0.23, text=L"(n)$\,$", align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-taukin/THERM-taukin$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["Tice"][:])
    lines!(ax, periods ./ 1e3, G, color=colors2[r], linewidth=2)
end

# Panel o
ax = Axis(fig[8, 1], ylabel=L"$f_{str}$", xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0.0, 0.5, 1.0], convert_strings_to_latex([0.0, 0.5, 1.0]))
ylims!(ax, (-0.1, 1.2))
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -7.95e2, 1, text=L"(o)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-taukin/THERM-taukin$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["fstr"], color=colors2[r], linewidth=2)
end

# Panel p
ax = Axis(fig[8, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.23, text=L"(p)$\,$", align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-taukin/THERM-taukin$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["fstr"][:])
    lines!(ax, periods ./ 1e3, G, color=colors2[r], linewidth=2)
end

# Colorbar(fig[5, 3], colormap=colors2,
#     label=L"$\tau_{kin} (yr)$",
#     limits=(0, length(runs2plot)),
#     ticks=(taukin_levels[1], convert_strings_to_latex(taukin_levels[2])),
#     ticklabelsize=35,
#     vertical=true, halign=0.5, valign=0.7
# )

Colorbar(fig[8, 1], width=Relative(1.8 / 6), height=Relative(1 / 10), colormap=cgrad(cmap, length(runs2plot), categorical=true),
    label=L"$\tau_{kin}$ (yr)",
    limits=(0, length(runs2plot)),
    ticks=(taukin_levels[1], convert_strings_to_latex(taukin_levels[2])),
    ticklabelsize=28,
    vertical=false, halign=0.5, valign=0.7
)

rowgap!(fig.layout, 5, 0.0)
rowgap!(fig.layout, 6, 0.0)
rowgap!(fig.layout, 7, -90.0)

colsize!(fig.layout, 2, Relative(2 / 6))

save("figures/fig13.pdf", fig)

# ===========================================================================================================
# Figure AGING-τalpha
# ===========================================================================================================

figname = "exp06_AGING-taualpha_PSD"
runs2plot = range(1, 11)
df1 = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW1.nc")
params = JLD2.load_object("data/runs//exp06/AGING-SLOW/AGING-SLOW1_params.jld2")

colormap = cgrad(cmap, length(runs2plot), categorical=true)
colors = collect(colormap)

xticks = Int.(-8e2:1e2:0)
xlims = (-8e2, 0)
ylims = (-600, 3000)
yticks = Int.(0:1000:3000)

fig = Figure(resolution=(1500, 750), fonts=(; regular="TeX"), fontsize=28)

# Panel a
ax = Axis(fig[1, 1], ylabel=temp_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks_temp, convert_strings_to_latex(yticks_temp))
ylims!(ax, (-30, 15))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b)
text!(ax, -7.95e2, 10, text=a_label, align=(:left, :center))
hlines!(ax, [params.insol_ref], color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[r])
end

dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[1]).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[1], linewidth=2)

strt_idx = findmin(abs.(bintanja2008_T[1, :] .+ 125e3))[2]
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :], color=:black, linewidth=2.5, linestyle=:dash, label=L"Bintanja and van de Wal (2008) $\,$")
strt_idx = findmin(abs.(barker2011_temp[1, :] .+ 125e3))[2]
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :], color=:grey25, linewidth=2.5, linestyle=:dash, label=L"Barker et al. (2011) $\,$")
axislegend(ax, framevisible=false, position=:rt, labelsize=25, nbanks=2)

# Panel b
ax = Axis(fig[1, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=b_label, align=(:left, :center))
G, periods = create_PSD(df1["time"][:], df1["I"][:])
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["T"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[1]).nc")
G, periods = create_PSD(dfr["time"][:], dfr["T"][:])
lines!(ax, periods ./ 1e3, G, color=colors[1], linewidth=2.5)

G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :])
lines!(ax, periods ./ 1e3, G, color=:black, linestyle=:dash, linewidth=2.5)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :])
lines!(ax, periods ./ 1e3, G, color=:grey25, linestyle=:dash, linewidth=2.5)

# Panel c
ax = Axis(fig[2, 1], ylabel=icethick_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks, convert_strings_to_latex(yticks))
ylims!(ax, ylims)
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b, :t)
text!(ax, -7.95e2, 2500, text=L"c)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[r])
end

dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[1]).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[1], linewidth=2)

# Panel d
ax = Axis(fig[2, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=d_label, align=(:left, :center))
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[1]).nc")
G, periods = create_PSD(dfr["time"][:], dfr["H"][:])
lines!(ax, periods ./ 1e3, G, color=colors[1], linewidth=2.5)

# Panel e
ax = Axis(fig[3, 1], ylabel=L"$\alpha$", xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0.2, 0.5, 0.9], convert_strings_to_latex([0.2, 0.5, 0.9]))
ylims!(ax, (0.0, 1.2))
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -7.95e2, 1, text=L"e)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[r]).nc")
    params = JLD2.load_object("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[r])_params.jld2")
    #alpha_eff = dfr["albedo"][:]
    #alpha_eff[dfr["m"][:] .> 0] .= NaN
    lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=colors[r], linewidth=2)
    #lines!(ax, dfr["time"] ./ 1e3, alpha_eff, color=colors[r], linewidth=2)
end

dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[1]).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=colors[1], linewidth=2)

#lines!(ax, df1["time"] ./ 1e3, df1["albedo"] .* NaN, linestyle=:dot, color=pacco_color, linewidth=4, label=L"\alpha")
#lines!(ax, df1["time"] ./ 1e3, df1["albedo"] .* NaN, color=pacco_color, linewidth=4, label=L"\alpha_\mathrm{eff}")
#axislegend(ax, framevisible=false, position=:rt, labelsize=30, nbanks=2)

# Panel f
ax = Axis(fig[3, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.23, text=f_label, align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

dfr = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(runs2plot[1]).nc")
G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:])
lines!(ax, periods ./ 1e3, G, color=colors[1], linewidth=2.5)

Colorbar(fig[2, 1], width=Relative(1.8 / 6), height=Relative(1 / 10), colormap=colormap,
    label=L"$\tau_{\alpha}$ (yr)",
    limits=(0, length(runs2plot)),
    ticks=(taualpha_levels[1], convert_strings_to_latex(taualpha_levels[2])),
    ticklabelsize=28,
    vertical=false, halign=0.5, valign=0.85 
)

rowgap!(fig.layout, 1, -90.0)
rowgap!(fig.layout, 2, 0.0)
colsize!(fig.layout, 2, Relative(2 / 6))

save("figures/fig15.pdf", fig)

# ===========================================================================================================
# Figure 10. AGING (best τα)
# ===========================================================================================================
figname = "exp06_AGING_best"
best_exp = 1
df = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(best_exp).nc")
params = JLD2.load_object("data/runs//exp06/AGING-SLOW/AGING-SLOW$(best_exp)_params.jld2")

# Interpole data for correlation:

# # Temperature
# data = [df["time"][:], bintanja2008_T[1, :], barker2011_temp[1, :]]
# tmin = argmin([size(tvec) for tvec in data])   
# #println(tmin)

# t1, x1 = df["time"][:], df["T"][:] .- df["Tref"][:]
# t2, x2 = bintanja2008_T[1, :], bintanja2008_T[2, :]
# t3, x3 = barker2011_temp[1, :], barker2011_temp[2, :]

# tmin = maximum([minimum(tvec) for tvec in data])
# tmax = minimum([maximum(tvec) for tvec in data])
# imin = argmin(abs.(t1 .- tmin))
# imax = argmin(abs.(t1 .- tmax))
# t = t1[imin:imax-1]

# x1itp = linear_interpolation(t1, x1)
# x2itp = linear_interpolation(t2, x2)
# x3itp = linear_interpolation(t3, x3)

# x1_new = x1itp.(t)
# x2_new = x2itp.(t)
# x3_new = x3itp.(t)

# print("AGING best experiment metrics")
# print("Temperature rmse --> $(rmsd(x1_new, x2_new)), $(rmsd(x1_new, x3_new))")
# print("Temperature corr --> $(cor(x1_new, x2_new)), $(cor(x1_new, x3_new))")

# # CO2
# data = [df["time"][:], luthi2008_co2[1, :], berends2021_co2[1, :], yamamoto2022_co2[1, :]]
# tmin = argmin([size(tvec) for tvec in data])   
# #println(tmin)

# t1, x1 = df["time"][:], df["T"][:] .- df["Tref"][:]
# t2, x2 = luthi2008_co2[1, :], luthi2008_co2[2, :]
# t3, x3 = berends2021_co2[1, :], berends2021_co2[2, :]
# t4, x4 = yamamoto2022_co2[1, :], yamamoto2022_co2[2, :]

# tmin = maximum([minimum(tvec) for tvec in data])
# tmax = minimum([maximum(tvec) for tvec in data])
# imin = argmin(abs.(t1 .- tmin))
# imax = argmin(abs.(t1 .- tmax))
# t = t1[imin:imax-1]

# x1itp = linear_interpolation(t1, x1)
# x2itp = linear_interpolation(t2, x2)
# x3itp = linear_interpolation(t3, x3)
# x4itp = linear_interpolation(t4, x4)

# x1_new = x1itp.(t)
# x2_new = x2itp.(t)
# x3_new = x3itp.(t)
# x4_new = x4itp.(t)

# print("CO2 rmse --> $(rmsd(x1_new, x2_new)), $(rmsd(x1_new, x3_new)), $(rmsd(x1_new, x4_new))")
# print("CO2 corr --> $(cor(x1_new, x2_new)), $(cor(x1_new, x3_new)), $(cor(x1_new, x4_new))")

# # Volume
# data = [df["time"][:], bintanja2008_vol[1, :], berends2021_vol[1, :], spratt2016_vol[1, :]]
# tmin = argmin([size(tvec) for tvec in data])   
# #println(tmin)

# t1, x1 = df["time"][:], df["T"][:] .- df["Tref"][:]
# t2, x2 = bintanja2008_vol[1, :], bintanja2008_vol[2, :]
# t3, x3 = berends2021_vol[1, :], berends2021_vol[2, :]
# t4, x4 = spratt2016_vol[1, :], spratt2016_vol[2, :]

# tmin = maximum([minimum(tvec) for tvec in data])
# tmax = minimum([maximum(tvec) for tvec in data])
# imin = argmin(abs.(t4 .- tmin))
# imax = argmin(abs.(t4 .- tmax))
# t = t4[imin:imax-1]

# x1itp = linear_interpolation(t1, x1)
# x2itp = linear_interpolation(t2, x2)
# x3itp = linear_interpolation(t3, x3)
# x4itp = linear_interpolation(t4, x4)

# x1_new = x1itp.(t)
# x2_new = x2itp.(t)
# x3_new = x3itp.(t)
# x4_new = x4itp.(t)

# print("Volume rmse --> $(rmsd(x1_new, x2_new)), $(rmsd(x1_new, x3_new)), $(rmsd(x1_new, x4_new))")
# print("Volume corr --> $(cor(x1_new, x2_new)), $(cor(x1_new, x3_new)), $(cor(x1_new, x4_new))")

# Plotting
fig = Figure(resolution=(1500, 1200), fonts=(; regular="TeX"), fontsize=28)
linewidth = 3

# Panel a. Insolation forcing
ax = Axis(fig[1, 1], ylabel=ins_label, xgridvisible=false, ygridvisible=false, yaxisposition=:left)
ax.xticks = xticks_time
ax.yticks = (yticks_insol, convert_strings_to_latex(yticks_insol))
ylims!(ax, ylims_ins)
xlims!(ax, xlims_time_800)
hidexdecorations!(ax)
hidespines!(ax, :b)
text!(ax, -792, 560, text=a_label, align=(:left, :center))
hlines!(ax, [params.insol_ref], color=background_color, linestyle=:solid)
lines!(ax, df["time"] ./ 1e3, df["I"], color=laskar2004color, linewidth=linewidth, label=L"Laskar (2004)$\,$")

# Panel b. Insolation PSD
ax = Axis(fig[1, 2], xgridvisible=false, ygridvisible=false, yaxisposition=:left)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=b_label, align=(:left, :center))
G, periods = create_PSD(cut_run(df["time"][:], df["time"][:], -8e5), cut_run(df["time"][:], df["I"][:], -8e5))
lines!(ax, periods ./ 1e3, G, color=laskar2004color, linewidth=linewidth)

# Panel c. Climatic Temperature
ax = Axis(fig[2, 1], ylabel=temp_label, xlabel=time_label, xgridvisible=false, ygridvisible=false, yaxisposition=:left)
ax.yticks = (yticks_temp, convert_strings_to_latex(yticks_temp))
ax.xticks = xticks_time
ylims!(ax, ylims_temp)
xlims!(ax, xlims_time_800)
hidexdecorations!(ax)
hidespines!(ax, :b, :t)
text!(ax, -792, 10, text=c_label, align=(:left, :center))
hlines!(ax, [0.0], color=background_color, linestyle=:solid)
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :], linewidth=3, color=bintanja2008color, label=L"Bintanja and van de Wal (2008)$\,$")
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :], linewidth=3, color=barker2011color, label=L"Barker et al. (2011)$\,$")
lines!(ax, df["time"] ./ 1e3, df["T"] .- df["Tref"], color=pacco_color, linewidth=linewidth)
axislegend(ax, framevisible=false, position=(0.15, 1.0), labelsize=25, nbanks=2)

display("Temperature correlations:")


# Panel d. Temperature PSD
ax = Axis(fig[2, 2], xgridvisible=false, ygridvisible=false, yaxisposition=:left)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=d_label, align=(:left, :center))
G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :])
lines!(ax, periods ./ 1e3, G, color=bintanja2008color, linewidth=linewidth)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :])
lines!(ax, periods ./ 1e3, G, color=barker2011color, linewidth=linewidth)
G, periods = create_PSD(cut_run(df["time"][:], df["time"][:], -8e5), cut_run(df["time"][:], df["T"][:], -8e5))
lines!(ax, periods ./ 1e3, G, color=pacco_color, linewidth=linewidth)

# Panel e. Carbon dioxide
ax = Axis(fig[3, 1], ylabel=carbon_label, xgridvisible=false, ygridvisible=false, yaxisposition=:left)
ax.yticks = (yticks_carbon, convert_strings_to_latex(yticks_carbon))
ax.xticks = xticks_time
xlims!(ax, xlims_time_800)
hidexdecorations!(ax)
hidespines!(ax, :b, :t)
hlines!(ax, [params.Cref], color=background_color, linestyle=:solid)
text!(ax, -792, 320, text=e_label, align=(:left, :center))
lines!(ax, luthi2008_co2[1, :] ./ 1e3, luthi2008_co2[2, :], linewidth=3, color=luthi2008color, label=L"Lüthi et al. (2008)$\,$")
lines!(ax, berends2021_co2[1, :] ./ 1e3, berends2021_co2[2, :], linewidth=3, color=berends2021color, label=L"Berends et al. (2021b)$\,$")
lines!(ax, yamamoto2022_co2[1, :] ./ 1e3, yamamoto2022_co2[2, :], linewidth=3, color=yamamoto2022color, label=L"Yamamoto et al. (2022)$\,$")
lines!(ax, df["time"] ./ 1e3, df["C"], color=pacco_color, linewidth=linewidth)
axislegend(ax, framevisible=false, position=(0.8, 1.1), labelsize=25, nbanks=3)

# Panel f. Carbon dioxide PSD
ax = Axis(fig[3, 2], xgridvisible=false, ygridvisible=false, yaxisposition=:left)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=f_label, align=(:left, :center))
G, periods = create_PSD(luthi2008_co2[1, :], luthi2008_co2[2, :])
lines!(ax, periods ./ 1e3, G, color=luthi2008color, linewidth=linewidth)
G, periods = create_PSD(berends2021_co2[1, :], berends2021_co2[2, :])
lines!(ax, periods ./ 1e3, G, color=berends2021color, linewidth=linewidth)
G, periods = create_PSD(yamamoto2022_co2[1, :], yamamoto2022_co2[2, :])
lines!(ax, periods ./ 1e3, G, color=yamamoto2022color, linewidth=linewidth)
G, periods = create_PSD(cut_run(df["time"][:], df["time"][:], -8e5), cut_run(df["time"][:], df["C"][:], -8e5))
lines!(ax, periods ./ 1e3, G, color=pacco_color, linewidth=linewidth)

# Panel g. Ice thickness
ax = Axis(fig[4, 1], ylabel=icethick_label, xgridvisible=false, ygridvisible=false, yaxisposition=:left)
ax.yticks = (yticks_ice, convert_strings_to_latex(yticks_ice))
ax.xticks = xticks_time
ylims!(ax, ylims_ice)
xlims!(ax, xlims_time_800)
hidexdecorations!(ax)
hidespines!(ax, :b, :t)
hlines!(ax, [0.0], color=background_color, linestyle=:solid)
text!(ax, -792, 2000, text=g_label, align=(:left, :center))
lines!(ax, df["time"] ./ 1e3, df["H"], color=pacco_color, linewidth=linewidth, label="Ice-sheet thickness, H")

# Panel h. Ice thickness PSD
ax = Axis(fig[4, 2], xgridvisible=false, ygridvisible=false, yaxisposition=:left)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=h_label, align=(:left, :center))
G, periods = create_PSD(cut_run(df["time"][:], df["time"][:], -8e5), cut_run(df["time"][:], df["H"][:], -8e5))
lines!(ax, periods ./ 1e3, G, color=pacco_color, linewidth=linewidth)

# Panel i. Ice volume
ax = Axis(fig[5, 1], ylabel=icevol_label, xlabel=time_label, xgridvisible=false, ygridvisible=false, yaxisposition=:left)
ax.yticks = (yticks_vol, convert_strings_to_latex(yticks_vol))
ax.xticks = xticks_time
ylims!(ax, ylims_vol)
xlims!(ax, xlims_time_800)
hidespines!(ax, :t)
hlines!(ax, [0.0], color=background_color, linestyle=:solid)
text!(ax, -792, 25, text=i_label, align=(:left, :center))
lines!(ax, df["time"] ./ 1e3, df["Vol"], color=pacco_color, linewidth=linewidth)
lines!(ax, bintanja2008_vol[1, :] ./ 1e3, bintanja2008_vol[2, :], linewidth=3, color=bintanja2008color, label=L"Bintanja and van de Wal (2008)$\,$")
lines!(ax, berends2021_vol[1, :] ./ 1e3, berends2021_vol[2, :], linewidth=3, color=berends2021color, label=L"Berends et al. (2021b)$\,$")
lines!(ax, spratt2016_vol[1, :] ./ 1e3, spratt2016_vol[2, :], linewidth=3, color=spratt2016color, label=L"Spratt and Lisiecki (2016)$\,$")
axislegend(ax, framevisible=false, position=(0.15, 1.0), labelsize=25, nbanks=2)

# Panel j. Ice volume PSD
ax = Axis(fig[5, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false, yaxisposition=:left)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.23, text=j_label, align=(:left, :center))
G, periods = create_PSD(bintanja2008_vol[1, :], bintanja2008_vol[2, :])
lines!(ax, periods ./ 1e3, G, color=bintanja2008color, linewidth=linewidth)
G, periods = create_PSD(berends2021_vol[1, :], berends2021_vol[2, :])
lines!(ax, periods ./ 1e3, G, color=berends2021color, linewidth=linewidth)
G, periods = create_PSD(spratt2016_vol[1, :], spratt2016_vol[2, :])
lines!(ax, periods ./ 1e3, G, color=spratt2016color, linewidth=linewidth)
G, periods = create_PSD(cut_run(df["time"][:], df["time"][:], -8e5), cut_run(df["time"][:], df["Vol"][:], -8e5))
lines!(ax, periods ./ 1e3, G, color=pacco_color, linewidth=linewidth)

rowgap!(fig.layout, 0.0)
colsize!(fig.layout, 2, Relative(2 / 6))
save("figures/fig16.pdf", fig)

# ===========================================================================================================
# Figure 11. AGING (best) thresholds
# ===========================================================================================================
figname = "exp06_AGING_best_thr"
best_exp = 1
df = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(best_exp).nc")
params = JLD2.load_object("data/runs//exp06/AGING-SLOW/AGING-SLOW$(best_exp)_params.jld2")

albedo_eff = df["albedo"][:]
albedo_eff[df["m"][:].>0] .= params.albedo_newice

absortion = 1 .- albedo_eff
sw = absortion .* (df["I"] .- params.insol_ref)
sw_eff = copy(sw)
sw_eff[df["I"][:].<=params.insol_ref] .= 0.0
sw_eff[df["H"][:].<=10.0] .= 0.0

xticks = Int.(-8e2:0.5e2:0)
xlims = (-3.5e2, 0)
ylims = (-2000, 3500)
yticks = Int.(-1500:1500:3000)

fig = Figure(resolution=(750, 750), fonts=(; regular="TeX"), fontsize=32)

# Panel a
ax = Axis(fig[1, 1], ylabel=L"$I - I_\mathrm{ref}$ ($\mathrm{W \cdot m^{-2}}$)", xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([-50, 0, 50], convert_strings_to_latex([-50, 0, 50]))
xlims!(ax, xlims)
ylims!(ax, (-80, 100))
hidespines!(ax, :b)
hidexdecorations!(ax)
text!(ax, -3.45e2, 80, text=L"(a)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
# band!(ax, df["time"] ./ 1e3, 0.0, sw_eff, color=:royalblue, label=L"SW_\mathrm{eff}")
lines!(ax, df["time"] ./ 1e3, df["I"] .- params.insol_ref, color=:black, linewidth=3)
# axislegend(ax, framevisible=false, position=:rt, labelsize=25, nbanks=3)

# Panel b
ax = Axis(fig[2, 1], ylabel=L"$I_\mathrm{eff}$ ($\mathrm{W \cdot m^{-2}}$)", xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 10, 20], convert_strings_to_latex([0, 10, 20]))
xlims!(ax, xlims)
ylims!(ax, (-5, 35))
hidespines!(ax, :b, :t)
hidexdecorations!(ax)
text!(ax, -3.45e2, 30, text=L"(b)$\,$", align=(:left, :center))
lines!(ax, df["time"] ./ 1e3, sw_eff, color=:black, linewidth=3)

# Panel c
ax = Axis(fig[3, 1], xlabel=time_label, ylabel=icethick_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks, convert_strings_to_latex(yticks))
ylims!(ax, ylims)
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -3.45e2, 2500, text=L"(c)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
lines!(ax, df["time"] ./ 1e3, df["H"], color=:black, linewidth=3)

rowgap!(fig.layout, 0.0)
save("figures/figA02.pdf", fig)

# ===========================================================================================================
# Figure 13. AGING phase space (H)
# ===========================================================================================================
var2plot="H"
figname = "exp06_AGING_PhaseSpace_$(var2plot)"
best_exp = 1
df = NCDataset("data/runs//exp06/AGING-SLOW/AGING-SLOW$(best_exp).nc")
params = JLD2.load_object("data/runs//exp06/AGING-SLOW/AGING-SLOW$(best_exp)_params.jld2")
Ianom = df["I"] .- params.insol_ref
t1, t2 = 120, 2

albedo_eff = df["albedo"][:]
albedo_eff[df["m"][:].>0] .= params.albedo_newice

absortion = 1 .- albedo_eff
sw = absortion .* (df["I"] .- params.insol_ref)
sw_eff = copy(sw)


fig = Figure(resolution=(600, 500), fonts=(; regular="TeX"), fontsize=24)
# Panel a
ax = Axis(fig[1, 1], ylabel= icethick_label, xlabel=L"$I - I_\mathrm{ref}$ ($\mathrm{W \cdot m^{-2}}$)", xgridvisible=false, ygridvisible=false)
ax.xticks = ([-50, 0, 50], convert_strings_to_latex([-50, 0, 50]))
ax.yticks = (yticks_ice, convert_strings_to_latex(yticks_ice))
ylims!(ax, ylims_ice)
xlims!(ax, (-60, 80))
vlines!(ax, 0.0, color=background_color, linestyle=:solid)

start_lgc = -Int(t1 - length(df["time"]))
end_lgc = -Int(t2 - length(df["time"]))
lines!(ax, Ianom, df["H"], color=:grey, linewidth=1)
#scatter!(ax, sw_eff[start_lgc:end_lgc], df["H"][start_lgc:end_lgc], color=:royalblue, linewidth=3)
# scatterlines!(ax, Ianom[start_lgc:end_lgc], df["H"][start_lgc:end_lgc], color=:black, linewidth=3)
lines!(ax, Ianom[start_lgc:end_lgc], df["H"][start_lgc:end_lgc], color=:black, linewidth=3)

rowgap!(fig.layout, 0.0)
save("figures/fig17.pdf", fig)

# ===========================================================================================================
# Figure AGING-Cs
# ===========================================================================================================

figname = "exp06_AGING-Cs_PSD"
runs2plot = range(1, 23)
df1 = NCDataset("data/runs//exp06/AGING-Cs/AGING-Cs1.nc")
params = JLD2.load_object("data/runs//exp06/AGING-Cs/AGING-Cs1_params.jld2")

colormap = cgrad(cmap, length(runs2plot), categorical=true)
colors = collect(colormap)

xticks = Int.(-8e2:1e2:0)
xlims = (-8e2, 0)
ylims = (-600, 3000)
yticks = Int.(0:1000:3000)

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
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=b_label, align=(:left, :center))
G, periods = create_PSD(df1["time"][:], df1["I"][:])
lines!(ax, periods ./ 1e3, G, color=pacco_color, linewidth=3)

# Panel c
ax = Axis(fig[2, 1], ylabel=icethick_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks, convert_strings_to_latex(yticks))
ylims!(ax, ylims)
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :b, :t)
text!(ax, -7.95e2, 2500, text=L"(c)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-Cs/AGING-Cs$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[r])
end

# Panel d
ax = Axis(fig[2, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=d_label, align=(:left, :center))
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-Cs/AGING-Cs$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

# Panel e
ax = Axis(fig[3, 1], ylabel=temp_label, xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = (yticks_temp, convert_strings_to_latex(yticks_temp))
ylims!(ax, ylims_temp)
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -7.95e2, 10, text=L"(e)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-Cs/AGING-Cs$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[r], linewidth=2)
end

strt_idx = findmin(abs.(bintanja2008_T[1, :] .+ 125e3))[2]
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :], color=:black, linewidth=2.5, linestyle=:dash, label=L"Bintanja and van de Wal (2008) $\,$")
strt_idx = findmin(abs.(barker2011_temp[1, :] .+ 125e3))[2]
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :], color=:grey25, linewidth=2.5, linestyle=:dash, label=L"Barker et al. (2011) $\,$")
axislegend(ax, framevisible=false, position=:rt, labelsize=25, nbanks=2)

# Panel f
ax = Axis(fig[3, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.23, text=f_label, align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-Cs/AGING-Cs$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["T"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :])
lines!(ax, periods ./ 1e3, G, color=:black, linestyle=:dash, linewidth=2.5)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :])
lines!(ax, periods ./ 1e3, G, color=:grey25, linestyle=:dash, linewidth=2.5)

Colorbar(fig[1, 2], width=Relative(1.2 / 3), height=Relative(1 / 10), colormap=colormap,
    label=L"$C_s$ ($\mathrm{m\cdot yr^{-1}\cdot Pa^{-2}}$)",
    limits=(0, length(runs2plot)),
    ticks=(Cs_levels[1], convert_strings_to_latex(Cs_levels[2])),
    ticklabelsize=28,
    vertical=false, halign=0.7, valign=0.4
)

rowgap!(fig.layout, 1, 0.0)
rowgap!(fig.layout, 2, 0.0)
colsize!(fig.layout, 2, Relative(2 / 6))

save("figures/figA01.pdf", fig)