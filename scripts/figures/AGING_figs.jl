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
taualpha_levels_2 = [[1, 2, 3, 4, 5] .- 0.5,
                    ["1", "10", "100", "130", "1000"]] # [1e3, 1e4, 1e5, 1.3e5, 1e6] 
hgeo_levels = [[1, 9, 19], ["10^{-3}", L"5\cdot10^{-3}", "10^{-2}"]] # collect(1:0.5:10) .* 1e-3

# Ticks
xticks_time = (xticks_time_800, convert_strings_to_latex(-1 .* xticks_time_800))

function calc_zthr(Tsl::Vector, Tref::Vector, p)
    return (p.sref .+ p.ks .* (Tsl .- Tref) .- p.lambda .* (Tsl .- p.Tthreshold)) / ((p.ks .- p.lambda) .* p.Γ)
end

# This code produce figures from AGING experiment after Andrey Ganopolski's comments
# and our decission to change albedo equation to one more simplistic and clear

exp_string = "/AGING-DIAG/AGING-DIAG"
selected_run=4

# ===========================================================================================================
# AGING-τalpha DIAG
# ===========================================================================================================
figname = "exp06_AGING-taualpha_PSD_2"

runs2plot = range(1, 5)
df1 = NCDataset("data/runs//exp06$(exp_string)1.nc")
params = JLD2.load_object("data/runs//exp06$(exp_string)1_params.jld2")
cmap = [:darkslateblue, :steelblue, :skyblue2, :black, :darkorange]
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
    dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[r])
end

dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[selected_run]).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[selected_run], linewidth=2)

strt_idx = findmin(abs.(bintanja2008_T[1, :] .+ 125e3))[2]
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :], color=:black, linewidth=2.5, linestyle=:dash, label=L"Bintanja and van de Wal (2008) $\,$")
strt_idx = findmin(abs.(barker2011_temp[1, :] .+ 125e3))[2]
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :], color=:grey25, linewidth=2.5, linestyle=:dash, label=L"Barker et al. (2011) $\,$")

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
    dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["T"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[selected_run]).nc")
G, periods = create_PSD(dfr["time"][:], dfr["T"][:])
lines!(ax, periods ./ 1e3, G, color=colors[selected_run], linewidth=2.5)

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
    dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[r]).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[r])
end

dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[selected_run]).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[selected_run], linewidth=2)

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
    dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[selected_run]).nc")
G, periods = create_PSD(dfr["time"][:], dfr["H"][:])
lines!(ax, periods ./ 1e3, G, color=colors[selected_run], linewidth=2.5)

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
    dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[r]).nc")
    params = JLD2.load_object("data/runs//exp06$(exp_string)$(runs2plot[r])_params.jld2")
    #alpha_eff = dfr["albedo"][:]
    #alpha_eff[dfr["m"][:] .> 0] .= NaN
    lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=colors[r], linewidth=2)
    #lines!(ax, dfr["time"] ./ 1e3, alpha_eff, color=colors[r], linewidth=2)
end

dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[selected_run]).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=colors[selected_run], linewidth=2)

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
    dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[r]).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:])
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

dfr = NCDataset("data/runs//exp06$(exp_string)$(runs2plot[selected_run]).nc")
G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:])
lines!(ax, periods ./ 1e3, G, color=colors[selected_run], linewidth=2.5)

ax = Axis(fig[0, 1])
hidedecorations!(ax)
hidespines!(ax)
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :] .* NaN, color=:black, linewidth=2.5, linestyle=:dash, label=L"Bintanja and van de Wal (2008) $\,$")
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :] .* NaN, color=:grey25, linewidth=2.5, linestyle=:dash, label=L"Barker et al. (2011) $\,$")
axislegend(ax, framevisible=false, position=(:center, :center), labelsize=28, nbanks=2)

Colorbar(fig[0, 2], width=Relative(2 / 3), colormap=colormap,
    label=L"$\tau_{\alpha}$ (kyr)",
    limits=(0, length(runs2plot)),
    ticks=(taualpha_levels_2[1], convert_strings_to_latex(taualpha_levels_2[2])),
    ticklabelsize=28,
    vertical=false)

rowgap!(fig.layout, 0.0)
colsize!(fig.layout, 2, Relative(2 / 6))
rowsize!(fig.layout, 0, Relative(0.1))
save("figures/fig15.pdf", fig)

# ===========================================================================================================
# AGING (best τα)
# ===========================================================================================================
figname = "exp06_AGING_best"
df = NCDataset("data/runs//exp06$(exp_string)$(selected_run).nc")
params = JLD2.load_object("data/runs//exp06$(exp_string)$(selected_run)_params.jld2")

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
# AGING phase space (H)
# ===========================================================================================================
var2plot="H"
figname = "exp06_AGING_PhaseSpace_$(var2plot)"
df = NCDataset("data/runs//exp06$(exp_string)$(selected_run).nc")
params = JLD2.load_object("data/runs//exp06$(exp_string)$(selected_run)_params.jld2")
t1, t2 = 120, 2
Ianom = df["I"] .- params.insol_ref

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

cmap = [:black, :purple, :royalblue, :olive, :darkorange, :darkred]
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

# ===========================================================================================================
# AGING (best) thresholds
# ===========================================================================================================
figname = "exp06_AGING_best_thr"
df = NCDataset("data/runs//exp06$(exp_string)$(selected_run).nc")
params = JLD2.load_object("data/runs//exp06$(exp_string)$(selected_run)_params.jld2")

Ieff = calc_Ieff(df, params)

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
lines!(ax, df["time"] ./ 1e3, Ieff, color=:black, linewidth=3)

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