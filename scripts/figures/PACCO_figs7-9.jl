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

Css = vcat(1e-10:1e-10:9e-10, 1e-9:1e-9:9e-9, 1e-8:1e-8:9e-8, 1e-7:1e-7:9e-7, 1e-6:1e-6:9e-6, 1e-5:1e-5:9e-5, 1e-4) 
taukins = [1e3, 5e3, 1e4, 2.5e4, 5e4, 7.5e4, 1e5, 2.5e5, 5e5, 7.5e5, 1e6]
taualphas = 10e3:5e3:120e3 # vcat(1e3, 10e3:5e3:120e3, 400e3, 1e6, 1e7)
Cs_levels = [[1, 27, length(Css)], ["10^{-10}", "10^{-7}", "10^{-4}"]]
taukin_levels = [[1, 3, 7, 11], ["1", "10^{1}", "10^{2}", "10^{3}"]]  # [10, 50, 100, 500, 1000, 2500, 5000, 7500, 10000]
taualpha_levels = [[1, 2, 12, 24, 25, 26, 27], ["1", "10", "60", "120", "400", "10^{3}", "10^{4}"]]   # tau_alpha = vcat(1e3, 10e3:5e3:120e3, 400e3, 1e6, 1e7)
cmap = [:purple, :royalblue, :olive, :darkorange, :darkred]
cmap2 = [:lightblue, :royalblue, :black, :goldenrod, :darkorange, :darkred]

# Ticks
xticks_time = (xticks_time_800, convert_strings_to_latex(-1 .* xticks_time_800))

function calc_zthr(Tsl::Vector, Tref::Vector, p)
    return (p.sref .+ p.ks .* (Tsl .- Tref) .- p.lambda .* (Tsl .- p.Tthreshold)) / ((p.ks .- p.lambda) .* p.Γ)
end

# ===========================================================================================================
# Figure BASE-MV-Cs
# ===========================================================================================================

figname = "exp04_BASE-MV-Cs_PSD"
runs2plot = range(1, length(Css))
df1 = NCDataset("data/runs//exp04/BASE-MV-Cs/BASE-MV-Cs01.nc")
params = JLD2.load_object("data/runs//exp04/BASE-MV-Cs/BASE-MV-Cs01_params.jld2")

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
text!(ax, -7.95e2, 2500, text=L"(c)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp04/BASE-MV-Cs/BASE-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
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
    dfr = NCDataset("data/runs//exp04/BASE-MV-Cs/BASE-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
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
    dfr = NCDataset("data/runs//exp04/BASE-MV-Cs/BASE-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
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
    dfr = NCDataset("data/runs//exp04/BASE-MV-Cs/BASE-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["T"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=:black, linestyle=:dash, linewidth=3.5)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=:grey25, linestyle=:dash, linewidth=3.5)

Colorbar(fig[2:3, 3], width=Relative(0.3), height=Relative(1 / 3), colormap=colormap,
    label = L"$C_s$ ($\mathrm{m\cdot yr^{-1}\cdot Pa^{-2}}$)",
    limits=(1, length(runs2plot)+1),
    ticks=(Cs_levels[1] .+ 0.5, convert_strings_to_latex(Cs_levels[2])),
    ticklabelsize=28,
    vertical=true, halign=0.0,
)

rowgap!(fig.layout, 0.0)
colsize!(fig.layout, 2, Relative(2 / 6))
colsize!(fig.layout, 3, Relative(0.04))
resize_to_layout!(fig)

save("figures/fig11.pdf", fig)

# ===========================================================================================================
# Figure THERM-taukin
# ===========================================================================================================
figname = "exp05_THERM-taukin_PSD"
runs2plot = range(1, length(taukins))
cmap13 = [:plum3, :purple, :slategray2, :lightsteelblue3, :cornflowerblue, :royalblue3, :royalblue4, :orange, :darkorange, :darkorange2, :darkorange4]
colors13 = collect(cgrad(cmap13, length(runs2plot), categorical=true))
xticks = Int.(-8e2:1e2:0)
xlims = (-8e2, 0)
ylims = (-600, 3000)
yticks = [0, 1000, 2000]

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
    dfr = NCDataset("data/runs//exp05/THERM-MV-taukin/THERM-MV-taukin$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors13[r])
end

strt_idx = findmin(abs.(bintanja2008_T[1, :] .+ 125e3))[2]
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :], color=:black, linewidth=2.5, linestyle=:dash, label=L"Bintanja and van de Wal (2008) $\,$")
strt_idx = findmin(abs.(barker2011_temp[1, :] .+ 125e3))[2]
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :], color=:grey25, linewidth=2.5, linestyle=:dash, label=L"Barker et al. (2011) $\,$")
fig[0, 1] = Legend(fig, ax, framevisible=false, valign=0, labelsize=25, nbanks=2)

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
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-MV-taukin/THERM-MV-taukin$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["T"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors13[r], linewidth=2)
end

G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=:black, linestyle=:dash, linewidth=2.5)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :], -8e5, 0, 1e3)
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
    dfr = NCDataset("data/runs//exp05/THERM-MV-taukin/THERM-MV-taukin$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors13[r])
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
    dfr = NCDataset("data/runs//exp05/THERM-MV-taukin/THERM-MV-taukin$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors13[r], linewidth=2)
end

# Panel e
ax = Axis(fig[3, 1], ylabel=L"$T_{ice}$ ($^o$C)", xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([-20, -10, 0], convert_strings_to_latex([-20, -10, 0]))
ylims!(ax, (-15, 6))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :t, :b)
text!(ax, -7.95e2, 4, text=L"(e)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-MV-taukin/THERM-MV-taukin$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["Tice"] .- params.degK, color=colors13[r], linewidth=2)
end

# Panel f
ax = Axis(fig[3, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hidexdecorations!(ax)
hideydecorations!(ax)
hidespines!(ax, :t, :b, :l, :r)
text!(ax, 5, 0.23, text=L"(f)$\,$", align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-MV-taukin/THERM-MV-taukin$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["Tice"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors13[r], linewidth=2)
end

# Panel g
ax = Axis(fig[4, 1], ylabel=L"$f_{str}$", xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0.0, 0.5, 1.0], convert_strings_to_latex([0.0, 0.5, 1.0]))
ylims!(ax, (-0.1, 0.8))
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -7.95e2, 0.7, text=L"(g)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-MV-taukin/THERM-MV-taukin$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["fstr"], color=colors13[r], linewidth=2)
end

# Panel h
ax = Axis(fig[4, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.23, text=L"(h)$\,$", align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-MV-taukin/THERM-MV-taukin$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["fstr"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors13[r], linewidth=2)
end

Colorbar(fig[0, 2], width=Relative(2/3), height=Relative(1 / 4), colormap=cgrad(cmap13, length(runs2plot), categorical=true),
    label=L"$\tau_{kin}$ (kyr)",
    limits=(1, length(runs2plot)+1),
    ticks=(taukin_levels[1] .+ 0.5, convert_strings_to_latex(taukin_levels[2])),
    ticklabelsize=28,
    vertical=false, valign=0
)

rowgap!(fig.layout, 0.0)
rowsize!(fig.layout, 0, Relative(0.1))
colsize!(fig.layout, 2, Relative(2 / 6))
resize_to_layout!(fig)

save("figures/fig13.pdf", fig)

# ===========================================================================================================
# Figure AGING-taualpha
# ===========================================================================================================
figname = "exp06_AGING-taualpha_PSD"
runs2plot = range(2, 24)
df1 = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha01.nc")
params = JLD2.load_object("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha01_params.jld2")
selected_AGING_run, selected_AGING_run_idx = "07", 7-1

c1kyr, c400kyr, c1000kyr, c10000kyr = :purple, :darkorange, :firebrick, :darkred
cmap15 = [:slategray2, :lightsteelblue3, :cornflowerblue, :royalblue3, :royalblue4, :darkseagreen1, :darkolivegreen2, :darkolivegreen3, :olive, :darkolivegreen]
colormap15 = cgrad(cmap15, length(runs2plot), categorical=true)
colors = collect(colormap15)

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

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha01.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=c1kyr, linewidth=2, label=L" $\tau_\alpha$ = 1 kyr")


for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(string(runs2plot[r], pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[r], linewidth=1)
    
end

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha25.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=c400kyr, linewidth=2, label=L" $\tau_\alpha$ = 400 kyr")

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha26.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=c1000kyr, linewidth=2, label=L" $\tau_\alpha$ = 1000 kyr")

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha27.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=c10000kyr, linewidth=2, label=L" $\tau_\alpha$ = 10000 kyr")

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(selected_AGING_run).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[selected_AGING_run_idx], linewidth=3)

strt_idx = findmin(abs.(bintanja2008_T[1, :] .+ 125e3))[2]
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :], color=:black, linewidth=2.5, linestyle=:dash, label=L"Bintanja and van de Wal (2008) $\,$")
strt_idx = findmin(abs.(barker2011_temp[1, :] .+ 125e3))[2]
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :], color=:grey25, linewidth=2.5, linestyle=:dash, label=L"Barker et al. (2011) $\,$")
fig[0, 1] = Legend(fig, ax, framevisible=false, valign=0, labelsize=25, nbanks=2)

# Panel b
ax = Axis(fig[1, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=b_label, align=(:left, :center))

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha01.nc")
G, periods = create_PSD(dfr["time"][:], dfr["T"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c1kyr, linewidth=2)


for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(string(runs2plot[r], pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["T"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=1)
end

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha25.nc")
G, periods = create_PSD(dfr["time"][:], dfr["T"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c400kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha26.nc")
G, periods = create_PSD(dfr["time"][:], dfr["T"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c1000kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha25.nc")
G, periods = create_PSD(dfr["time"][:], dfr["T"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c10000kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(selected_AGING_run).nc")
G, periods = create_PSD(dfr["time"][:], dfr["T"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=colors[selected_AGING_run_idx], linewidth=3)

G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=:black, linestyle=:dash, linewidth=2.5)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :], -8e5, 0, 1e3)
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

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha01.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=c1kyr, linewidth=2)


for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(string(runs2plot[r], pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[r], linewidth=1)
    
end

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha25.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=c400kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha26.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=c1000kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha27.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=c10000kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(selected_AGING_run).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[selected_AGING_run_idx], linewidth=3)

# Panel d
ax = Axis(fig[2, 2], xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidexdecorations!(ax)
hidespines!(ax, :b, :t, :l, :r)
text!(ax, 5, 0.23, text=d_label, align=(:left, :center))

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha01.nc")
G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c1kyr, linewidth=2)


for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(string(runs2plot[r], pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=1)
    
end

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha25.nc")
G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c400kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha26.nc")
G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c1000kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha25.nc")
G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c10000kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(selected_AGING_run).nc")
G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=colors[selected_AGING_run_idx], linewidth=3)

# Panel e
ax = Axis(fig[3, 1], ylabel=L"$\alpha$", xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0.2, 0.5, 0.9], convert_strings_to_latex([0.2, 0.5, 0.9]))
ylims!(ax, (0.0, 1.2))
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -7.95e2, 1, text=L"e)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha01.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=c1kyr, linewidth=2)


for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(string(runs2plot[r], pad=ndigits(length(runs2plot)))).nc")
    params = JLD2.load_object("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(string(runs2plot[r], pad=ndigits(length(runs2plot))))_params.jld2")
    lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=colors[r], linewidth=1)
    
end

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha25.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=c400kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha26.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=c1000kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha27.nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=c10000kyr, linewidth=2)


dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(selected_AGING_run).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=colors[selected_AGING_run_idx], linewidth=3)

# Panel f
ax = Axis(fig[3, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.23, text=f_label, align=(:left, :center))

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha01.nc")
G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c1kyr, linewidth=2)


for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(string(runs2plot[r], pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=1)
    
end

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha25.nc")
G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c400kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha26.nc")
G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c1000kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha25.nc")
G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=c10000kyr, linewidth=2)

dfr = NCDataset("data/runs//exp06/AGING-MV-taualpha/AGING-MV-taualpha$(selected_AGING_run).nc")
G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=colors[selected_AGING_run_idx], linewidth=3)

Colorbar(fig[0, 2], width=Relative(2/3), height=Relative(1 / 4), colormap=colormap15,
    label=L"$\tau_\alpha$ ($\mathrm{kyr}$)",
    limits=(1, length(runs2plot)+ 1),
    ticks=([1, 11, 23] .+ 0.5, convert_strings_to_latex(["10", "60", "120"])),
    ticklabelsize=28,
    vertical=false, valign=0, halign=0.5
)

rowgap!(fig.layout, 0.0)
rowsize!(fig.layout, 0, Relative(0.1))
colsize!(fig.layout, 2, Relative(2 / 6))
resize_to_layout!(fig)
save("figures/fig15.pdf", fig)

# ===========================================================================================================
# Figure 10. AGING (best τα)
# ===========================================================================================================
figname = "exp06_AGING_best"
best_exp = "25"
df = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(best_exp).nc")
params = JLD2.load_object("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(best_exp)_params.jld2")

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
G, periods = create_PSD(cut_run(df["time"][:], df["time"][:], -8e5), cut_run(df["time"][:], df["I"][:], -8e5), -8e5, 0, 1e3)
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
G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=bintanja2008color, linewidth=linewidth)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=barker2011color, linewidth=linewidth)
G, periods = create_PSD(cut_run(df["time"][:], df["time"][:], -8e5), cut_run(df["time"][:], df["T"][:], -8e5), -8e5, 0, 1e3)
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
G, periods = create_PSD(luthi2008_co2[1, :], luthi2008_co2[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=luthi2008color, linewidth=linewidth)
G, periods = create_PSD(berends2021_co2[1, :], berends2021_co2[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=berends2021color, linewidth=linewidth)
G, periods = create_PSD(yamamoto2022_co2[1, :], yamamoto2022_co2[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=yamamoto2022color, linewidth=linewidth)
G, periods = create_PSD(cut_run(df["time"][:], df["time"][:], -8e5), cut_run(df["time"][:], df["C"][:], -8e5), -8e5, 0, 1e3)
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
G, periods = create_PSD(cut_run(df["time"][:], df["time"][:], -8e5), cut_run(df["time"][:], df["H"][:], -8e5), -8e5, 0, 1e3)
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
G, periods = create_PSD(bintanja2008_vol[1, :], bintanja2008_vol[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=bintanja2008color, linewidth=linewidth)
G, periods = create_PSD(berends2021_vol[1, :], berends2021_vol[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=berends2021color, linewidth=linewidth)
G, periods = create_PSD(spratt2016_vol[1, :], spratt2016_vol[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=spratt2016color, linewidth=linewidth)
G, periods = create_PSD(cut_run(df["time"][:], df["time"][:], -8e5), cut_run(df["time"][:], df["Vol"][:], -8e5), -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=pacco_color, linewidth=linewidth)

rowgap!(fig.layout, 0.0)
colsize!(fig.layout, 2, Relative(2 / 6))
save("figures/fig16.pdf", fig)

# ===========================================================================================================
# Figure 13. AGING phase space (H)
# ===========================================================================================================
var2plot="H"
figname = "exp06_AGING_PhaseSpace_$(var2plot)"
best_exp = "25"
df = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(best_exp).nc")
params = JLD2.load_object("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(best_exp)_params.jld2")
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

# APPENDIX A

# ===========================================================================================================
# Figure THERM-Cs
# ===========================================================================================================
figname = "exp05_THERM-Cs_PSD"
runs2plot = range(1, length(Css))
colors = collect(cgrad(cmap, length(runs2plot), categorical=true))
xticks = Int.(-8e2:1e2:0)
xlims = (-8e2, 0)
ylims = (-600, 3000)
yticks = [0, 1000, 2000]

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
    dfr = NCDataset("data/runs//exp05/THERM-MV-Cs/THERM-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[r])
end

strt_idx = findmin(abs.(bintanja2008_T[1, :] .+ 125e3))[2]
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :], color=:black, linewidth=2.5, linestyle=:dash, label=L"Bintanja and van de Wal (2008) $\,$")
strt_idx = findmin(abs.(barker2011_temp[1, :] .+ 125e3))[2]
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :], color=:grey25, linewidth=2.5, linestyle=:dash, label=L"Barker et al. (2011) $\,$")
fig[0, 1] = Legend(fig, ax, framevisible=false, valign=0, labelsize=25, nbanks=2)

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
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-MV-Cs/THERM-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["T"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=:black, linestyle=:dash, linewidth=2.5)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :], -8e5, 0, 1e3)
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
    dfr = NCDataset("data/runs//exp05/THERM-MV-Cs/THERM-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
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
    dfr = NCDataset("data/runs//exp05/THERM-MV-Cs/THERM-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

# Panel e
ax = Axis(fig[3, 1], ylabel=L"$T_{ice}$ ($^o$C)", xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([-20, -10, 0], convert_strings_to_latex([-20, -10, 0]))
ylims!(ax, (-15, 6))
xlims!(ax, xlims)
hidexdecorations!(ax)
hidespines!(ax, :t, :b)
text!(ax, -7.95e2, 4, text=L"(e)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-MV-Cs/THERM-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
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
    dfr = NCDataset("data/runs//exp05/THERM-MV-Cs/THERM-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["Tice"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

# Panel g
ax = Axis(fig[4, 1], ylabel=L"$f_{str}$", xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0.0, 0.5, 1.0], convert_strings_to_latex([0.0, 0.5, 1.0]))
ylims!(ax, (-0.1, 0.8))
xlims!(ax, xlims)
hidespines!(ax, :t)
text!(ax, -7.95e2, 0.7, text=L"(g)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-MV-Cs/THERM-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["fstr"], color=colors[r], linewidth=2)
end

# Panel h
ax = Axis(fig[4, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :b, :r)
text!(ax, 5, 0.23, text=h_label, align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp05/THERM-MV-Cs/THERM-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["fstr"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=2)
end

Colorbar(fig[0, 2], width=Relative(2/3), height=Relative(1 / 4), colormap=cgrad(cmap, length(runs2plot), categorical=true),
    label = L"$C_s$ ($\mathrm{m\cdot yr^{-1}\cdot Pa^{-2}}$)",
    limits=(1, length(runs2plot)+1),
    ticks=(Cs_levels[1] .+ 0.5, convert_strings_to_latex(Cs_levels[2])),
    ticklabelsize=28,
    vertical=false, valign=0
)

rowgap!(fig.layout, 0.0)
rowsize!(fig.layout, 0, Relative(0.1))
colsize!(fig.layout, 2, Relative(2 / 6))
resize_to_layout!(fig)

save("figures/figA1.pdf", fig)

# ===========================================================================================================
# Figure AGING-Cs
# ===========================================================================================================

figname = "exp06_AGING-Cs_PSD"
runs2plot = range(1, length(Css))
df1 = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs01.nc")
params = JLD2.load_object("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs01_params.jld2")
selected_AGING_run = 25

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
    dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[r], linewidth=1)
end

dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(selected_AGING_run, pad=ndigits(length(runs2plot)))).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["T"] .- dfr["Tref"], color=colors[selected_AGING_run], linewidth=3)

strt_idx = findmin(abs.(bintanja2008_T[1, :] .+ 125e3))[2]
lines!(ax, bintanja2008_T[1, :] ./ 1e3, bintanja2008_T[2, :], color=:black, linewidth=2.5, linestyle=:dash, label=L"Bintanja and van de Wal (2008) $\,$")
strt_idx = findmin(abs.(barker2011_temp[1, :] .+ 125e3))[2]
lines!(ax, barker2011_temp[1, :] ./ 1e3, barker2011_temp[2, :], color=:grey25, linewidth=2.5, linestyle=:dash, label=L"Barker et al. (2011) $\,$")
fig[0, 1] = Legend(fig, ax, framevisible=false, valign=0, labelsize=25, nbanks=2)

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
for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["T"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=1)
end

dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(selected_AGING_run, pad=ndigits(length(runs2plot)))).nc")
G, periods = create_PSD(dfr["time"][:], dfr["T"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=colors[selected_AGING_run], linewidth=3)

G, periods = create_PSD(bintanja2008_T[1, :], bintanja2008_T[2, :], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=:black, linestyle=:dash, linewidth=2.5)
G, periods = create_PSD(barker2011_temp[1, :], barker2011_temp[2, :], -8e5, 0, 1e3)
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
    dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[r], linewidth=1)
end

dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(selected_AGING_run, pad=ndigits(length(runs2plot)))).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["H"], color=colors[selected_AGING_run], linewidth=3)

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
    dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=1)
end

dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(selected_AGING_run, pad=ndigits(length(runs2plot)))).nc")
G, periods = create_PSD(dfr["time"][:], dfr["H"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=colors[selected_AGING_run], linewidth=3)

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
    dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    params = JLD2.load_object("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(r, pad=ndigits(length(runs2plot))))_params.jld2")
    lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=colors[r], linewidth=1)
end

dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(selected_AGING_run, pad=ndigits(length(runs2plot)))).nc")
lines!(ax, dfr["time"] ./ 1e3, dfr["albedo"], color=colors[selected_AGING_run], linewidth=3)

# Panel f
ax = Axis(fig[3, 2], xlabel=periods_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (ticks_periods, convert_strings_to_latex(ticks_periods))
ylims!(ax, (-0.01, 0.30))
xlims!(ax, ylims_periods)
hideydecorations!(ax)
hidespines!(ax, :t, :l, :r)
text!(ax, 5, 0.23, text=f_label, align=(:left, :center))

for r in eachindex(runs2plot)
    dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(r, pad=ndigits(length(runs2plot)))).nc")
    G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:], -8e5, 0, 1e3)
    lines!(ax, periods ./ 1e3, G, color=colors[r], linewidth=1)
end

dfr = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(string(selected_AGING_run, pad=ndigits(length(runs2plot)))).nc")
G, periods = create_PSD(dfr["time"][:], dfr["albedo"][:], -8e5, 0, 1e3)
lines!(ax, periods ./ 1e3, G, color=colors[selected_AGING_run], linewidth=3)


Colorbar(fig[0, 2], width=Relative(2/3), height=Relative(1 / 4), colormap=cgrad(cmap, length(runs2plot), categorical=true),
    label = L"$C_s$ ($\mathrm{m\cdot yr^{-1}\cdot Pa^{-2}}$)",
    limits=(1, length(runs2plot)+1),
    ticks=(Cs_levels[1] .+ 0.5, convert_strings_to_latex(Cs_levels[2])),
    ticklabelsize=28,
    vertical=false, valign=0
)

rowgap!(fig.layout, 0.0)
rowsize!(fig.layout, 0, Relative(0.1))
colsize!(fig.layout, 2, Relative(2 / 6))
resize_to_layout!(fig)

save("figures/figA2.pdf", fig)

# ===========================================================================================================
# Figure 11. AGING (best) thresholds
# ===========================================================================================================
figname = "exp06_AGING_best_thr"
best_exp = "25"
df = NCDataset("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(best_exp).nc")
params = JLD2.load_object("data/runs//exp06/AGING-MV-Cs/AGING-MV-Cs$(best_exp)_params.jld2")

sw_eff = df["sw"][:] ./ params.kI

xticks = Int.(-8e2:0.5e2:0)
xlims = (-3.5e2, 0)
ylims = (-100, 3500)
yticks = Int.(-1500:1500:3000)

fig = Figure(resolution=(750, 750), fonts=(; regular="TeX"), fontsize=28)

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
lines!(ax, df["time"] ./ 1e3, df["I"] .- params.insol_ref, color=:black, linewidth=3)

# Panel b
ax = Axis(fig[2, 1], ylabel=L"$I_\mathrm{eff}$ ($\mathrm{W \cdot m^{-2}}$)", xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 25, 50], convert_strings_to_latex([0, 25, 50]))
xlims!(ax, xlims)
ylims!(ax, (-5, 60))
hidespines!(ax, :b, :t)
hidexdecorations!(ax)
text!(ax, -3.45e2, 55, text=L"(b)$\,$", align=(:left, :center))
hlines!(ax, 0.0, color=background_color, linestyle=:solid)
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
save("figures/figA3.pdf", fig)