using Pkg
Pkg.activate(".")
using JLD2, Interpolations, CairoMakie, LaTeXStrings, Colors

include(pwd() *"/scripts/misc/figure_definitions.jl")
include(pwd() * "/scripts/misc/spectral_tools.jl")
include(pwd() * "/scripts/misc/tools.jl")

function create_data(file::String, color::Symbol, year2start::Real, dt::Real)
    data = JLD2.load_object(file)
    t, x = data[1, :], data[2, :]

    if t[1] < year2start
        t, x = interp_time_series(t, x, year2start, 0.0, dt)
        data = [t x]'
    end

    x_trend = sum(x) ./ length(x) .+ (x[end] .- x[1]) ./ (t[end] .- t[1]) .* t

    colormap = cgrad([:white, color], 8)
    Wnorm, periods = create_WV(t, x .- x_trend)

    return [data, Wnorm, periods, colormap]
end



# Plotting
fig = Figure(resolution=(2050, 800), fonts=(; regular="TeX"), fontsize=36)

xticks_time = xticks_time_800
xticks_time = (xticks_time, convert_strings_to_latex(-1 .* xticks_time))
xlims_time = xlims_time_800

ticks_periods = (ticks_periods, convert_strings_to_latex(ticks_periods))

# Panel a
file, color, text_a, text_b = "data/insolation/solstice_insolation_65N170_10yr_5MyrBP-0.jld2", laskar2004color, (550, a_label), (100, b_label)
data = create_data(file, color, xlims_time_800[1] * 1e3, 1e3)

ax = Axis(fig[1, 1], ylabel=ins_label, xgridvisible=false, ygridvisible=false)
ax.xticks = xticks_time
ax.yticks = (yticks_insol, convert_strings_to_latex(yticks_insol))
ylims!(ax, ylims_ins)
xlims!(ax, xlims_time)
hidespines!(ax, :b)
hidexdecorations!(ax)
lines!(ax, data[1][1, :] ./ 1e3, data[1][2, :], color=color, linewidth=5)
text!(ax, -7.95e2, text_a[1], text=text_a[2])

# Panel b
file1, color1 = "data/paleo_records/lisiecki-raymo_2005_d18O_800.jld2", lisiecki2005color
data1 = create_data(file1, color1, xlims_time_800[1] * 1e3, 1e3)
file2, color2 = "data/paleo_records/hodell-etal_2023_d18O_800.jld2", hodell2023color
data2 = create_data(file2, color2, xlims_time_800[1] * 1e3, 1e3)

ax = Axis(fig[2, 1], xlabel=time_label, ylabel=d18O_label, xgridvisible=false, ygridvisible=false, yreversed=true)
ax.yticks = (yticks_d18O, convert_strings_to_latex(yticks_d18O))
ax.xticks = xticks_time
xlims!(ax, xlims_time)
hidespines!(ax, :t)
lines!(ax, data1[1][1, :] ./ 1e3, data1[1][2, :], color=color1, linewidth=5)
lines!(ax, data2[1][1, :] ./ 1e3, data2[1][2, :], color=color2, linewidth=5)

text!(ax, -7.95e2, 3.4, text=b_label)
# Panel f
files = ["data/paleo_records/" .* ["barker-etal_2011.jld2", "hodell-etal_2023_d18O_800.jld2",
        "lisiecki-raymo_2005_d18O_800.jld2",
        "luthi-etal_2008.jld2",
        "yamamoto-etal_C_800.jld2",
        "spratt-lisiecki_2016.jld2",
        "bintanja-vandewal_2008_800.jld2"]
    "data/insolation/" .* ["solstice_insolation_65N170_10yr_5MyrBP-0.jld2"]]
data = [JLD2.load_object(fi) for fi in files]

Gvector, Pvector = [], []
for d in eachindex(data)
    t, x = data[d][1, :], data[d][2, :]
    x_trend = sum(x) ./ length(x) .+ (x[end] .- x[1]) ./ (t[end] .- t[1]) .* t
    G, periods = create_PSD(t, x .- x_trend)

    push!(Gvector, G)
    push!(Pvector, periods)
end

ax = Axis(fig[2, 2], xlabel=periods_label, xgridvisible=false)
hideydecorations!(ax)
hidespines!(ax, :l, :r, :t)
ax.xticks = ticks_periods
xlims!(ax, (-10, 150))

lines!(ax, Pvector[end] ./ 1e3, Gvector[end], color=laskar2004color, label=L"Boreal SSI$\,$", linewidth=5)
lines!(ax, Pvector[3] ./ 1e3, Gvector[3], color=lisiecki2005color, label=L"Lisiecki and Raymo (2005)$\,$", linewidth=5)
# lines!(ax, Pvector[7] ./ 1e3, Gvector[7], color=bintanja2008color, label=L"BW2008$\,$", linewidth=5)
# lines!(ax, Pvector[4] ./ 1e3, Gvector[4], color=luthi2008color, label=L"Lü2008$\,$", linewidth=5)
# lines!(ax, Pvector[5] ./ 1e3, Gvector[5], color=yamamoto2022color, label=L"Ya2022$\,$", linewidth=5)
# lines!(ax, Pvector[1] ./ 1e3, Gvector[1], color=barker2011color, label=L"Ba2011$\,$", linewidth=5)
# lines!(ax, Pvector[6] ./ 1e3, Gvector[6], color=spratt2016color, label=L"SL2016$\,$", linewidth=5)
lines!(ax, Pvector[2] ./ 1e3, Gvector[2], color=hodell2023color, label=L"Hodell et al. (2023)$\,$", linewidth=5)
#text!(ax, 100, 0.3, text=L"The 100 kyr problem$\,$", align=(:center, :center), fontsize=25)
text!(ax, 10, 0.22, text=c_label)
fig[1, 2] = Legend(fig, ax, framevisible=false, titlesize=20, labelsize=32, nbanks=1, patchsize=(40, 20))

rowgap!(fig.layout, 0.1)
colsize!(fig.layout, 1, Relative(4 / 6))
colsize!(fig.layout, 2, Relative(2 / 6))
save("figures/fig01.pdf", fig)