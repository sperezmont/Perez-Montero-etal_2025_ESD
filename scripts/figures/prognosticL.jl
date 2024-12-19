using Pkg
Pkg.activate(".")
using JLD2, NCDatasets, CairoMakie, Statistics, LaTeXStrings

include(pwd() * "/scripts/misc/figure_definitions.jl")
include(pwd() * "/scripts/misc/spectral_tools.jl")
include(pwd() * "/scripts/misc/tools.jl")

#### AGING configuration using L = coeff * H^2 in response to M. Verbitsky CC (https://doi.org/10.5194/egusphere-2024-1842-CC1) 
# Loading data
main_path = "/home/sergio/entra/models/pacco_vers/pacco_v0.6/output/Perez-Montero-etal_2024_ESD/exp06/"
new_path = "$(main_path)/AGING-MV-c/AGING-MV-c"
df_old = NCDataset("$(main_path)/AGING-Cs/AGING-Cs1/pacco.nc")    # Cs = ??

calc_L(S) = sqrt(S)*1e3 # m
calc_c(S, H) = calc_L(S)/H^2  # m⁻¹

runs2plot = 8:1:8

min_c = calc_c(50e3, 1000)
max_c = calc_c(14e6, 2000)  # based on Antarctica
hrz_coeffs = range(round(min_c, digits=1), round(max_c, digits=1), length=8)
hrz_coeffs = hrz_coeffs[runs2plot]
xticks = Int.(-8e2:1e2:0)

cmap = [:purple, :royalblue, :olive, :darkorange, :darkred]
cmap2 = [:purple, :royalblue, :olive, :darkorange, :darkred, :black]
colormap = cgrad(cmap, 8, categorical=:true)
colors = collect(colormap)

fig = Figure(resolution=(1000, 750), fonts=(; regular="TeX"), fontsize=28)
ax = Axis(fig[1, 1], ylabel=L"$H$ (m)")
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 1000, 2000], convert_strings_to_latex([0, 1000, 2000]))
hidexdecorations!(ax)
lines!(ax, df_old["time"][:] ./ 1e3, df_old["H"][:], color=:black, linewidth=3)
for i in runs2plot
    df = NCDataset("$(new_path)$(i)/pacco.nc")
    lines!(ax, df["time"][:] ./ 1e3, df["H"][:], color=colors[i], linewidth=3)
end

ax = Axis(fig[2, 1], ylabel=L"$L$ (km)")
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([1000, 2000, 3000, 4000], convert_strings_to_latex([1000, 2000, 3000, 4000]))
hidexdecorations!(ax)
lines!(ax, df_old["time"][:] ./ 1e3, df_old["L"][:] ./ 1e3, color=:black, linewidth=3)
for i in runs2plot
    df = NCDataset("$(new_path)$(i)/pacco.nc")
    lines!(ax, df["time"][:] ./ 1e3, df["L"][:] ./ 1e3, color=colors[i])
end

ax = Axis(fig[3, 1], ylabel=L"$\tau_d$ (10$^4$ Pa)")
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 2, 4, 6], convert_strings_to_latex([0, 2, 4, 6]))
hidexdecorations!(ax)
lines!(ax, df_old["time"][:] ./ 1e3, df_old["taud"][:] ./ 1e4, color=:black, linewidth=3)
for i in runs2plot
    df = NCDataset("$(new_path)$(i)/pacco.nc")
    lines!(ax, df["time"][:] ./ 1e3, df["taud"][:] ./ 1e4, color=colors[i])
end

ax = Axis(fig[4, 1], ylabel=L"$q$ (m$\cdot$yr$^{-1}$)", xlabel=L"Time (kyr BP)$\,$", yscale=Makie.pseudolog10)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 0.004, 0.008], convert_strings_to_latex([0, 0.004, 0.008]))
lines!(ax, df_old["time"][:] ./ 1e3, df_old["v"][:] .* df_old["H"][:] ./ 1e6, color=:black, linewidth=3)
for i in runs2plot
    df = NCDataset("$(new_path)$(i)/pacco.nc")
    lines!(ax, df["time"][:] ./ 1e3, df["v"][:] .* df["H"][:] ./ 1e6, color=colors[i])
end

# ax = Axis(fig[0, 1])
# hidedecorations!(ax)
# hidespines!(ax)
# lines!(ax, df_old["time"][:] ./ 1e3, df_old["H"][:] .* NaN, color=:black, linewidth=3, label=L"fixed $L$")
# axislegend(ax, framevisible=false, position=(:right, :center), labelsize=28, nbanks=1)

# rowsize!(fig.layout, 0, Relative(0.1))

save("figures/prognosticL.png", fig)

