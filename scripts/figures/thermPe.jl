using Pkg
Pkg.activate(".")
using JLD2, NCDatasets, CairoMakie, Statistics, LaTeXStrings

include(pwd() * "/scripts/misc/figure_definitions.jl")
include(pwd() * "/scripts/misc/spectral_tools.jl")
include(pwd() * "/scripts/misc/tools.jl")

#### THERM configuration before and after M. Verbitsky CC2 (https://doi.org/10.5194/egusphere-2024-1842-CC2) 
# Loading data
main_path = "/home/sergio/entra/models/pacco_vers/pacco_v0.6/output/Perez-Montero-etal_2024_ESD/exp05/"
df_new = NCDataset("$(main_path)/THERM-MV-taukin/THERM-MV-taukin5/pacco.nc")
df_old = NCDataset("/home/sergio/entra/models/pacco_vers/pacco_preMV/Pacco.jl/output/paper-800kyr/exp05/THERM-taukin/THERM-taukin5/pacco.nc")    

xticks = Int.(-8e2:1e2:0)

fig = Figure(resolution=(1000, 750), fonts=(; regular="TeX"), fontsize=28)
ax = Axis(fig[1, 1], ylabel=L"$H$ (m)")
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 1000, 2000], convert_strings_to_latex([0, 1000, 2000]))
xlims!(-9e2, 0)
hidexdecorations!(ax)
lines!(ax, df_old["time"][:] ./ 1e3, df_old["H"][:], color=:royalblue, linewidth=3)
lines!(ax, df_new["time"][:] ./ 1e3, df_new["H"][:], color=:darkorange, linewidth=3)

ax = Axis(fig[2, 1], ylabel=L"$T$ ($^\mathrm{o}$C)")
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([-20, -10, 0], convert_strings_to_latex([-20, -10, 0]))
xlims!(-9e2, 0)
hidexdecorations!(ax)
lines!(ax, df_old["time"][:] ./ 1e3, df_old["T"][:] .- df_old["Tref"][:], color=:royalblue, linewidth=3)
lines!(ax, df_new["time"][:] ./ 1e3, df_new["T"][:] .- df_new["Tref"][:], color=:darkorange, linewidth=3)

ax = Axis(fig[3, 1], ylabel=L"$Pe$")
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 5, 10], convert_strings_to_latex([0, 5, 10]))
xlims!(-9e2, 0)
hidexdecorations!(ax)
lines!(ax, df_new["time"][:] ./ 1e3, df_new["Pe"][:], color=:darkorange, linewidth=3)

ax = Axis(fig[4, 1], ylabel=L"$T_\mathrm{ice}$ (ºC)")
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([-20, -10, 0], convert_strings_to_latex([-20, -10, 0]))
xlims!(-9e2, 0)
hidexdecorations!(ax)
lines!(ax, df_old["time"][:] ./ 1e3, df_old["Tice"][:] .- df_old["Tref"][:], color=:royalblue, linewidth=3)
lines!(ax, df_new["time"][:] ./ 1e3, df_new["Tice"][:] .- df_new["Tref"][:], color=:darkorange, linewidth=3)

ax = Axis(fig[5, 1], ylabel=L"$f_\mathrm{str}$", xlabel=time_label)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 0.25, 0.5], convert_strings_to_latex([0, 0.25, 0.5]))
xlims!(-9e2, 0)
lines!(ax, df_old["time"][:] ./ 1e3, df_old["fstr"][:], color=:royalblue, linewidth=3)
lines!(ax, df_new["time"][:] ./ 1e3, df_new["fstr"][:], color=:darkorange, linewidth=3)

save("figures/thermPe.png", fig)
