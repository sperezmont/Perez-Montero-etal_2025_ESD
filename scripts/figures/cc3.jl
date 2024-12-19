using Pkg
Pkg.activate(".")
using JLD2, NCDatasets, CairoMakie, Statistics, LaTeXStrings

include(pwd() * "/scripts/misc/figure_definitions.jl")
include(pwd() * "/scripts/misc/spectral_tools.jl")
include(pwd() * "/scripts/misc/tools.jl")

#### Response to M. Verbitsky CC3 (https://doi.org/10.5194/egusphere-2024-1842-CC3) 

# Define some functions
function calc_Hgrnd(H, zl, zb)
    hgrnd = H .- 1030 / 917 .* max.(zl .- zb, 0.0)     # based on Robinson et al. (2020, GMD)

    mask = similar(hgrnd)
    mask[hgrnd.>0] .= 1
    mask[hgrnd.<=0] .= 0

    hgrnd[mask.==0] .= -9999

    return hgrnd, mask
end

function perezmontero_discharge(constant, velocity, icethicknes)
    return constant .* velocity .* icethicknes   # constant = 1/Locn
end

function verbitsky_discharge(constant, velocity, icethicknes; thr=500)
    icethicknes_v = copy(icethicknes)
    icethicknes_v[icethicknes_v.<=thr] .= thr

    return constant .* velocity ./ (icethicknes_v .^ 3)
end

# Load data
path2data = "/home/sergio/entra/proyects/d05_paper-PACCO/Perez-Montero-etal_YYYY_ESD/data/yelmo_runs/cycle_from_yelmo_r20_1.4/yelmo2Dsm.nc"
path2data1D = "/home/sergio/entra/proyects/d05_paper-PACCO/Perez-Montero-etal_YYYY_ESD/data/yelmo_runs/cycle_from_yelmo_r20_1.4/yelmo1D.nc"
path2data_big = "/home/sergio/entra/proyects/d05_paper-PACCO/Perez-Montero-etal_YYYY_ESD/data/yelmo_runs/cycle_from_yelmo_r20_1.4/yelmo2D.nc"

df, df1D, dfbig = NCDataset(path2data), NCDataset(path2data1D), NCDataset(path2data_big)

# Compute mean ice thickness in the grounding lines from PD AIS
pathtoicedata = "/home/sergio/entra/ice_data/"
df_ant = NCDataset(pathtoicedata * "/Antarctica/ANT-32KM/ANT-32KM_TOPO-RTOPO-2.0.1.nc")
xc_ant, yc_ant = df_ant["xc"][:], df_ant["yc"][:]
mask_ant = df_ant["mask"][:, :]
mask_ant[mask_ant.!=2] .= 0
mask_ant[mask_ant.==2] .= 1

H_grnd_line_ant = df_ant["H_ice"][:, :] .- 1030 / 917 .* max.(0.0 .- df_ant["z_bed"][:, :], 0.0)
H_grnd_line_ant[H_grnd_line_ant.<=0] .= -9999   # floating ice

new_H_grnd_ant = copy(H_grnd_line_ant)
for i in 2:length(xc_ant)-1, j in 2:length(yc_ant)-1
    neighbours = [[i, j - 1], [i, j + 1], [i - 1, j], [i + 1, j],
        [i - 1, j - 1], [i - 1, j + 1], [i + 1, j - 1], [i + 1, j + 1]]
    maski = zeros(8)
    for e in eachindex(neighbours)
        if H_grnd_line_ant[neighbours[e][1], neighbours[e][2]] > -9999
            maski[e] = 1
        end
    end

    if count(==(0), maski) > 1
        if sum(maski) > 0
            mask_ant[i, j] = 1
            new_H_grnd_ant[i, j] = H_grnd_line_ant[i, j]
        else
            mask_ant[i, j] = 0
            new_H_grnd_ant[i, j] = -9999
        end
    else
        mask_ant[i, j] = 0
        new_H_grnd_ant[i, j] = -9999
    end
end

mask_ant[1, 1:191] .= 0
mask_ant[191, 1:191] .= 0
mask_ant[1:191, 1] .= 0
mask_ant[1:191, 191] .= 0

replace!(new_H_grnd_ant, -9999 => NaN32)

display(mean(filter(!isnan, new_H_grnd_ant)))

fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, new_H_grnd_ant)
save("kk.png", fig)

# Select relevant variables
time = df["time"][:]
resolution, dt = 32e3, 500 # m, yr

Hice = df["H_ice"][:, :, :]
zsl, zb = df["z_sl"][:, :, :], df["z_bed"][:, :, :]
xc, yc = df["xc"][:], df["yc"][:]
smb, bmb = df["smb"][:, :, :], df["bmb"][:, :, :]
maskbed = df["mask_bed"][:, :, :]
v = df["uxy_s"][:, :, :]

Vsle = df1D["V_sle"][:]
rsle = -1 .* (Vsle .- Vsle[end])

# Plot simulation snapshots
point1, point2, point3, point4 = 101, 161, 215, 225
xticks = Int.(-200e3:20e3:0)
fig = Figure(resolution=(1500, 750), fonts=(; regular="TeX"), fontsize=28)
ax = Axis(fig[1, 1:4], xlabel=time_label, ylabel=L"RSLE (m) $\,$", xaxisposition=:top, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([-120, -100, -80, -60, -40, -20, 0], convert_strings_to_latex([-120, -100, -80, -60, -40, -20, 0]))
xlims!(ax, (-120e3, 0))
lines!(ax, time, rsle, color=:black, linewidth=3)
scatter!(ax, time[point1], rsle[point1], color=:red, markersize=30, marker=:utriangle)
scatter!(ax, time[point2], rsle[point2], color=:red, markersize=30, marker=:diamond)
scatter!(ax, time[point3], rsle[point3], color=:red, markersize=30, marker=:pentagon)
scatter!(ax, time[point4], rsle[point4], color=:red, markersize=30, marker=:dtriangle)

function create_snapshot_panel!(f, pos, point, mark, x, y, vel, h, bed)
    ax = Axis(f[2, pos], xgridvisible=false, ygridvisible=false)
    hidedecorations!(ax)
    xlims!(ax, (x[1], x[end]))
    ylims!(ax, (y[1], y[end]))
    contour!(ax, x, y, bed[:, :, point], color=:grey, levels=[1], linewidth=0.5)
    contourf!(ax, x, y, log10.(vel[:, :, point]), colormap=:balance, levels=log10.([1e-1, 1, 1e1, 1e2, 1e3, 5e3]))
    contour!(ax, x, y, h[:, :, point], color=:black, levels=0:1000:5000)
    scatter!(ax, -4200, 2800, color=:red, marker=mark, markersize=30)
end

create_snapshot_panel!(fig, 1, point1, :utriangle, xc, yc, v, Hice, maskbed)
create_snapshot_panel!(fig, 2, point2, :diamond, xc, yc, v, Hice, maskbed)
create_snapshot_panel!(fig, 3, point3, :pentagon, xc, yc, v, Hice, maskbed)
create_snapshot_panel!(fig, 4, point4, :dtriangle, xc, yc, v, Hice, maskbed)

Colorbar(fig[2, 5], colormap=:balance,
    label=L"$v$ (m$\cdot$yr$^{-1}$)",
    limits=(1, 6),
    ticks=(1:1:6, [L"10^{-1}", L"1", L"10", L"10^{2}", L"10^{3}", L"5\cdot 10^{3}"]),
    ticklabelsize=28,
    vertical=true)

save("figures/cc3_yelmo.png", fig)

# Calculate modeled ice discharge
Hgrnd, mask_grnd = calc_Hgrnd(Hice, zsl, zb)
vgrnd = v

# Define GrIS masks
gris_mask_0 = Hgrnd[:, :, 1]
gris_mask_f = Hgrnd[:, :, end]

gris_mask_0[gris_mask_0.<=0] .= 0
gris_mask_0[gris_mask_0.>0] .= 1

gris_mask_f[gris_mask_f.<=0] .= 0
gris_mask_f[gris_mask_f.>0] .= 1

mask2 = similar(Hgrnd)
for t in eachindex(time)
    mask2[:, :, t] = gris_mask_0
end
gris_mask_0 = mask2

mask2 = similar(Hgrnd)
for t in eachindex(time)
    mask2[:, :, t] = gris_mask_f
end
gris_mask_f = mask2

# Calculate grounded variables
smb_g, bmb_g = smb, bmb
smb_g[mask_grnd.==0] .= -9999
bmb_g[mask_grnd.==0] .= -9999
vgrnd[mask_grnd.==0] .= -9999

# Mask variables with GrIS mask
Hgrnd[gris_mask_0.==1] .= -9999
smb_g[gris_mask_0.==1] .= -9999
bmb_g[gris_mask_0.==1] .= -9999
vgrnd[gris_mask_0.==1] .= -9999

# Replace with NaN
replace!(Hgrnd, -9999 => NaN32)
replace!(smb_g, -9999 => NaN32)
replace!(bmb_g, -9999 => NaN32)
replace!(vgrnd, -9999 => NaN32)

dHdt_g = (Hgrnd[:, :, 2:end] .- Hgrnd[:, :, 1:end-1]) ./ dt
smb_g, bmb_g = smb_g[:, :, 2:end], bmb_g[:, :, 2:end]
grnd_points = sum(mask_grnd[:, :, 2:end], dims=(1, 2))
max_grnd_points = maximum(grnd_points)

qyelmo = dHdt_g .- smb_g .- bmb_g # qyelmo = dHdt_g - smb_g - bmb_g
qyelmo = -1 .* resolution^2 .* sum(replace!(qyelmo, NaN32 => 0.0), dims=(1, 2))

mask_grline = df["mask_bed"][:, :, :]
mask_grline[mask_grline.!=4] .= 0.0
mask_grline[mask_grline.==4] .= 1.0

Hgrline = Hice .* mask_grline

qyelmo2 = resolution .* v .* Hgrline
qyelmo2 = sum(replace!(qyelmo2, NaN32 => 0.0), dims=(1, 2))

# Calculate diagnostice ice discharge
Hgrnd[gris_mask_f.==1] .= -9999
vgrnd[gris_mask_f.==1] .= -9999

replace!(Hgrnd, -9999 => NaN32)
replace!(vgrnd, -9999 => NaN32)

meanH = [mean(filter(!isnan, Hgrnd[:, :, t+1])) for t in eachindex(time[2:end])]
meanv = Vector{Float32}(undef, length(meanH))
for t in eachindex(time[2:end])
    numpoints = sum(mask_grnd[:, :, t+1], dims=(1, 2))[1, 1, 1]
    if numpoints > 0.1 * max_grnd_points
        meanv[t] = mean(filter(!isnan, vgrnd[:, :, t+1]))
    else
        meanv[t] = 0.0
    end
end

qpacco = perezmontero_discharge(2e13 * 1 / 1e6, meanv, meanH)
qmv_selected = verbitsky_discharge(1.5e13 * 1e6 / 0.9^2, meanv, meanH, thr=700)

qmv, constants_verb = [], [5e11, 2.5e12, 5e12, 7.5e12, 2.5e13, 5e13, 7.5e13, 2.5e14, 5e14] .* 1e6 / 0.9^2
for c in constants_verb
    push!(qmv, verbitsky_discharge(c, meanv, meanH, thr=500))
end

qmv2, thrs = [], [1, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600]
for t in thrs
    push!(qmv2, verbitsky_discharge(2.5e12 * 1e6 / 0.9^2, meanv, meanH, thr=t))
end

# Compute and compare mean driving stress
taud_big, H_big, mask_big, time_big = dfbig["taud"][:, :, :], dfbig["H_ice"][:, :, :], dfbig["mask_bed"][:, :, :], dfbig["time"][:]

mask_big[mask_big .<= 1.0] .= 0.0
mask_big[mask_big .>= 4.0] .= 0.0
mask_big[mask_big .!= 0.0] .= 1.0
numpoints_big = [sum(mask_big[:, :, t], dims=(1, 2))[1, 1, 1] for t in eachindex(time_big)]

taud_big[mask_big .== 0.0] .= NaN32 
H_big[mask_big .== 0.0] .= NaN32 
maskgrisbig = similar(mask_big)
for t in eachindex(time_big)
    maskgrisbig[:, :, t] = mask_big[:, :, 1]
end
gris_mask_big = maskgrisbig

taud_big[gris_mask_big .== 1.0] .= NaN32 
H_big[gris_mask_big .== 1.0] .= NaN32 

mean_taud = [mean(filter(!isnan, taud_big[:, :, t])) for t in eachindex(time_big)]
mean_Hbig = [mean(filter(!isnan, H_big[:, :, t])) for t in eachindex(time_big)]

z = df["z_srf"][:, :, :]
z[isnan.(Hgrnd)] .= NaN32
meanz = [mean(filter(!isnan, z[:, :, t+1])) for t in eachindex(time[2:end])]

taud_pacco = 910 .* 9.8 .* meanH .* meanz ./ 1e6 
taud_verbitsky = 910 .* 9.8 .* meanH .* meanz ./ (0.9 .* meanH .^ 2)

t_selected = 11
fig = Figure(resolution=(1500, 500))
ax = Axis(fig[1, 1], title="mask, 20 kyr BP")
heatmap!(ax, mask_big[:, :, t_selected])
ax = Axis(fig[1, 2], title="taud, 20 kyr BP")
hm = heatmap!(ax, taud_big[:, :, t_selected], colorrange=(0, 1e5), lowclip=:white, highclip=:red)
ax = Axis(fig[1, 3], ylabel="mean taud", xgridvisible=false, ygridvisible=false)
xlims!(ax, (-120e3, 0))
hidexdecorations!(ax)
lines!(ax, time_big, mean_taud, color=:black, label="Yelmo")
lines!(ax, time[2:end], taud_pacco, color=:royalblue, label="PACCO")
lines!(ax, time[2:end], taud_verbitsky, color=:darkorange, label="Verbitsky")
axislegend(ax, framevisible=false, nbanks=1, position=(:left, :center), patchsize=(40, 20))
ax = Axis(fig[1, 3], ylabel="mean H", yaxisposition=:right, xgridvisible=false, ygridvisible=false, ylabelcolor=:grey, yticklabelcolor=:grey)
xlims!(ax, (-120e3, 0))
lines!(ax, time_big, mean_Hbig, color=:grey, linestyle=:dash)
lines!(ax, time[2:end], meanH, color=:purple, linestyle=:dash)
scatter!(ax, time_big[t_selected], mean_Hbig[t_selected], color=:grey, markersize=10)
Colorbar(fig[2, 2], hm, vertical=false)
save("kk2.png", fig)

# Plot ice discharge
xticks = Int.(-200e3:20e3:0)
cmap = [:coral1, :darkorange, :darkred]
colormap = cgrad(cmap, length(constants_verb), categorical=:true)
colormap2 = cgrad(cmap, length(thrs), categorical=:true)
colors = collect(colormap)
colors2 = collect(colormap2)

fig = Figure(resolution=(1000, 500), fonts=(; regular="TeX"), fontsize=28)
ax = Axis(fig[1, 1], ylabel=L"$q$ (10$^{12}$ m$^3 \cdot $yr$^{-1}$)", xlabel=time_label, xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 1e12, 2e12, 3e12, 4e12], convert_strings_to_latex([0, 1, 2, 3, 4]))
ylims!(ax, (-1e10, 5e12))
xlims!(ax, (-120e3, 0))
lines!(ax, time, qyelmo2[1, 1, :], color=:black, linewidth=3, label=L"$q_\mathrm{Yelmo}$")
lines!(ax, time[2:end], qpacco, color=:royalblue, linewidth=3, label=L"$q=k_1 \cdot \bar{v} \cdot \bar{H}$")
lines!(ax, time[2:end], qmv_selected, color=:darkorange, linewidth=3, label=L"$q=k_2 \cdot \bar{v} \cdot \bar{H} / L^2$")
axislegend(ax, framevisible=false, nbanks=1, position=(:left, :top), patchsize=(40, 20))

save("figures/cc3_diagnostic.png", fig)

# Plot ice discharge
xticks = Int.(-200e3:20e3:0)
cmap = [:coral1, :darkorange, :darkred]
colormap = cgrad(cmap, length(constants_verb), categorical=:true)
colormap2 = cgrad(cmap, length(thrs), categorical=:true)
colors = collect(colormap)
colors2 = collect(colormap2)

fig = Figure(resolution=(1500, 750), fonts=(; regular="TeX"), fontsize=28)
ax = Axis(fig[1, 1], ylabel=L"$q$ (10$^{12}$ m$^3 \cdot $yr$^{-1}$)", xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 1e12, 2e12, 3e12, 4e12], convert_strings_to_latex([0, 1, 2, 3, 4]))
hidexdecorations!(ax)
ylims!(ax, (-1e10, 5e12))
xlims!(ax, (-120e3, 0))
lines!(ax, time, qyelmo2[1, 1, :], color=:black, linewidth=3, label=L"$q_\mathrm{Yelmo}$")
lines!(ax, time[2:end], qpacco, color=:royalblue, linewidth=3, label=L"$k_1 \cdot v \cdot H$")
for i in eachindex(constants_verb)
    lines!(ax, time[2:end], qmv[i], color=colors[i], linestyle=:dash, linewidth=3)
end

Colorbar(fig[1, 2], colormap=colormap,
    label=L"$k_2$ (10$^{20}$ m$^5$)",
    limits=(0, length(constants_verb)),
    ticks=(1:1:length(constants_verb), convert_strings_to_latex(round.(constants_verb ./ 1e20, digits=3))),
    # ticks=([1, length(constants_verb)] .- 0.5, convert_strings_to_latex([round(constants_verb[1] / 1e20,
    #         digits=3), round(constants_verb[end] / 1e20, digits=3)])),
    ticklabelsize=28,
    vertical=true)

ax = Axis(fig[2, 1], xlabel=L"Time (yr BP)$\,$", ylabel=L"$q$ (10$^{12}$ m$^3 \cdot $yr$^{-1}$)", xgridvisible=false, ygridvisible=false)
ax.xticks = (xticks, convert_strings_to_latex(-1 .* xticks))
ax.yticks = ([0, 1e12, 2e12, 3e12, 4e12], convert_strings_to_latex([0, 1, 2, 3, 4]))
ylims!(ax, (-1e10, 5e12))
xlims!(ax, (-120e3, 0))
lines!(ax, time, qyelmo2[1, 1, :], color=:black, linewidth=3, label=L"$q_\mathrm{Yelmo}$")
lines!(ax, time[2:end], qpacco, color=:royalblue, linewidth=3, label=L"$k_1 \cdot v \cdot H$")
for i in eachindex(thrs)
    lines!(ax, time[2:end], qmv2[i], color=colors2[i], linestyle=:dash, linewidth=3)
end
lines!(ax, time[2:end], qmv2[1] .* NaN, color=:grey20, linestyle=:dash, linewidth=3, label=L"k_2 \cdot v / H^3")

axislegend(ax, framevisible=false, nbanks=3, position=(:center, :top), patchsize=(40, 20))

Colorbar(fig[2, 2], colormap=colormap2,
    label=L"Min. threshold (m)$\,$",
    limits=(0, length(thrs)),
    ticks=(1:1:length(thrs), convert_strings_to_latex(round.(thrs, digits=0))),
    # ticks=([1, length(thrs)] .- 0.5, convert_strings_to_latex([round(thrs[1],
    #         digits=0), round(thrs[end], digits=0)])),
    ticklabelsize=28,
    vertical=true)

colsize!(fig.layout, 2, Relative(0.01))

save("figures/cc3_supplementary.png", fig)