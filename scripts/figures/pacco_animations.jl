using Pkg
Pkg.activate(".")
using NCDatasets, JLD2, Interpolations, DSP, LaTeXStrings, CairoMakie

include(pwd() *"/scripts/misc/figure_definitions.jl")
include(pwd() * "/scripts/misc/spectral_tools.jl")
include(pwd() * "/scripts/misc/tools.jl")

function calc_radius(data::NCDataset)
    Surfaceoc = 3.618e8          # [km²] Oceanic surface
    rhoi = 910.0                 # [kg/m³] Ice density 
    rhow = 1000.0                # [kg/m³] Water density
    
    Vol, H = data["Vol"][:], data["H"][:]   # eliminate extremely low values (bad values)
    Vol[abs.(Vol) .< 1] .= 0.0
    Vol[abs.(H) .< 1] .= 0.0

    Surface = -1 .* Vol ./ H .* (rhow * Surfaceoc / rhoi)   # Total surface
    Surface[Surface .<= 0] .= 0.0   # avoid small but negative values

    Radius = sqrt.(Surface ./ pi)
    Radius[isnan.(Radius)] .= 0.0

    return Radius   # km
end

function calc_ice_profile(Z, R, bR; n=10)
    # Assume profile dependent on Z, b(R) and R:   z(r) = A⋅r³ + Z
    # since z(R) = b(R)
    A = (bR .- Z) ./ (R^3) 

    if Z <= bR
        A = -1 .* (bR .- Z) ./ R^3
    elseif R < 0.1
        A = 0    
    end

    # r-axis
    r = range(0, R, length=n)

    return r, A .* r .^3 .+ Z    
end

function calc_hice_profile(H, R;n=10)
    # follows Vialov profile
    # h(r) = H[1 - (r/R)^a]^b
    # a = (m+1)/m
    # b = m/(2m + 2)
    m = 1
    a = (m+1) / m
    b = m / (2*m + 2)

    # r-axis
    r = range(0, R, length=n)

    return r, H .* (1 .- (r ./ R) .^ a) .^ b 
end

function calc_bed_profile(B, R, yR; n=10)
    # Assume a profile dependent on B and R:   b(r) = A⋅r³ + B
    # since b(R) = yR (local approximation), yR >= B always!!
    A = (yR .- B)  ./ (R^3)

    if R < 0.1
        A, B = 0, yR    
    end

    # r-axis
    r = range(0, R, length=n)
    return r, A .* r .^3 .+ B
end

function update_profiles(H, Z, B, R, SED, Brelax, xaxislimit, mode)
    extra_domain = R:xaxislimit # outside of the ice

    # 1. Bedrock profile
    xb, yb = calc_bed_profile(B, R, Brelax)
    new_xb = vcat(-1 .* reverse(extra_domain), -1 .* reverse(xb), xb, extra_domain)
    new_yb = vcat(Brelax .* ones(length(extra_domain)), reverse(yb), yb, Brelax .* ones(length(extra_domain))) 

    # 2. Ice surface profile
    if mode == "vialov"
        xz, h = calc_hice_profile(H, R)
        yz = h .+ SED .+ yb
    else
        xz, yz = calc_ice_profile(Z, R, yb[end])   
    end
    new_xz = vcat(-1 .* reverse(extra_domain), -1 .* reverse(xz), xz, extra_domain)
    new_yz = vcat(Brelax .* ones(length(extra_domain)), reverse(yz), yz, Brelax .* ones(length(extra_domain))) 

    # 3. Sediments profile
    xs = new_xb
    ys = new_yb .+ SED

    return new_xz, new_yz, new_xb, new_yb, xs, ys
end

function create_animation(experiment::String, frames::Any;where2save::String="./", fmt=".gif", framerate::Real=30, fps::Real=60, limsz::Tuple=(-600, 2000), limsR::Tuple=(-3500, 3500), path2proxy::String=pwd()*"/data/paleo_records/", plot_sun=false, plot_albedo=true)
    sun_sufix = "_s$(Int(plot_sun))"
    albedo_sufix = "_a$(Int(plot_albedo))"

    plotname = where2save * "/pacco_animation$(sun_sufix)$(albedo_sufix)$(fmt)"    

    df = NCDataset(experiment)
    
    R = calc_radius(df)
    z = df["z"][:]
    z[R .<= 0] .= 0.0 

    limst = (df["time"][frames[1]] ./ 1e3, df["time"][frames[end]] ./ 1e3)
    
    ice_color = collect(cgrad([:grey10, :dimgrey, :lightslategrey, :slategray3, :lightsteelblue1], 10, categorical=true))
    ice_color_ranges = range(0.2, 0.9, length=10)

    ice_vol_proxy1 = JLD2.load_object("$(path2proxy)/bintanja-vandewal_2008_3000.jld2")
    ice_vol_proxy2 = JLD2.load_object("$(path2proxy)/spratt-lisiecki_2016.jld2")

    # . Define figure
    fig = Figure(resolution=(1300, 1000), fonts=(; regular="TeX"), fontsize=28)

    # . Define axes
    ax1, ax2, ax3 = Axis(fig[1, 1], ylabel=L"$I$ (W$\cdot$m$^{-2}$)", xlabel=time_label, xgridvisible=false, ygridvisible=false, xaxisposition=:top),
                    Axis(fig[2, 1], ylabel=L"$V_\mathrm{ice}$ (m)", xgridvisible=false, ygridvisible=false, yaxisposition=:left),
                    Axis(fig[3, 1], ylabel=z_label, xlabel=L"Distance from the center (km) $\,$", xgridvisible=false, ygridvisible=false, yaxisposition=:left)   # Insolation, ice volume, ice sheet evolution
    
    ax1.xticks = (Int.(-3e3:100:100), convert_strings_to_latex(-1 .* Int.(-3e3:100:100)))
    ax1.yticks = ([450, 500, 550], convert_strings_to_latex([450, 500, 550]))

    ax2.xticks = (Int.(-3e3:100:100), convert_strings_to_latex(-1 .* Int.(-3e3:100:100)))
    ax2.yticks = ([-120, -60, 0], convert_strings_to_latex([-120, -60, 0]))

    ax3.xticks = (Int.(-3e3:1e3:3e3), convert_strings_to_latex(Int.(-3e3:1e3:3e3)))
    ax3.yticks = ([0, 500, 1000, 1500, 2000], convert_strings_to_latex([0, 500, 1000, 1500, 2000]))
    
    xlims!(ax1, limst)
    xlims!(ax2, limst)
    ylims!(ax2, (-150, 80))
    ylims!(ax3, limsz)
    xlims!(ax3, limsR)

    hidexdecorations!(ax2)

    hidespines!(ax1, :b)
    hidespines!(ax2, :t)
    hidespines!(ax3, :t, :r)

    rowgap!(fig.layout, 1, 0.0)

    if plot_albedo
        Colorbar(fig[3, 1], width=Relative(1/6), height=Relative(1/15), colormap=ice_color,
                        label=L"Albedo $\,$",
                        limits=(0, 1.0),
                        ticks=([0, 0.5, 1.0], convert_strings_to_latex([0.2, 0.5, 0.9])),
                        ticklabelsize=28, halign=0.95, valign=0.85,
                        vertical=false
                    )
        rowgap!(fig.layout, 2, 10.0)
        legendpos = (0.5, 1)
    else
        rowgap!(fig.layout, 2, 50.0)
        legendpos = (0.5, 1)
    end

    rowsize!(fig.layout, 3, Relative(4 / 7))

    # . Plot static curves
    lines!(ax1, df["time"] ./ 1e3, df["I"], color=:grey20)

    lines!(ax2, ice_vol_proxy1[1, :] ./ 1e3, ice_vol_proxy1[2, :], color=:royalblue, label=L"Bintanja and van de Wal, 2008$\,$")
    lines!(ax2, ice_vol_proxy2[1, :] ./ 1e3, ice_vol_proxy2[2, :], color=:darkorange, label=L"Spratt and Lisiecki, 2016$\,$")
    lines!(ax2, df["time"] ./ 1e3, df["Vol"] .- 30000, color=:black, label=L"PACCO $\,$")
    axislegend(ax2, framevisible=false, position=legendpos, labelsize=25, nbanks=3)

    # . Create animation
    p1 = scatter!(ax1, df["time"][1] ./ 1e3, df["I"][1], color=:grey20, markersize=1)
    s3_sun = scatter!(ax3, -100000, 1500, markersize=df["I"][1], color=:yellow)
    p2 = scatter!(ax2, df["time"][1] ./ 1e3, df["T"][1] .- df["Tref"][1], color=:grey20, markersize=1)
    b3_z = band!(ax3, -3000:3000, -10000, -8000, color=:grey)   
    b3_bed = band!(ax3, -3000:3000, -10000, -8000, color=:grey)   
    b3_sed = band!(ax3,-3000:3000, -10000, -8000, color=:sienna4)

    ice = false
    record(fig, plotname, frames, framerate=framerate) do frame
        # cleaning
        if frame > 1
            delete!(ax1, p1)
            if plot_sun
                delete!(ax3, s3_sun)
            end
            delete!(ax2, p2)
            if ice
                delete!(ax3, b3_z)
            end
            delete!(ax3, b3_bed)
            delete!(ax3, b3_sed)
        end
        
        # Insolation
        lines!(ax1, df["time"][1:frame] ./ 1e3, df["I"][1:frame], color=:black, linewidth=2)
        p1 = scatter!(ax1, df["time"][frame] ./ 1e3, df["I"][frame], color=:black, strokecolor=:grey, strokewidth=2, markersize=15)

        if plot_sun
            s3_sun = scatter!(ax3, -3000, 1500, markersize=df["I"][frame]- 420, color=:yellow, strokecolor=:grey, strokewidth=1)
        end

        # Ice volume
        lines!(ax2, df["time"][1:frame] ./ 1e3, df["Vol"][1:frame], color=:black, linewidth=2)
        p2 = scatter!(ax2, df["time"][frame] ./ 1e3, df["Vol"][frame], color=:black, strokecolor=:grey, strokewidth=2, markersize=15)

        # Ice sheet
        xz, yz, xb, yb, xs, ys = update_profiles(df["H"][frame], z[frame], df["B"][frame], R[frame], df["Hsed"][frame], 500, limsR[2], "")

        if df["H"][frame] > 0
            if plot_albedo
                color_idx = findmin(abs.(ice_color_ranges .- df["albedo"][frame]))[2]   # read color
                color = ice_color[color_idx]
            else
                color=:snow3
            end
            b3_z = band!(ax3, xz, yb, yz, color=color)   

            ice = true
        else
            ice = false
        end

        b3_bed = band!(ax3, xb, -1000, yb, color=:peachpuff4)
        b3_sed = band!(ax3, xb, yb, ys, color=:sienna4)

        sleep(1/fps) # refreshes the display!
    end
end

create_animation("data/runs/paper-800kyr/exp06/AGING-SLOW/AGING-SLOW1.nc", where2save=pwd()*"/figures/", 100:1:900, framerate=15, fps=60, plot_sun=true, plot_albedo=true)






