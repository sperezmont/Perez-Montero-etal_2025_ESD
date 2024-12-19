using Pkg
Pkg.activate(".")
using JLD2, NCDatasets, CairoMakie, Statistics, LaTeXStrings

function cut_run(t, x, time2cut1, time2cut2)
    minidx = findmin(abs.(t .- time2cut1))[2]
    maxidx = findmin(abs.(t .- time2cut2))[2]
    return x[minidx:maxidx]
end

include(pwd() * "/scripts/misc/figure_definitions.jl")
include(pwd() * "/scripts/misc/spectral_tools.jl")
include(pwd() * "/scripts/misc/tools.jl")

taus = vcat(1:1:9, 10:10:90, 1e2:1e2:9e2, 1e3:1e3:9e3, 1e4:5e3:150e3)#10e3:5e3:150e3
old_albedo = 0.2:0.01:0.8

figname = "exp06_AGING-ensemble"
vars2compare = ["Vol", "T"]
metrics2use = ["euclidean", "mape", "cor", "scor", "kcor", "rmsd",
               "rmsdfft",
               "maxperiod",
               "intperiod_80_120", "int_window_maxperiod", "weighted_period",
                "filtered_cor_60to120", "filtered_cor_80to120"]
metric_label = Dict("euclidean" => L"euclidean", "mape" => "mape", "cor" => L"cor", "scor" => L"scor", "kcor" => L"kcor", "rmsd" => "rmsd",
                     "rmsdfft" => "rmsdfft", "maxperiod" => "maxperiod", "weighted_period" => "weighted_period",
                     "intperiod_80_120" => "intperiod_80_120", "int_window_maxperiod" => "int_window_maxperiod",
                     "filtered_cor_60to120" => "filtered_cor_60to120", "filtered_cor_80to120" => "filtered_cor_80to120")

for v in eachindex(vars2compare), m in eachindex(metrics2use)

    var2compare = vars2compare[v]
    metric = metrics2use[m] 

    if var2compare == "T"
        barker2011_temp = JLD2.load_object("data/paleo_records/barker-etal_2011.jld2")
        ref_t, ref_x = interp_time_series(barker2011_temp[1, :], barker2011_temp[2, :], -800e3, 0, 1e3)
    elseif var2compare == "Vol"
        spratt2016_vol = JLD2.load_object("data/paleo_records/spratt-lisiecki_2016.jld2")
        ref_t, ref_x = interp_time_series(spratt2016_vol[1, :], spratt2016_vol[2, :], -800e3, 0, 1e3)
    end
    mint, maxt = ref_t[1], ref_t[end] # time span to compare with
    ref_G, ref_periods = create_PSD(ref_t, ref_x)

    if isfile("data/AGING_ensemble_$(var2compare)_$(metric).jld2")
        mat = JLD2.load_object("data/AGING_ensemble_$(var2compare)_$(metric).jld2")
    else
        mat = Matrix{Float32}(undef, length(taus), length(old_albedo))
        for i in eachindex(taus), j in eachindex(old_albedo)
            df = NCDataset("data/runs/exp06/AGING-DiagTauAlb/AGING-DiagTauAlb_$(string(i))_$(string(j)).nc")

            t, x = df["time"][:], df[var2compare][:]

            if var2compare == "T"
                x = x .- df["Tref"][:]
            end
            
            new_t = cut_run(t, t, mint, maxt)
            new_x = cut_run(t, x, mint, maxt)

            if metric == "euclidean"    # euclidean distance
                mat[i, j] = sqrt(sum((new_x .- ref_x) .^ 2))
            elseif metric == "mape"     # mean absolute percentage error
                mat[i, j] = 100 .* mean(abs.((ref_x .- new_x) ./ (max.(abs.(ref_x), 1e-3))))
            elseif metric == "cor"      # Pearson correlation coefficient
                mat[i, j] = cor(new_x, ref_x)
            elseif metric == "scor"      # Spearman correlation coefficient
                mat[i, j] = corspearman(new_x, ref_x)
            elseif metric == "kcor"      # Kendall tau, correlation coefficient
                mat[i, j] = corkendall(new_x, ref_x)
            elseif metric == "rmsd"
                mat[i, j] = rmsd(new_x, ref_x)
            elseif metric == "rmsdfft"
                G, periods = create_PSD(new_t, new_x)
                mat[i, j] = rmsd(G, ref_G)
            elseif metric == "maxperiod"
                G, periods = create_PSD(new_t, new_x)
                mat[i, j] = periods[findmax(G)[2]] / 1e3
            elseif  metric == "intperiod_80_120"
                G, periods = create_PSD(new_t, new_x)
                idx80, idx120 = findmin(abs.(periods .- 80e3))[2], findmin(abs.(periods .- 120e3))[2]
                mat[i, j] = sum(G[idx120:idx80])
            elseif metric == "filtered_cor_60to120"
                G, periods = create_PSD(new_t, new_x)
                if 60 <= periods[findmax(G)[2]] / 1e3 <= 120
                    mat[i, j] = cor(new_x, ref_x)
                else
                    mat[i, j] = NaN32
                end
            elseif metric == "filtered_cor_80to120"
                G, periods = create_PSD(new_t, new_x)
                if 80 <= periods[findmax(G)[2]] / 1e3 <= 120
                    mat[i, j] = cor(new_x, ref_x)
                else
                    mat[i, j] = NaN32
                end
            elseif metric == "int_window_maxperiod"
                G, periods = create_PSD(new_t, new_x)

                vals, windows = [], [(10, 30), (30, 50), (80, 120)]
                for w in eachindex(windows)
                    idx_lb, idx_ub = findmin(abs.(periods .- windows[w][1] * 1e3))[2], findmin(abs.(periods .- windows[w][2]* 1e3))[2]
                    push!(vals, sum(G[idx_ub:idx_lb]))
                end
                mat[i, j] = findmax(vals)[2]
            elseif metric == "weighted_period"
                G, periods = create_PSD(new_t, new_x)
                weights = ref_G ./ maximum(ref_G)
                mat[i, j] = periods[findmax(G .* weights)[2]] / 1e3
            end

        end
        JLD2.save_object("data/AGING_ensemble_$(var2compare)_$(metric).jld2", mat)
    end

    # Plot
    fontsize = 28
    fig = Figure(resolution=(750, 750), fonts=(; regular="TeX"), fontsize=fontsize)
    linewidth = 2.3

    colormap = cgrad(:bamako, Int(floor(length(old_albedo)/2)), categorical=true)
    colormap = cgrad(:tab20b, Int(floor(length(old_albedo)/2)), categorical=true)

    ax = Axis(fig[2, 1], ylabel=L"$\alpha_\mathrm{oi}$", xlabel=L"$\tau_{\alpha}$ (kyr)",
              xgridvisible=false, ygridvisible=false, yaxisposition=:left, yminorticksvisible=true, xminorticksvisible=true, aspect=1)
    ax2 = Axis(fig[2, 1], ylabel=L"$\alpha_\mathrm{ni}$ - $\alpha_\mathrm{oi}$", yaxisposition=:right, yreversed=true,
              xgridvisible=false, ygridvisible=false, yminorticksvisible=true, xminorticksvisible=true, aspect=1)
    ylims!(ax, (0.2, 0.8))
    ylims!(ax2, 0.9 .- (0.2, 0.8))
    ax.yticks = ([0.2, 0.4, 0.6, 0.8], convert_strings_to_latex([0.2, 0.4, 0.6, 0.8]))
    ax2.yticks = (0.9 .- [0.2, 0.4, 0.6, 0.8], convert_strings_to_latex(round.(0.9 .- [0.2, 0.4, 0.6, 0.8], digits=1)))
    ax.yminorticks = old_albedo
    ax2.yminorticks = 0.9 .- old_albedo

    ax.xticks = ([1, 10, 100, 150], convert_strings_to_latex([1, 10, 100, 150]))
    hidexdecorations!(ax2)
    ax.xminorticks = taus ./ 1e3

    # ax.yticks = (log10.([1e-7, 1e-6]), convert_strings_to_latex(["10^{-7}", "10^{-6}"]))
    # ax.yminorticks = log10.(vcat(5e-8:1e-8:9e-8, 1e-7:1e-7:1e-6))
    # ax.xticks = (log10.([1e-7, 1e-6, 1e-5, 1e-4, 1e-3]), convert_strings_to_latex(["10^{-7}", "10^{-6}", "10^{-5}", "10^{-4}", "10^{-3}"]))
    # ax.xminorticks = log10.(vcat(1e-7:1e-7:9e-7, 1e-6:1e-6:9e-6, 1e-5:1e-5:9e-5, 1e-4:1e-4:1e-3))
    hm = heatmap!(ax, taus ./ 1e3, old_albedo, mat, colormap=colormap)
    heatmap!(ax2, taus ./ 1e3, 0.9 .- old_albedo, mat .* NaN)


    Colorbar(fig[1, 1], hm, width=Relative(2.3 / 3),
        ticklabelsize=fontsize, labelsize=fontsize,
        label=metric_label[metric],
        #ticks=([-1000, -500, 0, +500, +1000], [L"-1000$\,$", L"-500$\,$", L"0$\,$", L"500$\,$", L"1000$\,$"]),
        vertical=false
    )

    CairoMakie.trim!(fig.layout)
    save("figures/aging_ensemble_$(var2compare)_$(metric).png", fig)
    i_max, j_max = findmax(mat)[2][1], findmax(mat)[2][2]
    i_min, j_min = findmin(mat)[2][1], findmin(mat)[2][2]
    display((var2compare, metric, findmax(mat), findmin(mat)))
end



# dfmax = NCDataset("data/runs/exp06/AGING-DiagTauAlb/AGING-DiagTauAlb_$(string(i_max))_$(string(j_max)).nc")
# paramsmax = JLD2.load_object("data/runs/exp06/AGING-DiagTauAlb/AGING-DiagTauAlb_$(string(i_max))_$(string(j_max))_params.jld2")
# dfmin = NCDataset("data/runs/exp06/AGING-DiagTauAlb/AGING-DiagTauAlb_$(string(i_min))_$(string(j_min)).nc")
# paramsmin = JLD2.load_object("data/runs/exp06/AGING-DiagTauAlb/AGING-DiagTauAlb_$(string(i_min))_$(string(j_min))_params.jld2")
# display((paramsmax.tau_albedo, paramsmax.albedo_oldice, paramsmin.tau_albedo, paramsmin.albedo_oldice))

# fig = Figure()
# ax = Axis(fig[1, 1])
# lines!(ax, ref_t, ref_x)
# lines!(ax, dfmax["time"][:], dfmax[var2compare][:])
# lines!(ax, dfmin["time"][:], dfmin[var2compare][:])
# fig