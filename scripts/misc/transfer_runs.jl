using Pkg
Pkg.activate(".")

locpacco = "/home/sergio/entra/models/pacco_vers/pacco_v0.6/"
locdir = "/home/sergio/entra/proyects/d05_paper-PACCO/Perez-Montero-etal_YYYY_ESD/"

locexps = "$(locpacco)/output/Perez-Montero-etal_2024_ESD/"

function get_runs(exp_path::String, subdirs::Vector)
    runs = []
    for s in subdirs
        push!(runs, readdir("$(exp_path)/$(s)"))
    end
    return runs
end

function format_runs(runs)

    return new_runs, new_params
end

function move_runs!(runs, path2runs::String, path2store::String)
    for r in eachindex(runs)
        run = runs[r]
        cp("$(path2runs)/$(run)/pacco.nc", "$(path2store)/$(run).nc", force=true)
        cp("$(path2runs)/$(run)/params.jld2", "$(path2store)/$(run)_params.jld2", force=true)
    end

    return nothing
end

function transfer_runs_from_pacco!(exp_path::String, subdirs::Vector, path2data::String)
    runs = get_runs(exp_path, subdirs)

    # Check if experiment directory exists
    isdir(path2data) || mkdir(path2data)

    for s in eachindex(subdirs)
        # Create if necessary the subdirectories
        isdir("$(path2data)/$(subdirs[s])/") || mkdir("$(path2data)/$(subdirs[s])/")

        # Move runs from PACCO local folder to the desired directory
        move_runs!(runs[s], "$(exp_path)/$(subdirs[s])/", "$(path2data)/$(subdirs[s])/")
    end

    return nothing
end

# exp01, LIN 
transfer_runs_from_pacco!("$(locexps)/exp01/", ["LIN-MV"], "$(locdir)/data/runs/exp01/")

# exp02, NONLIN 
transfer_runs_from_pacco!("$(locexps)/exp02/", ["NONLIN-MV"], "$(locdir)/data/runs/exp02/")

# exp03, ISOS 
transfer_runs_from_pacco!("$(locexps)/exp03/", ["ISOS-MV", "RISOS-MV"], "$(locdir)/data/runs/exp03/")

# exp04, BASE 
transfer_runs_from_pacco!("$(locexps)/exp04/", ["BASE-MV-Cs"], "$(locdir)/data/runs/exp04/")

# exp05, THERM 
transfer_runs_from_pacco!("$(locexps)/exp05/", ["THERM-MV-Cs", "THERM-MV-taukin"], "$(locdir)/data/runs/exp05/")

# exp06, AGING
transfer_runs_from_pacco!("$(locexps)/exp06/", ["AGING-MV-Cs", "AGING-MV-taualpha"], "$(locdir)/data/runs/exp06/")