using Pkg
Pkg.activate(".")
using NCDatasets, JLD2
include("/home/sergio/entra/proyects/d05_paper-PACCO/scripts/misc/tools.jl")

## Berends et al., 2020
convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/Vol/berends-etal_2020.nc", -8.0e5, 0.0, 1e3)
mv("data/paleo_records/berends-etal_2020.jld2", "data/paleo_records/berends-etal_2020_Vol_800.jld2", force=true)

convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/C/berends-etal_2020.nc", -8.0e5, 0.0, 1e3)
mv("data/paleo_records/berends-etal_2020.jld2", "data/paleo_records/berends-etal_2020_C_800.jld2", force=true)

## Yamamoto et al., 2022
convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/C/yamamoto-etal_2022.nc", -8.0e5, 0.0, 1e3)
mv("data/paleo_records/yamamoto-etal_2022.jld2", "data/paleo_records/yamamoto-etal_C_800.jld2", force=true)

## Lisiecki and Raymo, 2005
convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/d18O/lisiecki-raymo_2005.nc", -8.0e5, 0.0, 1e3)
mv("data/paleo_records/lisiecki-raymo_2005.jld2", "data/paleo_records/lisiecki-raymo_2005_d18O_800.jld2", force=true)

convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/T/lisiecki-raymo_2005.nc", -8.0e5, 0.0, 1e3)
mv("data/paleo_records/lisiecki-raymo_2005.jld2", "data/paleo_records/lisiecki-raymo_2005_T_800.jld2", force=true)

## Bintanja and van de Wal, 2008
convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/Vol/bintanja-vandewal_2008.nc", -8.0e5, 0.0, 1e3)
mv("data/paleo_records/bintanja-vandewal_2008.jld2", "data/paleo_records/bintanja-vandewal_2008_800.jld2", force=true)

convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/T/bintanja-vandewal_2008.nc", -8.0e5, 0.0, 1e3)
mv("data/paleo_records/bintanja-vandewal_2008.jld2", "data/paleo_records/bintanja-vandewal_2008_T_800.jld2", force=true)

## Lüthi et al., 2008
convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/C/luthi-etal_2008.nc", -8.0e5, 0.0, 1e3)

## Barker et al., 2011
convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/T/barker-etal_2011.nc", -8.0e5, 0.0, 1e3)

## Spratt and Lisiecki, 2016
convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/Vol/spratt-lisiecki_2016.nc", -8.0e5, 0.0, 1e3)

## Ahn et al., 2017
convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/d18O/ahn-etal_2017.nc", -8.0e5, 0.0, 1e3)
mv("data/paleo_records/ahn-etal_2017.jld2", "data/paleo_records/ahn-etal_2017_800.jld2", force=true)

## Hodell et al., 2023
convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/T/hodell-etal_2023.nc", -8.0e5, 0.0, 1e3)
mv("data/paleo_records/hodell-etal_2023.jld2", "data/paleo_records/hodell-etal_2023_T_800.jld2", force=true)

convert_nc_to_jld2("/home/sergio/entra/ice_data/dataForPACCO/d18O/hodell-etal_2023.nc", -8.0e5, 0.0, 1e3)
mv("data/paleo_records/hodell-etal_2023.jld2", "data/paleo_records/hodell-etal_2023_d18O_800.jld2", force=true)