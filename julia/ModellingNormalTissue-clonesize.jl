using StemCellModels
using DataFrames
using ProgressMeter
using StatsBase
using Statistics
using Random
using StatsBase
using ArgParse
using CSV

s = ArgParseSettings()
@add_arg_table s begin
    "--resultsfile"
        help = "File to save results"
end

parsed_args = parse_args(ARGS, s)


#function to simulate population of cells with different paramters
function simulatepopulation(;Δ = 0.0, ρ = 100.0, Amin = 1 ./ ρ, gap = 0.02,
        tend = 20.0, Ncells = 1000, r = 0.5, λ = 1.0,
        sample = false, nsample = 10, Amax = 0.0, Nits = 50,
        μp = 0.001, μd = 0.001)

    rangeN = 1:10^4
    Cntd = zeros(length(rangeN), Nits)
    Cntp = zeros(length(rangeN), Nits)
    freq_sel = Float64[]
    freq_neut = Float64[]
    SMtog = SkinStemCellModel(Ncells, Δmut = Δ, μp = μp, μd = μd, tend = tend, r = r, λ = λ)
    @showprogress for i in 1:Nits
        scst = runsimulation(SMtog, progress = false, finish = "time", onedriver = true, restart = true)
        #println(scst)
        append!(freq_sel, filter!(x -> x > 0, scst.mutationsize_d))
        append!(freq_neut, filter!(x -> x > 0, scst.mutationsize_p))
    end

    return freq_sel
end

myDF = DataFrame([Float64, Float64, Float64, Float64, Float64, Int64, String, Float64], [:f, :delta, :t, :rlam, :mu, :Nsims, :gene, :N0], 0)

for Δ in [0.0, 0.01, 0.015, 0.02, 0.025, 0.05]
    for t in [10.0, 20.0, 30.0, 40.0]
        println("Δ = $Δ, t = $t")
        Nsims = 1000
        mu = 0.001
        N0 = 1000
        freq = simulatepopulation(;Δ = Δ, tend = t, Amin = 0.05,
                                    ρ = 100.0, λ = 0.5, Amax = 25.0, r = 1.0,
                                    Nits = Nsims, μd = mu, μp = 0.0,
                                    Ncells = N0)
        DFtemp = DataFrame(f = freq, delta = Δ, t = t, rlam = 0.5, mu = mu, Nsims = Nsims, gene = "gene_$(Δ)_$(t)", N0 = N0)
        append!(myDF, DFtemp)
    end
end

CSV.write(parsed_args[:resultsfile], myDF)
