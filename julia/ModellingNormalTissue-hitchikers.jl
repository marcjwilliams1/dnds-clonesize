using StemCellModels
using DataFrames
using ProgressMeter
using StatsBase
using Statistics
using Random
using StatsBase
using ArgParse
using CSV
using Distributions

s = ArgParseSettings()
@add_arg_table s begin
    "--resultsfile"
        help = "File to save results"
end

parsed_args = parse_args(ARGS, s)
print(parsed_args)

function simulatepopulation(;Δ = 0.0, ρ = 100.0, Amin = 1 ./ ρ, gap = 0.02,
        tend = 20.0, N0 = 1000, r = 0.5, λ = 1.0,
        sample = false, nsample = 10, Amax = 0.0, Nits = 50,
        μp = 0.001, μd = 0.001)

    rangeN = 1:10^4
    Cntd = zeros(length(rangeN), Nits)
    Cntp = zeros(length(rangeN), Nits)
    freq_sel = Float64[]
    vaf_sel = Float64[]
    freq_neut = Float64[]
    freq_neutall = Float64[]
    freq_hitchikers = Float64[]
    vaf_neut = Float64[]
    Nend = Int64[]
    SMtog = StemCellModel(N0, Δmut = Δ, μp = μp, μd = μd, tend = tend, r = r, λ = λ)
    for i in 1:Nits
        scst = runsimulation(SMtog, progress = false, finish = "time", onedriver = true, restart = true)
        #println(scst)
        append!(freq_sel, filter!(x -> x > 0, scst.mutationsize_d))
        append!(vaf_sel, filter!(x -> x > 0, scst.mutationfrequencies_d))
        append!(vaf_neut, filter!(x -> x > 0, scst.mutationfrequencies_p))
        append!(freq_neutall, filter!(x -> x > 0, scst.mutationsize_p))
        x = StemCellModels.gethitchikers(scst)
        append!(freq_neut, x[1])
        append!(freq_hitchikers, x[2])
        Cntd[:, i] = counts(sort(convert(Array{Int64, 1}, scst.mutationsize_d)), rangeN)
        push!(Nend, length(scst.stemcells))
    end

    return [freq_sel, freq_neutall], Nend, Cntd, [vaf_sel, vaf_neut], [freq_neut, freq_hitchikers]
end


myDF = DataFrame([Int64, Int64, Int64, Int64, Int64, Int64, Int64, Float64, Float64, Float64], [:CnNS, :CnS, :CnSI, :CnSH, :n, :Nsims, :N0, :muNS, :muS, :tend], 0)


for m in [0.0001, 0.001, 0.01]
    println("Mutation rate = $m")
    for t in [21.5, 25.5, 37.5, 45.5, 49.5, 53.5, 57.5, 69.5, 73.5]
        println("Time = $t")
        delta = 0.1
        Nsims = 10^4
        tend = t
        μd = m
        μp = m ./ 3
        N0 = 10^3
        freq, nend, Cnt, vaf, h = simulatepopulation(Nits = Nsims, Δ = delta, tend = tend, μd = μd, μp = μp, N0 = N0);
        nmax = 10000
        DF = DataFrame(CnNS = counts(convert(Array{Int64, 1}, freq[1]), nmax),
                       CnS = counts(convert(Array{Int64, 1}, freq[2]), nmax),
                       CnSI = counts(convert(Array{Int64, 1}, h[1]), nmax),
                       CnSH = counts(convert(Array{Int64, 1}, h[2]), nmax),
                       n = 1:nmax)
        DF[:Nsims] = Nsims
        DF[:N0] = N0
        DF[:muNS] = μd
        DF[:muS] = μp
        DF[:tend] = tend
        append!(myDF, DF)
    end
end

CSV.write(parsed_args["resultsfile"], myDF)
