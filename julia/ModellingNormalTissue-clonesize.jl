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
using Dates


println("Starting simulations")

s = ArgParseSettings()
@add_arg_table s begin
    "--resultsfile"
        help = "File to save results"
    "--resultsfile2"
        help = "File to save results"
end

parsed_args = parse_args(ARGS, s)
print(parsed_args)

Random.seed!(123)

#function to simulate population of cells with different paramters
function simulatepopulation(;Δ = 0.0,
        tend = 20.0, Ncells = 1000, r = 0.5, λ = 1.0,
        sample = false, nsample = 10, Nits = 50,
        μp = 0.001, μd = 0.001)

    rangeN = 1:10^4
    Cntd = zeros(length(rangeN), Nits)
    Cntp = zeros(length(rangeN), Nits)
    freq_mut = Float64[]
    SMtog = StemCellModel(Ncells, Δmut = Δ, μp = μp, μd = μd, tend = tend, r = r, λ = λ)
    for i in 1:Nits
        scst = runsimulation(SMtog, progress = false, finish = "time", onedriver = true, restart = true, maxcells = 10^6)
        append!(freq_mut, filter!(x -> x > 0, scst.mutationsize_d))
    end

    return freq_mut
end

myDF = DataFrame([Float64, Float64, Float64, Float64, Float64, Int64, String, Float64], [:f, :delta, :t, :rlam, :mu, :Nsims, :gene, :N0], 0)

for Δ in [0.0, 0.025, 0.05, 0.1]
    for t in [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
        println("$(Dates.format(now(), "HH:MM")  ): Δ = $Δ, t = $t")
        flush(stdout)
        flush(stderr)
        Nsims = 1000
        mu = rand(TruncatedNormal(0.001, 0.0005, 0.0, 10.0))
        N0 = 1000
        freq = simulatepopulation(;Δ = Δ, tend = t,
                                    λ = 0.5, r = 1.0,
                                    Nits = Nsims, μd = mu,
                                    μp = 0.0,
                                    Ncells = N0)
        DFtemp = DataFrame(f = freq, delta = Δ, t = t, rlam = 0.5, mu = mu, Nsims = Nsims, gene = "gene_$(randstring(5))", N0 = N0)
        append!(myDF, DFtemp)
    end
end
#CSV.write("file.csv", myDF)

CSV.write(parsed_args["resultsfile"], myDF)




#function to simulate population of cells with different paramters
function simulatepopulation(;Δ = 0.0,
        tend = 20.0, Ncells = 1000, r = 0.5, λ = 1.0,
        sample = false, nsample = 10, Nits = 50,
        μp = 0.001, μd = 0.001)

    rangeN = 1:10^4
    Cntd = zeros(length(rangeN), Nits)
    Cntp = zeros(length(rangeN), Nits)
    freq_mut = Float64[]
    Δvec = Float64[]
    for i in 1:Nits
        if Δ == 0.0
            Δsample = 0.0
        else
            #Δsample = rand(Beta(0.3, 0.3 * ((1 - Δ) / Δ)))
            Δtemp = rand([0.05])
            if Δtemp == 0.0
                Δsample = 0.0
            else
                Δsample = rand(Exponential(Δtemp))
                #Δsample = Δtemp
            end
        end
        println("Simulation $(i), Δ = $(Δsample)")
        flush(stdout)
        if Δsample > 0.7
            continue
        end
        SMtog = StemCellModel(Ncells, Δmut = Δsample, μp = μp, μd = μd, tend = tend, r = r, λ = λ)
        @time scst = runsimulation(SMtog, progress = false, finish = "time", onedriver = true, restart = true, maxcells = 5 * 10^6)
        append!(freq_mut, filter!(x -> x > 0, scst.mutationsize_d))
        push!(Δvec, Δsample)
    end

    return freq_mut, mean(Δvec)
end

myDF = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64, Int64, String, Float64], [:f, :delta, :deltasample, :t, :rlam, :mu, :Nsims, :gene, :N0], 0)

for Δ in [0.0, 0.1]
    for t in [20.0, 30.0, 40.0, 50.0, 60.0, 70.0]
        println("$(Dates.format(now(), "HH:MM")  ): Δ = $Δ, t = $t\n")
        flush(stdout)
        flush(stderr)
        Nsims = 1000
        mu = rand(TruncatedNormal(0.001, 0.0005, 0.0, 10.0))
        N0 = 1000
        freq, del = simulatepopulation(;Δ = Δ, tend = t,
                                    λ = 0.5, r = 1.0,
                                    Nits = Nsims, μd = mu,
                                    μp = 0.0,
                                    Ncells = N0)
        DFtemp = DataFrame(f = freq, delta = Δ, deltasample = del, t = t, rlam = 0.5, mu = mu, Nsims = Nsims, gene = "gene_$(randstring(5))", N0 = N0)
        append!(myDF, DFtemp)
    end
end
#CSV.write("file.csv", myDF)

CSV.write(parsed_args["resultsfile2"], myDF)
