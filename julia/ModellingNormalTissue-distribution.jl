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
    "--resultsfile1"
        help = "File 1 to save results"
    "--resultsfile2"
        help = "File 2 to save results"
end

Random.seed!(123)

parsed_args = parse_args(ARGS, s)
print(parsed_args)
include("optim.jl")

function intervaldnds_data(Msel, Mneut, mmin, Nmax; mgap = 5.0, mup = 1.0, mud = 1.0)
    #calculate cumulative dn/ds across cohort, normalizing for mutation rates
    dndscumsum = Float64[]
    dnvec = Float64[]
    dsvec = Float64[]
    driverVAFt = filter(x -> x >= mmin, Msel)
    passengerVAFt = filter(x -> x >= mmin, Mneut)
    freqs = collect(mmin:mgap:Nmax)
    for f in freqs
        dn = length(filter(x -> x< f, driverVAFt))/mud
        ds = length(filter(x -> x< f, passengerVAFt))/mup
        push!(dndscumsum, dn/ds)
        push!(dnvec, dn)
        push!(dsvec, ds)
    end
    return dndscumsum[2:end], freqs[2:end], dnvec[2:end], dsvec[2:end]
end

#function to simulate population of cells with different paramters
function simulatepopulation(;Δ = 0.0, ρ = 100.0, Amin = 1 ./ ρ, gap = 0.02, tend = 20.0,
        Ncells = 1000, r = 0.5, λ = 1.0,
        sample = false, nsample = 10, Amax = 0.0)

    rangeN = 1:10^4
    Nits = 100
    Cntd = zeros(length(rangeN), Nits)
    Cntp = zeros(length(rangeN), Nits)
    Msel = Int64[]
    Mneut = Int64[]
    Δvec = Float64[]
    SMtog = StemCellModel(Ncells, Δmut = Δ, μp = 0.001, μd = 0.001, tend = tend, r = r, λ = λ)
    for i in 1:Nits
        Δsample = rand(Exponential(Δ))
        push!(Δvec, Δsample)
        SMtog = StemCellModel(Ncells, Δmut = Δsample, μp = 0.001, μd = 0.001, tend = tend, r = r, λ = λ)
        scst = runsimulation(SMtog, progress = false, finish = "time", onedriver = true, restart = true, maxcells = 10^6)
        append!(Msel, filter!(x -> x > 0, counts(sort(StemCellModels.cellsconvert(scst.stemcells)[2]), rangeN)))
        append!(Mneut, filter!(x -> x > 0, counts(sort(StemCellModels.cellsconvert(scst.stemcells)[1]), rangeN)))
    end

    Asel = Msel ./ ρ
    Aneut = Mneut ./ ρ
    Asel = filter(x -> x >= Amin, Asel)
    Aneut = filter(x -> x >= Amin, Aneut)

    if sample
        prevlength = length(Asel)
        Asel = Asel[StatsBase.sample(1:length(Asel), nsample, replace = false)]
        Aneut = Aneut[StatsBase.sample(1:length(Aneut),
                Int64(round(nsample * (length(Aneut) / prevlength))), replace = false)]
    end

    if Amax == 0.0
        Amax = maximum(Asel)
    end

    x1 = intervaldnds_data(Asel, Aneut, Amin, Amax; mgap = gap, mup = SMtog.μp, mud = SMtog.μd)

    DFinterval = DataFrame(dnds = x1[1], A = x1[2], dn = x1[3] * SMtog.μd / Nits,
    ds = x1[4] * SMtog.μd / Nits)

    return DFinterval, SMtog, [DFinterval[:dn][end] * Nits, DFinterval[:ds][end] * Nits], mean(Δvec)
end


println("Generating some example data with fits")

#create empty dataframe to store data
myDF = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Int64],
[:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda, :rsq, :deltatrue,
    :deltatruesample, :lambdartrue, :simindex], 0)

global j = 0
for i in 1:1000
    print(i)
    for Δ in [0.05, 0.1, 0.2]
        println(Δ)
        for rlam in [0.25]
            global j = j + 1
            println(rlam)
            rlambda = rlam
            DF1, SM, x1, Δsample = simulatepopulation(;Δ = Δ, tend = 30.0, Amin = 0.05, ρ = 100.0, λ = rlambda, Amax = 25.0)
            x = LLoptimizationresults(DF1[:dnds], DF1[:A]; t = 30.0, Amin = 0.05, ρ = 100.0)
            x.DF[:deltatrue] = Δ
            x.DF[:deltatruesample] = Δsample
            x.DF[:lambdartrue] = rlambda * SM.r
            x.DF[:simindex] = j
            append!(myDF, x.DF)
        end
    end
end

println("Writing data to file")
@rput myDF
R"""
res <- dplyr::distinct(myDF, deltafit, lambdarfit, deltatrue, deltatruesample, simindex, rsq)
readr::write_csv(res, $(parsed_args["resultsfile1"]))
""";




#function to simulate population of cells with different paramters
function simulatepopulation(;Δ = 0.0, ρ = 100.0, Amin = 1 ./ ρ, gap = 0.02, tend = 20.0,
        Ncells = 1000, r = 0.5, λ = 1.0,
        sample = false, nsample = 10, Amax = 0.0)

    rangeN = 1:10^4
    Nits = 100
    Cntd = zeros(length(rangeN), Nits)
    Cntp = zeros(length(rangeN), Nits)
    Msel = Int64[]
    Mneut = Int64[]
    Δvec = Float64[]
    SMtog = StemCellModel(Ncells, Δmut = Δ, μp = 0.001, μd = 0.001, tend = tend, r = r, λ = λ)
    for i in 1:Nits
        Δsample = rand(Beta(0.3, 0.3 * ((1 - Δ) / Δ)))
        Δsample = Δ
        push!(Δvec, Δsample)
        SMtog = StemCellModel(Ncells, Δmut = Δsample, μp = 0.001, μd = 0.001, tend = tend, r = r, λ = λ)
        scst = runsimulation(SMtog, progress = false, finish = "time", onedriver = true, restart = true, maxcells = 10^6)
        append!(Msel, filter!(x -> x > 0, counts(sort(StemCellModels.cellsconvert(scst.stemcells)[2]), rangeN)))
        append!(Mneut, filter!(x -> x > 0, counts(sort(StemCellModels.cellsconvert(scst.stemcells)[1]), rangeN)))
    end

    Asel = Msel ./ ρ
    Aneut = Mneut ./ ρ
    Asel = filter(x -> x >= Amin, Asel)
    Aneut = filter(x -> x >= Amin, Aneut)

    if sample
        prevlength = length(Asel)
        Asel = Asel[StatsBase.sample(1:length(Asel), nsample, replace = false)]
        Aneut = Aneut[StatsBase.sample(1:length(Aneut),
                Int64(round(nsample * (length(Aneut) / prevlength))), replace = false)]
    end

    if Amax == 0.0
        Amax = maximum(Asel)
    end

    x1 = intervaldnds_data(Asel, Aneut, Amin, Amax; mgap = gap, mup = SMtog.μp, mud = SMtog.μd)

    DFinterval = DataFrame(dnds = x1[1], A = x1[2], dn = x1[3] * SMtog.μd / Nits,
    ds = x1[4] * SMtog.μd / Nits)

    return DFinterval, SMtog, [DFinterval[:dn][end] * Nits, DFinterval[:ds][end] * Nits], mean(Δvec)
end


println("Generating some example data with fits")

#create empty dataframe to store data
myDF = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Int64],
[:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda, :rsq, :deltatrue,
    :deltatruesample, :lambdartrue, :simindex], 0)

global j = 0
for i in 1:1000
    print(i)
    for Δ in [0.05, 0.1, 0.2]
        println(Δ)
        for rlam in [0.25]
            global j = j + 1
            println(rlam)
            rlambda = rlam
            DF1, SM, x1, Δsample = simulatepopulation(;Δ = Δ, tend = 30.0, Amin = 0.05, ρ = 100.0, λ = rlambda, Amax = 25.0)
            x = LLoptimizationresults(DF1[:dnds], DF1[:A]; t = 30.0, Amin = 0.05, ρ = 100.0)
            x.DF[:deltatrue] = Δ
            x.DF[:deltatruesample] = Δsample
            x.DF[:lambdartrue] = rlambda * SM.r
            x.DF[:simindex] = j
            append!(myDF, x.DF)
        end
    end
end

println("Writing data to file")
@rput myDF
R"""
res <- dplyr::distinct(myDF, deltafit, lambdarfit, deltatrue, deltatruesample, simindex, rsq)
readr::write_csv(res, $(parsed_args["resultsfile2"]))
""";
