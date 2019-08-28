using Pkg
Pkg.build("RCall")
using CSV
using DataFrames
using ArgParse
using RCall

@rimport readr as readr

# load functions to fit data via least squares
include("optim.jl")

#parse arguments
s = ArgParseSettings()
@add_arg_table s begin
    "--oesophagusdndsdata"
        help = "Oesophagus dnds values"
    "--oesophagusmetadata"
        help = "Information on patient ages etc"
    "--oesophagusfit"
        help = "Fits for oesophagus missense mutations"
end

parsed_args = parse_args(ARGS, s)
println(parsed_args)

println("Read in data...")
println(parsed_args["oesophagusmetadata"])

DFdonor = rcopy(readr.read_csv(parsed_args["oesophagusmetadata"]))
DFcohort = rcopy(readr.read_csv(parsed_args["oesophagusdndsdata"]))
DFcohort[:A] = 2 * DFcohort[:areacutoff];

println("Reading in per gene data...")
DF = rcopy(readr.read_csv(parsed_args["oesophagusdndsdatagenes"]))
DF[:A] = 2 .* DF[:cutoff];
println("Data read in")

delta = Float64[]
rlam = Float64[]
genevec = String[]
mutationtypevec = String[]
rsqvec = Float64[]
patientvec = String[]
myDFmiss = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, String, Float64, String, Float64],
    [:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :Age2, :gene, :nmutations], 0)
myDFnon = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, String, Float64, String, Float64],
    [:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :Age2, :gene, :nmutations], 0)

println("Analyse per patient per gene data...")
for gene in unique(DF[:gene])
    println("Analysing gene: $gene")
    DFgene = filter(row -> row[:gene] == gene, DF);
    println("Analysing missense mutations")
    for p in unique(DF[:ID])
        DFpatient = filter(row -> row[:patient] == p, DFgene);
        age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]

        nmuts = DFpatient[:muts][end]
        if nmuts < 5
            println("Patient $p, only $(nmuts) mutations, moving to next patient")
            continue
        end
        println("Analysing patient $p, $(nmuts) mutations")

        x = LLoptimizationresults(DFpatient[:dnds], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
        x.DF[:patient] = p
        x.DF[:Age] = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
        x.DF[:Age2] = age
        x.DF[:gene] = gene
        x.DF[:nmutations] = nmuts
        append!(myDFmiss, x.DF)
        push!(delta, x.Δ)
        push!(rlam, x.rλ)
        push!(genevec, gene)
        push!(mutationtypevec, "All")
        push!(rsqvec, x.DF[:rsq][1])
        push!(patientvec, p)
    end
end

println("Number of fits: $(length(unique(myDFmiss[:gene])))")

@rput myDFmiss
@rput myDFnon
R"""
library(readr)
write_csv(myDFmiss, $(parsed_args["oesophagusfits"]))
""";
