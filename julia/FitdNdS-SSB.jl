using Pkg
#Pkg.build("RCall")
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
    "--oesophagusdndsdata_miss"
        help = "Oesophagus dnds values - missense"
    "--oesophagusdndsdata_non"
        help = "Oesophagus dnds values - nonsense"
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

println("Reading in per gene data...")
DF = rcopy(readr.read_csv(parsed_args["oesophagusdndsdata"]))
DF[:A] = 2 .* DF[:cutoff];
DFnon = rcopy(readr.read_csv(parsed_args["oesophagusdndsdata_non"], guess_max = 10^5))
DFnon[:A] = 2 .* DFnon[:cutoff];
DFmis = rcopy(readr.read_csv(parsed_args["oesophagusdndsdata_miss"], guess_max = 10^5))
DFmis[:A] = 2 .* DFmis[:cutoff];
println("Data read in")

myDFall = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, String, Float64, String, Float64],
    [:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :Age2, :gene, :nmutations], 0)
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

println("Analyse per patient per gene data for all mutations...")
for gene in unique(DF[:gene])
    println("Analysing gene: $gene")
    DFgene = filter(row -> row[:gene] == gene, DF);
    println("Analysing all mutations")
    for p in unique(DF[:ID])
        DFpatient = filter(row -> row[:ID] == p, DFgene);
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
        append!(myDFall, x.DF)
    end
end

println("Number of fits: $(length(unique(myDFall[:gene])))")


println("Analyse per patient per gene data for nonsense mutations...")
for gene in unique(DFnon[:gene])
    println("Analysing gene: $gene")
    DFgene = filter(row -> row[:gene] == gene, DFnon);
    println("Analysing nonsense mutations")
    println(DFgene)
    for p in unique(DFnon[:ID])
        DFpatient = filter(row -> row[:ID] == p, DFgene);
        age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]

        if length(DFpatient[:muts]) == 0
            println("Missing fit: $p, $gene")
            continue
        end

        nmuts = DFpatient[:muts][end]
        println(nmuts)
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
        append!(myDFnon, x.DF)
    end
end

println("Number of fits: $(length(unique(myDFnon[:gene])))")


println("Analyse per patient per gene data for missense mutations...")
for gene in unique(DF[:gene])
    println("Analysing gene: $gene")
    DFgene = filter(row -> row[:gene] == gene, DFmis);
    println("Analysing missense mutations")
    for p in unique(DFmis[:ID])
        DFpatient = filter(row -> row[:ID] == p, DFgene);
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
    end
end

println("Number of fits: $(length(unique(myDFmiss[:gene])))")

myDFall[:mutation_class] = "All"
myDFmiss[:mutation_class] = "Missense"
myDFnon[:mutation_class] = "Nonsense"

append!(myDFall, myDFmiss)
append!(myDFall, myDFnon)

@rput myDFall
R"""
library(readr)
write_csv(myDFall, $(parsed_args["oesophagusfit"]))
""";
