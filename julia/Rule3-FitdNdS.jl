using CSV
using DataFrames
using ArgParse

# load functions to fit data via least squares
include("optim.jl")

#parse arguments
s = ArgParseSettings()
@add_arg_table s begin
    "--oesophagusdndsdata"
        help = "Oesophagus dnds values"
    "--oesophagusdndsdatagenes"
        help = "Oesophagus dnds values per gene"
    "--oesophagusdndsneutral"
        help = "Oesophagus dnds values"
    "--oesophagusmetadata"
        help = "Information on patient ages etc"
    "--skindndsdata"
        help = "Oesophagus dnds values"
    "--skindndsdatagenes"
        help = "Oesophagus dnds values per gene"
    "--skinmetadata"
        help = "Information on patient ages etc"
    "--oesophagusfitmissense"
        help = "Fits for oesophagus missense mutations"
    "--oesophagusfitnonsense"
        help = "Fits for oesophagus nonsense mutations"
    "--oesophagusfitmissensepergene"
        help = "Fits for oesophagus missense mutations per gene"
    "--oesophagusfitnonsensepergene"
        help = "Fits for oesophagus nonsense mutations per gene"
    "--skinfitmissense"
        help = "Fits for oesophagus missense mutations"
    "--skinfitnonsense"
        help = "Fits for skin nonsense mutations"
    "--skinfitmissensepergene"
        help = "Fits for oesophagus missense mutations per gene"
    "--skinfitnonsensepergene"
        help = "Fits for skin nonsense mutations per gene"
    "--oesophagusfitneutral"
        help = "Fits for oesophagus neutral genes"
end

parsed_args = parse_args(ARGS, s)

DFdonor = CSV.read(parsed_args["oesophagusmetadata"]);
DFcohort = CSV.read(parsed_args["oesophagusdndsdata"]);
DFcohort[:A] = 2 * DFcohort[:areacutoff];

myDFmiss = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, String, Float64, Float64],
[:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :Age2, :nmutations], 0)

for p in DFdonor[:patient]
    println(p)
    DFpatient = filter(row -> row[:patient] == p, DFcohort);
    DFpatient = filter(row -> row[:name] == "wmis", DFpatient)
    if DFpatient[:mle][1] == 0.0
        continue
    end
    age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]
    x = LLoptimizationresults(DFpatient[:mle], DFpatient[:A]; t = age, Amin = 2 * 0.0075, ρ = 5000.0)
    x.DF[:patient] = p
    x.DF[:Age] = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
    x.DF[:Age2] = age
    x.DF[:nmutations] = DFpatient[:nmutations]
    append!(myDFmiss, x.DF)
end

myDFnon = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, String, Float64, Float64],
[:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :Age2, :nmutations], 0)

for p in DFdonor[:patient]
    println(p)
    DFpatient = filter(row -> row[:patient] == p, DFcohort);
    DFpatient = filter(row -> row[:name] == "wnon", DFpatient)
    if DFpatient[:mle][1] == 0.0
        continue
    end
    age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]
    x = LLoptimizationresults(DFpatient[:mle], DFpatient[:A]; t = age, Amin = 2 * 0.0075, ρ = 5000.0)
    x.DF[:patient] = p
    x.DF[:Age] = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
    x.DF[:Age2] = age
    x.DF[:nmutations] = DFpatient[:nmutations]
    append!(myDFnon, x.DF)
end

myDFall = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, String, Float64, Float64],
[:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :Age2, :nmutations], 0)

for p in DFdonor[:patient]
    println(p)
    DFpatient = filter(row -> row[:patient] == p, DFcohort);
    DFpatient = filter(row -> row[:name] == "wall", DFpatient)
    if DFpatient[:mle][1] == 0.0
        continue
    end
    age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]
    x = LLoptimizationresults(DFpatient[:mle], DFpatient[:areacutoff]; t = age, Amin = 0.01, ρ = 5000.0)
    x.DF[:patient] = p
    x.DF[:Age] = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
    x.DF[:Age2] = age
    x.DF[:nmutations] = DFpatient[:nmutations]
    append!(myDFall, x.DF)
end


@rput myDFmiss
@rput myDFnon
R"""
library(readr)
write_csv(myDFmiss, $(parsed_args["oesophagusfitmissense"]))
write_csv(myDFnon, $(parsed_args["oesophagusfitnonsense"]))
""";

DF = CSV.File(parsed_args["oesophagusdndsdatagenes"]) |> DataFrame
DF[:A] = 2 * DF[:areacutoff];


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

for gene in unique(DF[:gene_name])
#for gene in ["TP53", "NOTCH1"]
    println(gene)
    DFgene = filter(row -> row[:gene_name] == gene, DF);

    for p in DFdonor[:patient]
        println(p)
        DFpatient = filter(row -> row[:patient] == p, DFgene);
        age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]

        nmuts = DFpatient[:n_syn][end] + DFpatient[:n_mis][end]
        if nmuts < 5
            continue
        end

        x = LLoptimizationresults(DFpatient[:wmis_cv], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
        x.DF[:patient] = p
        x.DF[:Age] = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
        x.DF[:Age2] = age
        x.DF[:gene] = gene
        x.DF[:nmutations] = nmuts
        append!(myDFmiss, x.DF)
        push!(delta, x.Δ)
        push!(rlam, x.rλ)
        push!(genevec, gene)
        push!(mutationtypevec, "Missense")
        push!(rsqvec, x.DF[:rsq][1])
        push!(patientvec, p)
    end

    for p in DFdonor[:patient]
        println(p)
        DFpatient = filter(row -> row[:patient] == p, DFgene);
        age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]

        nmuts = DFpatient[:n_syn][end] + DFpatient[:n_non][end]
        if nmuts < 5
            continue
        end

        x = LLoptimizationresults(DFpatient[:wnon_cv], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
        x.DF[:patient] = p
        x.DF[:Age] = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
        x.DF[:Age2] = age
        x.DF[:gene] = gene
        x.DF[:nmutations] = nmuts
        append!(myDFnon, x.DF)
        push!(delta, x.Δ)
        push!(rlam, x.rλ)
        push!(genevec, gene)
        push!(mutationtypevec, "Nonsense")
        push!(rsqvec, x.DF[:rsq][1])
        push!(patientvec, p)
    end
end

@rput myDFmiss
@rput myDFnon
R"""
library(readr)
write_csv(myDFmiss, $(parsed_args["oesophagusfitmissensepergene"]))
write_csv(myDFnon, $(parsed_args["oesophagusfitnonsensepergene"]))
""";

DF = CSV.File(args["oesophagusdndsneutral"]) |> DataFrame
DF[:A] = 2 * DF[:areacutoff];

DFpatient = filter(row -> row[:name] == "wall", DF)
DFpatient = filter(row -> row[:A] .> 0.1, DFpatient)

x = LLoptimizationresults(DFpatient[:mle], DFpatient[:A]; t = 70.0, Amin = 0.14, ρ = 5000.0)
allmutations = x.DF
@rput allmutations
R"""
write_csv(myDFmiss, $(parsed_args["oesophagusfitneutral"]))
""";


#########################################################
# Skin data
#########################################################

DFdonor = CSV.read(parsed_args["skinmetadata"]);
DFcohort = CSV.read(parsed_args["skindndsdata"]);
DFcohort[:A] = DFcohort[:areacutoff];
patient = unique(DFcohort[:patient])

myDFmiss = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, Float64, Float64],
    [:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :nmutations], 0)

for p in patient
    println(p)
    DFpatient = filter(row -> row[:patient] == p, DFcohort);
    DFpatient = filter(row -> row[:name] == "wmis", DFpatient)
    if DFpatient[:mle][1] == 0.0
        continue
    end
    age = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
    x = LLoptimizationresults(DFpatient[:mle], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
    x.DF[:patient] = p
    x.DF[:Age] = age
    x.DF[:nmutations] = DFpatient[:nmutations][1]
    append!(myDFmiss, x.DF)
end

myDFnon = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, Float64, Float64],
    [:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :nmutations], 0)

for p in patient
    println(p)
    DFpatient = filter(row -> row[:patient] == p, DFcohort);
    DFpatient = filter(row -> row[:name] == "wnon", DFpatient)
    if DFpatient[:mle][1] == 0.0
        continue
    end
    age = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
    x = LLoptimizationresults(DFpatient[:mle], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
    x.DF[:patient] = p
    x.DF[:Age] = age
    x.DF[:nmutations] = DFpatient[:nmutations][1]
    append!(myDFnon, x.DF)
end

myDFall = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, Float64, Float64],
    [:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :nmutations], 0)

for p in patient
    println(p)
    DFpatient = filter(row -> row[:patient] == p, DFcohort);
    DFpatient = filter(row -> row[:name] == "wall", DFpatient)
    if DFpatient[:mle][1] == 0.0
        continue
    end
    age = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
    x = LLoptimizationresults(DFpatient[:mle], DFpatient[:areacutoff]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
    x.DF[:patient] = p
    x.DF[:Age] = age
    x.DF[:nmutations] = DFpatient[:nmutations][1]
    append!(myDFall, x.DF)
end

@rput myDFmiss
@rput myDFnon
R"""
library(readr)
write_csv(myDFmiss, $(parsed_args["skinfitmissense"]))
write_csv(myDFnon, $(parsed_args["skinfitnonsense"]))
""";

DF = CSV.File("data/skin/dnds_patient_genes_combined.csv") |> DataFrame
DF[:A] = DF[:areacutoff];

delta = Float64[]
rlam = Float64[]
genevec = String[]
mutationtypevec = String[]
rsqvec = Float64[]
patientvec = String[]
myDFmiss= DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, Float64, String, Float64],
    [:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :gene, :nmutations], 0)
myDFnon = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, Float64, String, Float64],
    [:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :gene, :nmutations], 0)

for gene in unique(DF[:gene_name])
    println(gene)
    DFgene = filter(row -> row[:gene_name] == gene, DF);

    for p in DFdonor[:patient]
        println(p)
        DFpatient = filter(row -> row[:patient] == p, DFgene);
        age = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
        x = LLoptimizationresults(DFpatient[:wmis_cv], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
        x.DF[:patient] = p
        x.DF[:Age] = age
        x.DF[:gene] = gene
        x.DF[:nmutations] = DFpatient[:n_syn][end] + DFpatient[:n_mis][end] + DFpatient[:n_non][end]
        println(length(x.DF[:A]))
        append!(myDFmiss, x.DF)
        push!(delta, x.Δ)
        push!(rlam, x.rλ)
        push!(genevec, gene)
        push!(mutationtypevec, "Missense")
        push!(rsqvec, x.DF[:rsq][1])
        push!(patientvec, p)
    end

    for p in DFdonor[:patient]
        println(p)
        DFpatient = filter(row -> row[:patient] == p, DFgene);

        age = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
        x = LLoptimizationresults(DFpatient[:wnon_cv], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
        x.DF[:patient] = p
        x.DF[:Age] = age
        x.DF[:gene] = gene
        x.DF[:nmutations] = DFpatient[:n_syn][1] + DFpatient[:n_mis][1] + DFpatient[:n_non][1]
        append!(myDFnon, x.DF)
        push!(delta, x.Δ)
        push!(rlam, x.rλ)
        push!(genevec, gene)
        push!(mutationtypevec, "Nonsense")
        push!(rsqvec, x.DF[:rsq][1])
        push!(patientvec, p)
    end
end

@rput myDFmiss
@rput myDFnon
R"""
write_csv(myDFmiss, $(parsed_args["skinfitmissensepergene"]))
write_csv(myDFnon, $(parsed_args["skinfitnonsensepergene"]))
"""
