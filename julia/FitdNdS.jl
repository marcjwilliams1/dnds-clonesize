using CSV
using DataFrames
using RCall

@rimport readr as readr

# load functions to fit data via least squares
println(snakemake)
include(string(snakemake.scriptdir, "/optim.jl"))

println("Read in data...")
println(snakemake.input["oesophagusmetadata"])
DFdonor = rcopy(readr.read_csv(snakemake.input["oesophagusmetadata"]))
DFcohort = rcopy(readr.read_csv(snakemake.input["oesophagusdnds"]))
DFcohort[:A] = 2 * DFcohort[:areacutoff];

println("Reading in per gene data...")
DF = rcopy(readr.read_csv(snakemake.input["oesophagusdndsgenes"]))
DF[:A] = 2 .* DF[:areacutoff];
println("Data read in")

myDFmiss = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64, Float64,
    Float64, Float64, String, String, Float64, Float64],
[:dnds, :A, :dndsfit, :dndsfitlq, :dndsfituq, :deltafit, :lambdarfit,
    :deltafitlq, :lambdarfitlq, :deltafituq, :lambdarfituq, :sedelta, :selambda,
    :rsq, :patient, :Age, :Age2, :nmutations], 0)

println("Analyse per patient data for missense mutations...")
for p in DFdonor[:patient]
    println("Analysing patient $p")
    DFpatient = filter(row -> row[:patient] == p, DFcohort);
    DFpatient = filter(row -> row[:name] == "wmis", DFpatient)
    if DFpatient[:mle][1] == 0.0
        continue
    end
    age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]
    x = LLoptimizationresults(DFpatient[:mle], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
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

println("Analyse per patient data for nonsense mutations...")
for p in DFdonor[:patient]
    println("Analysing patient $p")
    DFpatient = filter(row -> row[:patient] == p, DFcohort);
    DFpatient = filter(row -> row[:name] == "wnon", DFpatient)
    if DFpatient[:mle][1] == 0.0
        continue
    end
    age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]
    x = LLoptimizationresults(DFpatient[:mle], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
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

println("Analyse per patient data for all mutations...")
for p in DFdonor[:patient]
    println("Analysing patient $p")
    DFpatient = filter(row -> row[:patient] == p, DFcohort);
    DFpatient = filter(row -> row[:name] == "wall", DFpatient)
    if DFpatient[:mle][1] == 0.0
        continue
    end
    age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]
    x = LLoptimizationresults(DFpatient[:mle], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
    x.DF[:patient] = p
    x.DF[:Age] = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
    x.DF[:Age2] = age
    x.DF[:nmutations] = DFpatient[:nmutations]
    append!(myDFall, x.DF)
end

println("Writing data to file....")
@rput myDFmiss
@rput myDFnon
@rput myDFall
R"""
library(readr)
write_csv(myDFmiss, $(snakemake.output["oesophagusfitmissense"]))
write_csv(myDFnon, $(snakemake.output["oesophagusfitnonsense"]))
write_csv(myDFall, $(snakemake.output["oesophagusfitall"]))
""";

println("Reading in per gene data...")
DF = rcopy(readr.read_csv(snakemake.input["oesophagusdndsgenes"]))
DF[:A] = 2 .* DF[:areacutoff];


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

println("Analyse per patient per gene data for missense mutations...")
for gene in unique(DF[:gene_name])
#for gene in ["TP53", "NOTCH1"]
    println("Analysing gene: $gene")
    DFgene = filter(row -> row[:gene_name] == gene, DF);
    println("Analysing missense mutations")
    if gene == "SPHKAP"
        println("Errors with this gene, moving to next patient")
        continue
    end
    for p in DFdonor[:patient]
        DFpatient = filter(row -> row[:patient] == p, DFgene);
        age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]

        nmuts = DFpatient[:n_syn][end] + DFpatient[:n_mis][end]
        if nmuts < 5
            println("Patient $p, only $(nmuts) mutations, moving to next patient")
            continue
        end
        println("Analysing patient $p, $(nmuts) mutations")

        try 
            global x = LLoptimizationresults(DFpatient[:wmis_cv], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
        catch e
            println("Patient $p failed, moving to next patient")
            continue
        end
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

    println("Analysing nonsense mutations")
    for p in DFdonor[:patient]
        println("Analysing patient $p")
        DFpatient = filter(row -> row[:patient] == p, DFgene);
        age = filter(row -> row[:patient] == p, DFdonor)[:Age2][1]

        nmuts = DFpatient[:n_syn][end] + DFpatient[:n_non][end]
        if nmuts < 5
            println("Only $(nmuts) mutations, moving to next patient")
            continue
        end

        println("Analysing patient $p, $(nmuts) mutations")

        try
            global x = LLoptimizationresults(DFpatient[:wnon_cv], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
        catch e
            println("Patient $p failed, moving to next patient")
            continue
        end
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

println("Number of fits missense: $(length(unique(myDFmiss[:gene])))")
println("Number of fits nonsense: $(length(unique(myDFnon[:gene])))")

@rput myDFmiss
@rput myDFnon
R"""
library(readr)
write_csv(myDFmiss, $(snakemake.output["oesophagusfitmissensepergene"]))
write_csv(myDFnon, $(snakemake.output["oesophagusfitnonsensepergene"]))
""";

#DF = CSV.File(args["oesophagusdndsneutral"]) |> DataFrame
DF = rcopy(readr.read_csv((snakemake.input["oesophagusdndsneutral"])))
DF[:A] = 2 * DF[:areacutoff];

DFpatient = filter(row -> row[:name] == "wall", DF)
DFpatient = filter(row -> row[:A] .> 0.1, DFpatient)

x = LLoptimizationresults(DFpatient[:mle], DFpatient[:A]; t = 70.0, Amin = 0.14, ρ = 5000.0)
allmutations = x.DF
@rput allmutations
R"""
write_csv(allmutations, $(snakemake.output["oesophagusfitneutral"]))
""";


#########################################################
# Skin data
#########################################################
println("Analysing Skin data")
DFdonor = rcopy(readr.read_csv(snakemake.input["skinmetadata"]));
DFcohort = rcopy(readr.read_csv(snakemake.input["skindnds"]));
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
write_csv(myDFmiss, $(snakemake.output["skinfitmissense"]))
write_csv(myDFnon, $(snakemake.output["skinfitnonsense"]))
""";


#DF = CSV.File("data/skin/dnds_patient_genes_combined.csv") |> DataFrame
DF = rcopy(readr.read_csv(snakemake.input["skindndsgenes"], guess_max = 10^5))
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
    if gene == "DICER1"
        continue
        println("Skipping DICER1")
    end

    if gene == "KCNH5"
        continue
        println("Skipping KCNH5")
    end

    for p in DFdonor[:patient]
        println(p)
        DFpatient = filter(row -> row[:patient] == p, DFgene);
        age = filter(row -> row[:patient] == p, DFdonor)[:Age][1]

        nmuts = DFpatient[:n_syn][end] + DFpatient[:n_mis][end]
        if nmuts < 5
            println("Patient $p, only $(nmuts) mutations, moving to next patient")
            continue
        end
        println("Analysing patient $p, $(nmuts) mutations")

        try 
            global x = LLoptimizationresults(DFpatient[:wmis_cv], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
        catch e
            println("Patient $p failed, moving to next patient")
            continue
        end
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

        nmuts = DFpatient[:n_syn][end] + DFpatient[:n_mis][end]
        if nmuts < 5
            println("Patient $p, only $(nmuts) mutations, moving to next patient")
            continue
        end
        println("Analysing patient $p, $(nmuts) mutations")

        age = filter(row -> row[:patient] == p, DFdonor)[:Age][1]
        try 
            global x = LLoptimizationresults(DFpatient[:wnon_cv], DFpatient[:A]; t = age, Amin = 2 * 0.01, ρ = 5000.0)
        catch e
            println("Patient $p failed, moving to next patient")
            continue
        end
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
write_csv(myDFmiss, $(snakemake.output["skinfitmissensepergene"]))
write_csv(myDFnon, $(snakemake.output["skinfitnonsensepergene"]))
"""
