rule fitdNdSnormal:
    input:
        oesophagusdnds="results/oesophagus/dnds.csv",
        oesophagusdndsgenes="results/oesophagus/dnds_genes.csv",
        oesophagusdndsneutral="results/oesophagus/dnds_neutral.csv",
        oesophagusmetadata="data/oesophagus/donorinfo.csv",
        skindnds="results/skin/dnds.csv",
        skindndsgenes="results/skin/dnds_genes.csv",
        skinmetadata="data/skin/donorinfo.csv",
    output:
        oesophagusfitall = "results/dataforfigures/oesophagusfitall.csv",
        oesophagusfitmissense = "results/dataforfigures/oesophagusfitmissense.csv",
        oesophagusfitnonsense = "results/dataforfigures/oesophagusfitnonsense.csv",
        skinfitmissense = "results/dataforfigures/skinfitmissense.csv",
        skinfitnonsense = "results/dataforfigures/skinfitnonsense.csv",
        oesophagusfitmissensepergene = "results/dataforfigures/oesophagusfitmissensepergene.csv",
        oesophagusfitnonsensepergene = "results/dataforfigures/oesophagusfitnonsensepergene.csv",
        skinfitmissensepergene = "results/dataforfigures/skinfitmissensepergene.csv",
        skinfitnonsensepergene = "results/dataforfigures/skinfitnonsensepergene.csv",
        oesophagusfitneutral = "results/dataforfigures/oesophagusneutral.csv"
    singularity: "docker://marcjwilliams1/julia_r_dnds:v1.1"
    script: "../julia/FitdNdS.jl"

rule fitdNdSnormalSNV:
    input:
        oesophagusdnds="results/oesophagus/dnds_snv.csv",
        oesophagusdndsgenes="results/oesophagus/dnds_genes_snv.csv",
        oesophagusdndsneutral="results/oesophagus/dnds_neutral_snv.csv",
        oesophagusmetadata="data/oesophagus/donorinfo.csv",
        skindnds="results/skin/dnds_snv.csv",
        skindndsgenes="results/skin/dnds_genes_snv.csv",
        skinmetadata="data/skin/donorinfo.csv",
    output:
        oesophagusfitall = "results/dataforfigures/oesophagusfitall_snv.csv",
        oesophagusfitmissense = "results/dataforfigures/oesophagusfitmissense_snv.csv",
        oesophagusfitnonsense = "results/dataforfigures/oesophagusfitnonsense_snv.csv",
        skinfitmissense = "results/dataforfigures/skinfitmissense_snv.csv",
        skinfitnonsense = "results/dataforfigures/skinfitnonsense_snv.csv",
        oesophagusfitmissensepergene = "results/dataforfigures/oesophagusfitmissensepergene_snv.csv",
        oesophagusfitnonsensepergene = "results/dataforfigures/oesophagusfitnonsensepergene_snv.csv",
        skinfitmissensepergene = "results/dataforfigures/skinfitmissensepergene_snv.csv",
        skinfitnonsensepergene = "results/dataforfigures/skinfitnonsensepergene_snv.csv",
        oesophagusfitneutral = "results/dataforfigures/oesophagusneutral_snv.csv"
    singularity: "docker://marcjwilliams1/julia_r_dnds:v1.1"
    script: "../julia/FitdNdS.jl"

rule formatresultsSSB:
    input:
        all = expand("data/SSBresults/allpatients_dnds_{ssb_genes}.txt", ssb_genes=SSB_genes),
        missense = expand("data/SSBresults/allpatients_dnds_{ssb_genes}_mis.txt", ssb_genes=SSB_genes),
        nonsense = expand("data/SSBresults/allpatients_dnds_{ssb_genes}_non.txt", ssb_genes=SSB_genes),
    output:
        all = "results/oesophagus/SSBresults/SSBdnds_results.csv",
        missense = "results/oesophagus/SSBresults/SSBdnds_results_missense.csv",
        nonsense = "results/oesophagus/SSBresults/SSBdnds_results_nonsense.csv"
    params:
        singlepatient=config["patient"],
        step=config["idndslimits"]["step"],
        minarea=config["idndslimits"]["minarea"],
        maxarea=config["idndslimits"]["maxarea"]
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    shell:
        """
        Rscript R/formatSSB.R \
            --inputfile {input.all} \
            --outputfile {output.all} \
            --step {params.step} \
            --minarea {params.minarea} \
            --maxarea {params.maxarea}

        Rscript R/formatSSB.R \
            --inputfile {input.missense} \
            --outputfile {output.missense} \
            --step {params.step} \
            --minarea {params.minarea} \
            --maxarea {params.maxarea}

        Rscript R/formatSSB.R \
            --inputfile {input.nonsense} \
            --outputfile {output.nonsense} \
            --step {params.step} \
            --minarea {params.minarea} \
            --maxarea {params.maxarea}

        """

rule fitdNdSnormalSSB:
    input:
        all = "results/oesophagus/SSBresults/SSBdnds_results.csv",
        missense = "results/oesophagus/SSBresults/SSBdnds_results_missense.csv",
        nonsense = "results/oesophagus/SSBresults/SSBdnds_results_nonsense.csv",
        oesophagusmetadata = "data/oesophagus/donorinfo.csv",
    output: 
        oesophagusfit = "results/dataforfigures/oesophagusfit-SSB.csv",
    singularity: "docker://marcjwilliams1/julia_r_dnds:v1.1"
    script: "../julia/FitdNdS-SSB.jl"
