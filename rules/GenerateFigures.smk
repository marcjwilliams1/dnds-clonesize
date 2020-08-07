rule Figure1:
    output:
        figure="Figures/Figure1.pdf"
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/Figure1.R"

rule Figure2:
    input:
        stemcellexamplefit="results/dataforfigures/stemcell_simulation_examplefits.csv",
        stemcellpower="results/dataforfigures/stemcell_simulation_power.csv",
        oesophagusfitmissense = "results/dataforfigures/oesophagusfitmissense.csv",
        oesophagusfitnonsense = "results/dataforfigures/oesophagusfitnonsense.csv",
        skinfitmissense = "results/dataforfigures/skinfitmissense.csv",
        skinfitnonsense = "results/dataforfigures/skinfitnonsense.csv",
        oesophagusfitmissensepergene = "results/dataforfigures/oesophagusfitmissensepergene.csv",
        oesophagusfitnonsensepergene = "results/dataforfigures/oesophagusfitnonsensepergene.csv",
        singlepatientdnds="results/dataforfigures/singlepatient_bins.csv",
        oesophagusfitneutral = "results/dataforfigures/oesophagusneutral.csv"
    output:
        figure="Figures/Figure2.pdf",
        suppfigures=expand("Figures/Figure1-S{S}.pdf", S = [1])
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/Figure2.R"

rule Figure3:
    input:
        stemcellexamplefit="results/dataforfigures/stemcell_simulation_examplefits.csv",
        stemcellpower="results/dataforfigures/stemcell_simulation_power.csv",
        oesophagusfitmissense = "results/dataforfigures/oesophagusfitmissense.csv",
        oesophagusfitnonsense = "results/dataforfigures/oesophagusfitnonsense.csv",
        oesophagusfitmissensepergene = "results/dataforfigures/oesophagusfitmissensepergene.csv",
        oesophagusfitnonsensepergene = "results/dataforfigures/oesophagusfitnonsensepergene.csv",
    params:
        mutationcutoff=config["mutationcutoff"],
        rsqcutoff=config["rsqcutoff"]
    output:
        figure="Figures/Figure3.pdf",
        suppfigures=["Figures/Figure2-S2.pdf", "Figures/Figure2-S3.pdf", "Figures/Figure2-S4.pdf", "Figures/Figure3-S1.pdf"]
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/Figure3.R"

rule Figure4:
    input:
        skinfitmissensepergene = "results/dataforfigures/skinfitmissensepergene.csv",
        skinfitnonsensepergene = "results/dataforfigures/skinfitnonsensepergene.csv",
        skinfitmissense = "results/dataforfigures/skinfitmissense.csv",
        skinfitnonsense = "results/dataforfigures/skinfitnonsense.csv",
    output:
        figure="Figures/Figure4.pdf",
        suppfigures=expand("Figures/Figure4-S{S}.pdf", S = [1])
    params:
        mutationcutoff=config["mutationcutoff"],
        rsqcutoff=config["rsqcutoff"]
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/FigureS6.R"

rule Figurecomparednds:
    input:
        alldndscv = "results/dataforfigures/oesophagusfitall_snv.csv",
        dndscvfitmissense = "results/dataforfigures/oesophagusfitmissense_snv.csv",
        dndscvfitnonsense = "results/dataforfigures/oesophagusfitnonsense_snv.csv",
        dndscvfitmissensepergene = "results/dataforfigures/oesophagusfitmissensepergene_snv.csv",
        dndscvfitnonsensepergene = "results/dataforfigures/oesophagusfitnonsensepergene_snv.csv",
        SSB = "results/dataforfigures/oesophagusfit-SSB.csv"
    output:
        suppfigures=["Figures/Figure5-S5.pdf","Figures/Extra/Figure-Extra-SSBfits.pdf"]
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/Figure-comparednds.R"

rule Figurebinsize:
    input:
        binsizesims = "results/dataforfigures/stemcell_simulation_differentbins.csv"
    output:
        suppfigures=["Figures/Extra/Figure-Extra-binsize.pdf"]
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/Figure-binsize.R"

rule Figure5:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
        fits="results/dataforfigures/data-clonesizefit.Rdata",
        modelfits="results/dataforfigures/data-clonesizefit-models.Rdata"
    output:
        figure="Figures/Figure5.pdf",
        suppfigures=expand("Figures/Figure5-S{S}.pdf", S = [3]),
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/Figure-clonesizedata.R"

rule FigureCloneSizeData2:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        fits="results/dataforfigures/brmsfit.Rdata",
        oesophagusfitmissense = "results/dataforfigures/oesophagusfitmissensepergene.csv",
    output:
        suppfigures=["Figures/Figure5-S4.pdf","Figures/Extra/Figure-Extra-PPcheckregression.pdf"]
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/Figure-clonesizedata2.R"

rule FigureCloneSizeSims:
    input:
        data="results/simulations/clonesize_overtime.csv",
        fits="results/dataforfigures/simulation-clonesizefit.Rdata",
    output:
        suppfigures=["Figures/Figure5-S2.pdf","Figures/Extra/Figure-Extra-PPcheck.pdf"]
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/Figure-clonesizesims.R"

rule FiguredNdSsimsdist:
    input:
        exp="results/dataforfigures/stemcell_simulation_examplefits_distribution-exp.csv",
        beta="results/dataforfigures/stemcell_simulation_examplefits_distribution-beta.csv"
    output:
        "Figures/Figure2-S1.pdf"
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/FiguredNdSsimsdist.R"

rule FigureCloneSizeSimsDist:
    input:
        data="results/simulations/clonesize_overtime-dist.csv",
        fits="results/dataforfigures/simulation-clonesizefit-dist.Rdata",
    output:
        suppfigures="Figures/Extra/Figure-Extra-SimsDistributionRegression.pdf"
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/Figure-clonesizesimsdist.R"

rule FigureHitchikers:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
        oesophagusdata_all="data/oesophagus/aau3879_TableS2.xlsx",
        simulationdata="results/simulations/clonesize_hitchikers.csv"
    output:
        suppfigures=["Figures/Figure5-S1.pdf"]
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    script: "../R/Figure-hitchikers.R"
