rule Figure1:
    output:
        figure="Figures/Figure1.pdf"
    shell:
        """
        module load gcc
        module load R/3.5.3
        Rscript R/Figure1.R \
            --figure {output.figure}
        """

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
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [1])
    shell:
        """
        module load gcc
        module load R/3.5.3
        Rscript R/Figure2.R \
            --figure {output.figure} \
            --suppfigures {output.suppfigures} \
            --oesophagusfitmissense {input.oesophagusfitmissense} \
            --oesophagusfitnonsense {input.oesophagusfitnonsense} \
            --skinfitmissense {input.skinfitmissense} \
            --skinfitnonsense {input.skinfitnonsense} \
            --oesophagusfitmissensepergene {input.oesophagusfitmissensepergene} \
            --oesophagusfitnonsensepergene {input.oesophagusfitnonsensepergene} \
            --oesophagusfitneutral {input.oesophagusfitneutral} \
            --stemcellexamplefit {input.stemcellexamplefit} \
            --stemcellpower {input.stemcellpower} \
            --singlepatient {input.singlepatientdnds}
        """

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
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [2,3,4,5])
    shell:
        """
        module load gcc
        module load R/3.5.3
        Rscript R/Figure3.R \
            --figure {output.figure} \
            --suppfigures {output.suppfigures} \
            --oesophagusfitmissense {input.oesophagusfitmissense} \
            --oesophagusfitnonsense {input.oesophagusfitnonsense} \
            --oesophagusfitmissensepergene {input.oesophagusfitmissensepergene} \
            --oesophagusfitnonsensepergene {input.oesophagusfitnonsensepergene} \
            --mutationcutoff {params.mutationcutoff} \
            --rsqcutoff {params.rsqcutoff}
         """

rule Figure4:
    input:
        skinfitmissensepergene = "results/dataforfigures/skinfitmissensepergene.csv",
        skinfitnonsensepergene = "results/dataforfigures/skinfitnonsensepergene.csv",
        skinfitmissense = "results/dataforfigures/skinfitmissense.csv",
        skinfitnonsense = "results/dataforfigures/skinfitnonsense.csv",
    output:
        figure="Figures/Figure4.pdf",
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [10])
    params:
        mutationcutoff=config["mutationcutoff"],
        rsqcutoff=config["rsqcutoff"]
    shell:
        """
        module load gcc
        module load R/3.5.3
        Rscript R/FigureS6.R \
            --figure {output.figure} \
            --suppfigures {output.suppfigures} \
            --skinfitmissense {input.skinfitmissense} \
            --skinfitnonsense {input.skinfitnonsense} \
            --skinfitmissensepergene {input.skinfitmissensepergene} \
            --skinfitnonsensepergene {input.skinfitnonsensepergene} \
            --mutationcutoff {params.mutationcutoff} \
            --rsqcutoff {params.rsqcutoff}
        """

rule Figurecomparednds:
    input:
        alldndscv = "results/dataforfigures/oesophagusfitall_snv.csv",
        dndscvfitmissense = "results/dataforfigures/oesophagusfitmissense_snv.csv",
        dndscvfitnonsense = "results/dataforfigures/oesophagusfitnonsense_snv.csv",
        dndscvfitmissensepergene = "results/dataforfigures/oesophagusfitmissensepergene_snv.csv",
        dndscvfitnonsensepergene = "results/dataforfigures/oesophagusfitnonsensepergene_snv.csv",
        SSB = "results/dataforfigures/oesophagusfit-SSB.csv"
    output:
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [12,13])
    shell:
        """
        module load gcc
        module load R/3.5.3
        Rscript R/Figure-comparednds.R \
            --suppfigures {output.suppfigures} \
            --SSB {input.SSB} \
            --alldndscv {input.alldndscv} \
            --dndscvmissense {input.dndscvfitmissense} \
            --dndscvnonsense {input.dndscvfitnonsense} \
            --dndscvmissensepergene {input.dndscvfitmissensepergene} \
            --dndscvnonsensepergene {input.dndscvfitnonsensepergene} \
        """

rule Figurebinsize:
    input:
        binsizesims = "results/dataforfigures/stemcell_simulation_differentbins.csv"
    output:
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [14])
    shell:
        """
        module load gcc
        module load R/3.5.3
        Rscript R/Figure-binsize.R \
            --suppfigures {output.suppfigures} \
            --binsizesims {input.binsizesims}
        """

rule Figure5:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
        fits="results/dataforfigures/data-clonesizefit.Rdata",
        modelfits="results/dataforfigures/data-clonesizefit-models.Rdata"
    output:
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [16]),
        figure="Figures/Figure5.pdf"
    params:
        singularityimage=config["stansingularity"]
    shell:
        """
        module unload python
        module load singularity
        singularity exec {params.singularityimage} \
        Rscript R/Figure-clonesizedata.R \
            --datafits {input.fits} \
            --datamodelfits {input.modelfits} \
            --figure {output.figure} \
            --oesophagusdata {input.oesophagusdata} \
            --oesophagusmetadata {input.oesophaguspatientinfo} \
            --suppfigures {output.suppfigures}

        """

rule FigureCloneSizeData2:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        fits="results/dataforfigures/brmsfit.Rdata",
        oesophagusfitmissense = "results/dataforfigures/oesophagusfitmissensepergene.csv",
    output:
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [17, 18]),
    params:
        singularityimage=config["stansingularity"]
    shell:
        """
        module unload python
        module load singularity
        singularity exec {params.singularityimage} \
        Rscript R/Figure-clonesizedata2.R \
            --datafits {input.fits} \
            --oesophagusfitmissense {input.oesophagusfitmissense} \
            --oesophagusmetadata {input.oesophaguspatientinfo} \
            --suppfigures {output.suppfigures}

        """

rule FigureCloneSizeSims:
    input:
        data="results/simulations/clonesize_overtime.csv",
        fits="results/dataforfigures/simulation-clonesizefit.Rdata",
    output:
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [15, 25]),
    params:
        singularityimage=config["stansingularity"]
    shell:
        """
        module unload python
        module load singularity
        singularity exec {params.singularityimage} \
        Rscript R/Figure-clonesizesims.R \
            --simulationfits {input.fits} \
            --simulationdata {input.data} \
            --suppfigures {output.suppfigures}

        """

rule FigureHitchikers:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
        oesophagusdata_all="data/oesophagus/aau3879_TableS2.xlsx",
        simulationdata="results/simulations/clonesize_hitchikers.csv"
    output:
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [26]),
    params:
        singularityimage=config["stansingularity"]
    shell:
        """
        module unload python
        module load singularity
        singularity exec {params.singularityimage} \
        Rscript R/Figure-hitchikers.R \
            --simulationdata {input.simulationdata} \
            --oesophagusdata {input.oesophagusdata} \
            --oesophagusdata_all {input.oesophagusdata_all} \
            --oesophagusmetadata {input.oesophaguspatientinfo} \
            --suppfigures {output.suppfigures}
        """
