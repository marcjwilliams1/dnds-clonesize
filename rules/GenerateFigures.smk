rule Figure1:
    output:
        figure="Figures/Figure1.pdf"
    shell:
        """
        module load R
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
        module load R
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
        module load R
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


rule FigureS6:
    input:
        skinfitmissensepergene = "results/dataforfigures/skinfitmissensepergene.csv",
        skinfitnonsensepergene = "results/dataforfigures/skinfitnonsensepergene.csv",
        skinfitmissense = "results/dataforfigures/skinfitmissense.csv",
        skinfitnonsense = "results/dataforfigures/skinfitnonsense.csv",
    output:
        figure="Figures/FigureS6.pdf",
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [10])
    params:
        mutationcutoff=config["mutationcutoff"],
        rsqcutoff=config["rsqcutoff"]
    shell:
        """
        module load R
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
        alldndscv = "results/dataforfigures/oesophagusfitall.csv",
        dndscvfitmissense = "results/dataforfigures/oesophagusfitmissense.csv",
        dndscvfitnonsense = "results/dataforfigures/oesophagusfitnonsense.csv",
        dndscvfitmissensepergene = "results/dataforfigures/oesophagusfitmissensepergene.csv",
        dndscvfitnonsensepergene = "results/dataforfigures/oesophagusfitnonsensepergene.csv",
        SSB = "results/dataforfigures/oesophagusfit-SSB.csv"
    output:
        figure = "Figures/FigureS13.pdf"
    shell:
        """
        module load R
        Rscript R/Figure-comparednds.R \
            --figure {output.figure} \
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
        figure = "Figures/FigureS12.pdf"
    shell:
        """
        module load R
        Rscript R/Figure-binsize.R \
            --figure {output.figure} \
            --binsizesims {input.binsizesims}
        """
