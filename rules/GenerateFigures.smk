rule Figure1:
    log:
        out = "logs/fitdNdS-normal.out",
        err = "logs/fitdNdS-normal.err"
    output:
        figure="Figures/Figure1.pdf"
    shell:
        """
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
        #skinfitmissensepergene = "results/dataforfigures/skinfitmissensepergene.csv",
        #skinfitnonsensepergene = "results/dataforfigures/skinfitnonsensepergene.csv",
        oesophagusfitneutral = "results/dataforfigures/oesophagusneutral.csv"
    log:
        out = "logs/fitdNdS-normal.out",
        err = "logs/fitdNdS-normal.err"
    output:
        figure="Figures/Figure2.pdf"
    shell:
        """
        Rscript R/Figure2.R \
            --figure {output.figure} \
            --oesophagusfitmissense {input.oesophagusfitmissense} \
            --oesophagusfitnonsense {input.oesophagusfitnonsense} \
            --skinfitmissense {input.skinfitmissense} \
            --skinfitnonsense {input.skinfitnonsense} \
            --oesophagusfitmissensepergene {input.oesophagusfitmissensepergene} \
            --oesophagusfitnonsensepergene {input.oesophagusfitnonsensepergene} \
            --oesophagusfitneutral {input.oesophagusfitneutral} \
            --stemcellexamplefit {input.stemcellexamplefit} \
            --stemcellpower {input.stemcellpower}
        """

rule Figure3:
    input:
        stemcellexamplefit="results/dataforfigures/stemcell_simulation_examplefits.csv",
        stemcellpower="results/dataforfigures/stemcell_simulation_power.csv",
        oesophagusfitmissense = "results/dataforfigures/oesophagusfitmissense.csv",
        oesophagusfitnonsense = "results/dataforfigures/oesophagusfitnonsense.csv",
        skinfitmissense = "results/dataforfigures/skinfitmissense.csv",
        skinfitnonsense = "results/dataforfigures/skinfitnonsense.csv",
        oesophagusfitmissensepergene = "results/dataforfigures/oesophagusfitmissensepergene.csv",
        oesophagusfitnonsensepergene = "results/dataforfigures/oesophagusfitnonsensepergene.csv",
    log:
        out = "logs/fitdNdS-normal.out",
        err = "logs/fitdNdS-normal.err"
    output:
        figure="Figures/Figure3.pdf",
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [1,2,3,4])
    shell:
        """
        Rscript R/Figure3.R \
            --figure {output.figure} \
            --suppfigures {output.suppfigures} \
            --oesophagusfitmissense {input.oesophagusfitmissense} \
            --oesophagusfitnonsense {input.oesophagusfitnonsense} \
            --oesophagusfitmissensepergene {input.oesophagusfitmissensepergene} \
            --oesophagusfitnonsensepergene {input.oesophagusfitnonsensepergene}
        """

rule Figure4:
    input:
        tcgadata="data/TCGA-combined-hg38.csv",
        dndsclonality_percancertype="results/TCGA/dndsclonality_percancertype.csv",
        dndsclonality="results/TCGA/dndsclonality.csv",
        baseline="results/TCGA/baseline.csv",
    output:
        figure="Figures/Figure4.pdf",
        suppfigures=expand("Figures/FigureS{S}.pdf", S = [5])
    shell:
        """
        Rscript R/Figure4.R \
            --figure {output.figure} \
            --suppfigures {output.suppfigures} \
            --tcgadata {input.tcgadata} \
            --baseline {input.baseline} \
            --dndsclonality_percancertype {input.dndsclonality_percancertype} \
            --dndsclonality {input.dndsclonality}
        """
