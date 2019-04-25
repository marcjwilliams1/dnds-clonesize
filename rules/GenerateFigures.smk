rule GenerateFigures:
    input:
        stemcellexamplefit="results/dataforfigures/stemcell_simulation_examplefits.csv",
        stemcellpower="results/dataforfigures/stemcell_simulation_power.csv"
    output:
        fig1 = "Figures/Figure1.pdf",
        fig2 = "Figures/Figure2.pdf"
    shell:
        """
        touch {output.fig1}
        touch {output.fig2}
        """
