rule brms:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
    output:
        "results/dataforfigures/brmsfit.Rdata"
    threads: 4
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    shell:
        """
        Rscript R/fitbrms.R \
            --oesophagusdata {input.oesophagusdata} \
            --oesophagusmetadata {input.oesophaguspatientinfo} \
            --output {output}
        """

rule brmssites:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
    output:
        fits = "results/dataforfigures/brmsfit-sites.Rdata",
        coefficients = "results/dataforfigures/brmsfit-sites-coef.Rdata"
    threads: 4
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    shell:
        """
        Rscript R/fitbrms-sites.R \
            --oesophagusdata {input.oesophagusdata} \
            --oesophagusmetadata {input.oesophaguspatientinfo} \
            --output {output.fits} \
            --outputcoef {output.coefficients} \
            --threads {threads}
        """

rule clonesizesims:
    input:
        data="results/simulations/clonesize_overtime.csv",
        datahitchike="results/simulations/clonesize_hitchikers.csv"
    output:
        "results/dataforfigures/simulation-clonesizefit.Rdata"
    threads: 4
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    shell:
        """
        Rscript R/fitclonesize-sims.R \
            --simulationdata {input.data} \
            --simulationdatahitchike {input.datahitchike} \
            --output {output} \
            --threads {threads}
        """

rule clonesizesimsdist:
    input:
        data="results/simulations/clonesize_overtime-dist.csv",
        datahitchike="results/simulations/clonesize_hitchikers.csv"
    output:
        "results/dataforfigures/simulation-clonesizefit-dist.Rdata"
    threads: 4
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    shell:
        """
        Rscript R/fitclonesize-sims.R \
            --simulationdata {input.data} \
            --simulationdatahitchike {input.datahitchike} \
            --output {output} \
            --threads {threads}
        """

rule clonesizedata:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
    output:
        "results/dataforfigures/data-clonesizefit.Rdata"
    threads: 4
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    shell:
        """
        Rscript R/fitclonesize-data.R \
            --oesophagusdata {input.oesophagusdata} \
            --oesophagusmetadata {input.oesophaguspatientinfo} \
            --output {output} \
            --threads {threads}
        """

rule clonesizedatacompare:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
    output:
        "results/dataforfigures/data-clonesizefit-models.Rdata"
    threads: 4
    singularity: "shub://marcjwilliams1/dnds-clonesize-R-container"
    shell:
        """
        Rscript R/fitclonesize-data-compare.R \
            --oesophagusdata {input.oesophagusdata} \
            --oesophagusmetadata {input.oesophaguspatientinfo} \
            --output {output} \
            --threads {threads}
        """
