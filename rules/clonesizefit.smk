rule brms:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
    output:
        "results/dataforfigures/brmsfit.Rdata"
    threads: 4
    shell:
        """
        module load R/3.5.3
        module load gcc
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
        "results/dataforfigures/brmsfit-sites.Rdata"
    threads: 4
    params:
        singularityimage=config["stansingularity"]
    shell:
        """
        #module load gcc
        #module load R/3.5.3
        module unload python
        module load singularity
        singularity exec {params.singularityimage} \
            Rscript R/fitbrms-sites.R \
            --oesophagusdata {input.oesophagusdata} \
            --oesophagusmetadata {input.oesophaguspatientinfo} \
            --output {output} \
            --threads {threads}
        """

rule clonesizesims:
    input:
        data="results/simulations/clonesize_overtime.csv"
    output:
        "results/dataforfigures/simulation-clonesizefit.Rdata"
    threads: 4
    params:
        singularityimage=config["stansingularity"]
    shell:
        """
        module unload python
        module load singularity
        singularity exec {params.singularityimage} \
            Rscript R/fitclonesize-sims.R \
            --simulationdata {input.data} \
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
    params:
        singularityimage=config["stansingularity"]
    shell:
        """
        module unload python
        module load singularity
        singularity exec {params.singularityimage} \
            Rscript R/fitclonesize-data.R \
            --oesophagusdata {input.oesophagusdata} \
            --oesophagusmetadata {input.oesophaguspatientinfo} \
            --output {output} \
            --threads {threads}
        """
