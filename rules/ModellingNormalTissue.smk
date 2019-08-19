rule ModellingNormalTissue:
    output:
        stemcellexamplefit="results/dataforfigures/stemcell_simulation_examplefits.csv",
        stemcellpower="results/dataforfigures/stemcell_simulation_power.csv"
    params:
        nsamples = 5
    log:
        out = "logs/Rule1-ModellingNormalTissue.out",
        err = "logs/Rule1-ModellingNormalTissue.err"
    shell:
        """
        julia julia/ModellingNormalTissue.jl \
        --examplefitsout {output.stemcellexamplefit} \
        --powerout {output.stemcellpower} \
        --nsamples {params.nsamples} 2>> {log.out} 1>> {log.err}
        """
