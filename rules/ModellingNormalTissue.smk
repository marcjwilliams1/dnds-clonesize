rule ModellingNormalTissue:
    output:
        stemcellexamplefit="results/dataforfigures/stemcell_simulation_examplefits.csv",
        stemcellpower="results/dataforfigures/stemcell_simulation_power.csv"
    params:
        nsamples = 5
    shell:
        """
        module load R
        module load julia
        julia julia/ModellingNormalTissue.jl \
        --examplefitsout {output.stemcellexamplefit} \
        --powerout {output.stemcellpower} \
        --nsamples {params.nsamples}
        """
