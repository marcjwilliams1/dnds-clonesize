rule ModellingNormalTissue:
    output:
        stemcellexamplefit="results/dataforfigures/stemcell_simulation_examplefits.csv",
        stemcellpower="results/dataforfigures/stemcell_simulation_power.csv"
    params:
        nsamples = config["nsamplesnormal"]
    shell:
        """
        module load R
        module load julia
        export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:`R RHOME`/lib"
        julia julia/ModellingNormalTissue.jl \
        --examplefitsout {output.stemcellexamplefit} \
        --powerout {output.stemcellpower} \
        --nsamples {params.nsamples}
        """
