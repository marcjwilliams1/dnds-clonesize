rule ModellingNormalTissue:
    output:
        stemcellexamplefit="results/dataforfigures/stemcell_simulation_examplefits.csv",
        stemcellexampledifferentbins="results/dataforfigures/stemcell_simulation_differentbins.csv",
    params:
        nsamples = config["nsamplesnormal"]
    shell:
        """
        module load R/3.5.3
        module load julia
        export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:`R RHOME`/lib"
        julia julia/ModellingNormalTissue.jl \
            --examplefitsout {output.stemcellexamplefit} \
            --exampledifferentbins {output.stemcellexampledifferentbins} \
            --nsamples {params.nsamples}
        """

rule ModellingNormalTissuePower:
    output:
        stemcellpower="results/dataforfigures/stemcell_simulation_power.csv"
    params:
        nsamples = config["nsamplesnormal"]
    shell:
        """
        module load R/3.5.3
        module load julia
        export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:`R RHOME`/lib"
        julia julia/ModellingNormalTissuePower.jl \
            --powerout {output.stemcellpower} \
            --nsamples {params.nsamples}
        """

rule ModellingNormalTissueCloneSize:
    output:
        resultsfile="results/simulations/clonesize_overtime.csv"
    shell:
        """
        module load julia
        julia julia/ModellingNormalTissue-clonesize.jl \
            --resultsfile {output.resultsfile}
        """

rule ModellingNormalTissueHitchikers:
    output:
        resultsfile="results/simulations/clonesize_hitchikers.csv"
    shell:
        """
        module load julia
        julia julia/ModellingNormalTissue-hitchikers.jl \
            --resultsfile {output.resultsfile}
        """
