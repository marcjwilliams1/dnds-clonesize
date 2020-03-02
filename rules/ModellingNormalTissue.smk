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

rule ModellingNormalTissueDistribution:
    output:
        exp="results/dataforfigures/stemcell_simulation_examplefits_distribution-exp.csv",
        beta="results/dataforfigures/stemcell_simulation_examplefits_distribution-beta.csv"
    shell:
        """
        module load R/3.5.3
        module load julia
        export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:`R RHOME`/lib"
        julia julia/ModellingNormalTissue-distribution.jl \
            --resultsfile1 {output.exp} \
            --resultsfile2 {output.beta}
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
        resultsfile="results/simulations/clonesize_overtime.csv",
        resultsfile2="results/simulations/clonesize_overtime-dist.csv"
    shell:
        """
        module load julia
        julia julia/ModellingNormalTissue-clonesize.jl \
            --resultsfile {output.resultsfile} \
            --resultsfile2 {output.resultsfile2}
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
