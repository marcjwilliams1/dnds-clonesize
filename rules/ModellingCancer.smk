rule ModellingCancer:
    output:
        syntheticcohort="results/dataforfigures/syntheticcohort.csv",
        syntheticcohort_diffmu="results/dataforfigures/syntheticcohort_diffmu.csv",
        syntheticcohort_power="results/dataforfigures/syntheticcohort_power.csv",
        syntheticcohort_fmin="results/dataforfigures/syntheticcohort_fmin.csv",
        syntheticcohort_inferreds="results/dataforfigures/syntheticcohort_inferreds.csv"
    shell:
        """
        module load R/3.5.3
        module load julia
        export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:`R RHOME`/lib"
        julia julia/ModellingCancer.jl \
        --syntheticcohort {output.syntheticcohort} \
        --syntheticcohort_diffmu {output.syntheticcohort_diffmu} \
        --syntheticcohort_power {output.syntheticcohort_power} \
        --syntheticcohort_fmin {output.syntheticcohort_fmin} \
        --syntheticcohort_inferreds {output.syntheticcohort_inferreds}
        """
