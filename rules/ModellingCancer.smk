rule ModellingCancer:
    output:
        syntheticcohort="results/dataforfigures/syntheticcohort.csv,
        syntheticcohort_diffmu="results/dataforfigures/syntheticcohort_diffmu.csv,
        syntheticcohort_power="results/dataforfigures/syntheticcohort_power.csv,
        syntheticcohort_fmin="results/dataforfigures/syntheticcohort_fmin.csv,
    log:
        out = "logs/ModellingCancer.out",
        err = "logs/ModellingCancer.err"
    shell:
        """
        julia julia/ModellingCancer.jl \
        --syntheticcohort {output.syntheticcohort} \
        --syntheticcohort_diffmu {output.syntheticcohort_diffmu} \
        --syntheticcohort_power {output.syntheticcohort_power} \
        --syntheticcohort_fmin {output.syntheticcohort_fmin} 2>> {log.out} 1>> {log.err}
        """
