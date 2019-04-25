rule fitdNdSnormal:
    input:
        oesophagusdnds="results/oesophagus/dnds.csv",
        oesophagusdndsgenes="results/oesophagus/dnds_genes.csv",
        skindnds="results/skin/dnds.csv",
        skindndsgenes="results/skin/dnds_genes.csv",
        oesophagusdndsneutral="results/esophagus/dnds_neutral.csv",
        oesophagusdndsgenesneutral="results/oesophagus/dnds_genes_neutral.csv",
        singlepatientdnds="results/dataforfigures/{input.singlepatient}_bins.csv",
        singlepatientdndsgenes="results/dataforfigures/{input.singlepatient}_bins_genes.csv"
    output:
    params:
    log:
        out = "logs/fitdNdS-normal.out",
        err = "logs/fitdNdS-normal.err"
    shell:
    """
    julia julia/fitdNdS.jl \
     2>> {log.out} 1>> {log.err}
    """
