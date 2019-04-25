rule CalculatedNdSNormal:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
        skindata="data/skin/skin_mutation.csv"
    output:
        oesophagusdnds="results/oesophagus/dnds.csv",
        oesophagusdndsgenes="results/oesophagus/dnds_genes.csv",
        skindnds="results/skin/dnds.csv",
        skindndsgenes="results/skin/dnds_genes.csv",
        oesophagusdndsneutral="results/esophagus/dnds_neutral.csv",
        oesophagusdndsgenesneutral="results/oesophagus/dnds_genes_neutral.csv",
        singlepatientdnds="results/dataforfigures/{input.singlepatient}_bins.csv",
        singlepatientdndsgenes="results/dataforfigures/{input.singlepatient}_bins_genes.csv"
    params:
        singlepatient=config["patient"],
        step=config["idndslimits"]["step"],
        minarea=config["idndslimits"]["minarea"],
        maxarea=config["idndslimits"]["maxarea"]
    log:
        out = "logs/Rule2-CalculatedNdS-normal.out",
        err = "logs/Rule2-CalculatedNdS-normal.err"
    shell:
        """
        Rscript R/Rule2-CalculatedNdS-normal.R \
        --oesophaguspatientinfo {input.oesophaguspatientinfo} \
        --oesophagusdata {input.oesophagusdata} \
        --skindata {input.skindata} \
        --oesophagusdnds {output.oesophagusdnds} \
        --oesophagusdndsgenes {output.oesophagusdndsgenes} \
        --skindnds {output.skindnds} \
        --skindndsgenes {output.skindndsgenes} \
        --oesophagusdndsneutral {output.oesophagusdndsneutral} \
        --oesophagusdndsgenesneutral {output.oesophagusdndsgenesneutral} \
        --singlepatient {params.singlepatient} \
        --singlepatientdnds {output.singlepatientdnds} \
        --singlepatientdndsgenes {output.singlepatientdndsgenes} \
        --step {params.step} \
        --minarea {params.minarea} \
        --maxarea {params.maxarea}  2>> {log.out} 1>> {log.err}
        """
