rule CalculatedNdSTCGA:
    input:
        tcgadata="data/TCGA-combined-hg38.csv",
        drivergenelist="data/genelists/Driver_gene_list_198_Science_Review.txt",
        essentialgenelist="data/genelists/blomen_essential.csv",
        allgenes="data/genelists/refcdsgenes.txt",
        dndscvref="data/dndscv/RefCDS_human_GRCh38.p12.rda"
    output:
        dndsclonality="results/TCGA/dndsclonality.csv",
        vafclonality="results/TCHA/VAFclonality.csv",
        dndsclonality_percancertype="results/TCGA/dndsclonality_percancertype.csv",
        intervaldnds="results/TCGA/intervaldnds.csv",
        nmutations_gene="results/TCGA/nmutations_gene.csv",
        nmutations_gene_percancertype="results/TCGA/nmutations_gene_percancertype.csv",
        baseline="results/TCGA/baseline.csv",
        baseline_validation="results/TCGA/baseline_validate.csv",
    params:
        step=config["idndslimits"]["step"],
        minarea=config["idndslimits"]["minarea"],
        maxarea=config["idndslimits"]["maxarea"]
    log:
        out = "logs/CalculatedNdS-TCGA.out",
        err = "logs/CalculatedNdS-TCGA.err"
    shell:
        """
        Rscript R/CalculatedNdS-TCGA.R \
            --data {input.tcgadata} \
            --drivergenelist {input.drivergenelist} \
            --essentialgenelist {input.essentialgenelist} \
            --allgenes {input.allgenes} \
            --dndscvref {input.dndscvref} \
            --baseline {output.baseline} \
            --baseline_validation {output.baseline_validation} \
            --vafclonality {output.vafclonality} \
            --dndsclonality {output.dndsclonality} \
            --dndsclonality_percancertype {output.dndsclonality_percancertype} \
            --intervaldnds {output.intervaldnds} \
            --nmutations_gene {output.nmutations_gene} \
            --nmutations_gene_percancertype {output.nmutations_gene_percancertype} #2>> {log.out} 1>> {log.err}
        """
