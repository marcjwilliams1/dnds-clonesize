rule CalculatedNdSTCGA:
    input:
        tcgadata="data/TCGA-combined-hg38.csv",
        drivergenelist="data/genelists/Driver_gene_list_198_Science_Review.txt",
        essentialgenelist="data/genelists/blomen_essential.csv",
        allgenes="data/genelists/refcdsgenes.txt",
        dndscvref="data/dndscv/RefCDS_human_GRCh38.p12.rda"
    output:
        dndsclonality="results/TCGA/dndsclonality.csv",
        dndsclonality_percancertype="results/TCGA/dndsclonality_percancertype.csv",
        intervaldnds="results/TCGA/intervaldnds.csv",
        numberofgenes="results/TCGA/numberofgenes.csv",
        skindnds="results/skin/dnds.csv",
        skindndsgenes="results/skin/dnds_genes.csv",
        oesophagusdndsneutral="results/oesophagus/dnds_neutral.csv",
        oesophagusdndsgenesneutral="results/oesophagus/dnds_genes_neutral.csv",
        singlepatientdnds="results/dataforfigures/singlepatient_bins.csv",
        singlepatientdndsgenes="results/dataforfigures/singlepatient_bins_genes.csv"
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
        --maxarea {params.maxarea}  #2>> {log.out} 1>> {log.err}
        """
