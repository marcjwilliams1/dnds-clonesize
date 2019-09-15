rule brms:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
    output:
        "results/dataforfigures/brmsfit.Rdata"
    threads: 4
    shell:
        """
        module load R/3.5.3
        module load gcc
        Rscript R/fitbrms.R \
            --oesophagusdata {input.oesophagusdata} \
            --oesophagusmetadata {input.oesophaguspatientinfo} \
            --output {output}
        """
