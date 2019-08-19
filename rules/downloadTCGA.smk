rule downloadSTCGA:
    input:
        ascatfile="data/ascat_acf_ploidy.tsv"
    output:
        tcgadata="data/TCGA-combined-hg38.csv",
        maffile="data/TCGA-MAF.csv",
        cnvfile="data/TCGA-CNV.csv"
    log:
        out = "logs/downloadTCGA.out",
        err = "logs/downloadTCGA.err"
    shell:
        """
        Rscript R/downloadTCGA.R \
            --MAFfile {output.maffile} \
            --CNVfile {output.cnvfile} \
            --ascatcellularity {input.ascatfile} \
            --outputfile {output.tcgadata}  #2>> {log.out} 1>> {log.err}
        """
