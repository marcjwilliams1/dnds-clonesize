rule downloadTCGASNV:
    output:
        maffile="data/TCGA-MAF.csv",
    shell:
        """
        module load R
        Rscript R/downloadTCGASNV.R \
            --MAFfile {output.maffile}
        """

rule downloadTCGACNV:
    output:
        cnvfile="data/TCGA-CNV.csv"
    shell:
        """
        module load R
        Rscript R/downloadTCGACNV.R \
            --CNVfile {output.cnvfile}
        """

rule combineTCGA:
    input:
        ascatfile="data/ascat_acf_ploidy.tsv",
        maffile="data/TCGA-MAF.csv",
        cnvfile="data/TCGA-CNV.csv"
    output:
        tcgadata="data/TCGA-combined-hg38.csv",
    shell:
        """
        module load R
        Rscript R/combineTCGA.R \
            --MAFfile {input.maffile} \
            --CNVfile {input.cnvfile} \
            --ascatcellularity {input.ascatfile} \
            --outputfile {output.tcgadata}
        """
