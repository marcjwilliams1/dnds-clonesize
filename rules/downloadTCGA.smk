rule downloadTCGA:
    input:
        ascatfile="data/ascat_acf_ploidy.tsv"
    output:
        tcgadata="data/TCGA-combined-hg38.csv",
        maffile="data/TCGA-MAF.csv",
        cnvfile="data/TCGA-CNV.csv"
    shell:
        """
        module load R
        module load julia
        Rscript R/downloadTCGA.R \
            --MAFfile {output.maffile} \
            --CNVfile {output.cnvfile} \
            --ascatcellularity {input.ascatfile} \
            --outputfile {output.tcgadata}
        """
