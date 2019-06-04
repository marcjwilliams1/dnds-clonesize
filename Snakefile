# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


configfile: "config.yaml"
report: "report/workflow.rst"


# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

shell.executable("/bin/bash")
shell.prefix("source ~/.bash_profile; ")

rule all:
    input:
        fig1 = "Figures/Figure1.pdf",
        fig2 = "Figures/Figure2.pdf"


include: "rules/Rule1-ModellingNormalTissue.smk"
include: "rules/Rule2-CalculatedNdS-normal.smk"
include: "rules/fitdNdSnormal.smk"
include: "rules/GenerateFigures.smk"
include: "rules/CalculatedNdS-tcga.smk"
