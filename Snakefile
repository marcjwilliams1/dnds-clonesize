# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


configfile: "config.yaml"
report: "report/workflow.rst"


# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

shell.executable("/bin/bash")
shell.prefix("source ~/.bash_profile; ")

figs=[1,2,3,4]

rule all:
    input:
        expand("Figures/Figure{FIG}.pdf", FIG = figs)


include: "rules/Rule1-ModellingNormalTissue.smk"
include: "rules/Rule2-CalculatedNdS-normal.smk"
include: "rules/fitdNdSnormal.smk"
include: "rules/GenerateFigures.smk"
include: "rules/CalculatedNdS-tcga.smk"
