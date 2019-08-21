# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


configfile: "config.yaml"
report: "report/workflow.rst"


# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

shell.executable("/bin/bash")
shell.prefix("source ~/.bash_profile; ")

figs=[1,2,3,4,5]
suppfigs=[1,2,3,4,5,6,7,8,9,10,11]

rule all:
    input:
        expand("Figures/Figure{FIG}.pdf", FIG = figs),
        expand("Figures/FigureS{FIG}.pdf", FIG = suppfigs)

include: "rules/downloadTCGA.smk"
include: "rules/ModellingNormalTissue.smk"
include: "rules/CalculatedNdS-normal.smk"
include: "rules/fitdNdSnormal.smk"
include: "rules/CalculatedNdS-tcga.smk"
include: "rules/ModellingCancer.smk"
include: "rules/GenerateFigures.smk"
