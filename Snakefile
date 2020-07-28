# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


configfile: "config.yaml"
report: "report/workflow.rst"


# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

shell.executable("/bin/bash")
shell.prefix("source ~/.bash_profile; ")

figs=[1,2,3,4,5]

#list of oesophagus samples
OES_SAMPLES=["PD36806","PD36712","PD30272","PD30986","PD30987","PD30274", "PD30988", "PD30273","PD31182"]
SSB_genes=["global", "notch1", "tp53"]

rule all:
    input:
        expand("Figures/Figure{FIG}.pdf", FIG = figs),
        expand("Figures/Figure1-S{FIG}.pdf", FIG = [1]),
        expand("Figures/Figure2-S{FIG}.pdf", FIG = [1,2,3,4]),
        expand("Figures/Figure3-S{FIG}.pdf", FIG = [1]),
        expand("Figures/Figure5-S{FIG}.pdf", FIG = [1,2,3,4,5]),
        "Figures/Extra/Figure-Extra-binsize.pdf",
        "Figures/Extra/Figure-Extra-SimsDistributionRegression.pdf",
        directory(expand("results/oesophagus/SSBfiles/{oes_sample}/", oes_sample=OES_SAMPLES)),
        "results/oesophagus/SSBresults/SSBdnds_results.csv",
        "results/dataforfigures/oesophagusfit-SSB.csv",
        "results/dataforfigures/brmsfit.Rdata",
        "results/oesophagus/sitednds_genes_hotspots.csv",
        "results/dataforfigures/brmsfit-sites.Rdata",
        "results/dataforfigures/simulation-clonesizefit.Rdata",
        "results/dataforfigures/data-clonesizefit.Rdata",
        "results/dataforfigures/data-clonesizefit-models.Rdata",
        "results/simulations/clonesize_hitchikers.csv",
        "results/dataforfigures/stemcell_simulation_examplefits_distribution-exp.csv",
        "results/simulations/clonesize_overtime-dist.csv"

include: "rules/ModellingNormalTissue.smk"
include: "rules/CalculatedNdS-normal.smk"
include: "rules/fitdNdSnormal.smk"
include: "rules/GenerateFigures.smk"
include: "rules/clonesizefit.smk"
