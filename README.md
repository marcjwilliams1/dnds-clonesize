# Snakemake workflow: dnds-clonesize

This is a snakemake workflow to reproduce figures from "Measuring the distribution of fitness effects in somatic evolution by combining clonal dynamics with dN/dS ratios"

## Authors

* Marc J Williams (@marcjwilliams1)

## Usage

1. Clone this repo ```git clone https://github.com/marcjwilliams1/dnds-clonesize```
2. Install [snakemake](https://snakemake.readthedocs.io/en/stable/) if you don't have it installed
3. Make a dry-run to ensure that the pipeline will compile ```snakemake -n```
4. Run the pipeline with ```snakemake``` adding additional options if you want to run on a cluster (recommended). See the snakemake docs for this.

## Singularity
Singularity image is available [here](https://singularity-hub.org/collections/3462). The snakemake pipeline will automatically pull the container or you can download it as follows:
```
singularity pull shub://marcjwilliams1/dnds-clonesize-R-container
```

## Jupyter notebooks
Also included are jupyter notebooks in the `notebooks/` folder for each of the main steps. These are perhaps easier to parse than the pipeline if you're interested in looking at how a particular analysis was done.
