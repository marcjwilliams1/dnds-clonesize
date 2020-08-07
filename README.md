# Snakemake workflow: dnds-clonesize

This is a snakemake workflow to reproduce figures from [Measuring the distribution of fitness effects in somatic evolution by combining clonal dynamics with dN/dS ratios](https://elifesciences.org/articles/48714).

## Authors

* Marc J Williams (@marcjwilliams1)

## Usage

1. Clone this repo ```git clone https://github.com/marcjwilliams1/dnds-clonesize```
2. Install [snakemake](https://snakemake.readthedocs.io/en/stable/) if you don't have it installed
3. Make a dry-run to ensure that the pipeline will compile ```snakemake -n```
4. Run the pipeline with ```snakemake``` adding additional options if you want to run on a cluster (recommended).

Two scripts are provided that will run the pipeline on either a SGE or LSF cluster. If you are using a different job submission system you can use these as a starting point.

## Containers
2 containerized images are provided with all the software dependencies. The first includes all relevent R packages and is available [here](https://singularity-hub.org/collections/3462). The second includes the julia dependencies and is available [here](https://hub.docker.com/r/marcjwilliams1/julia_r_dnds). The snakemake pipeline will automatically pull these container or you can download them as follows:
```
singularity pull shub://marcjwilliams1/dnds-clonesize-R-container
singularity pull docker://marcjwilliams1/julia_r_dnds:v1.1
```

To run the pipeline you will need singularity installed on your cluster.

## Jupyter notebooks
Also included are jupyter notebooks in the `notebooks/` folder for each of the main steps. These are perhaps easier to parse than the pipeline if you're interested in looking at how a particular analysis was done although are by now a bit out of date.
