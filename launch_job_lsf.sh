#!/bin/sh
#BSUB -J snakemake
#BSUB -n 1
#BSUB -R rusage[mem=4]
#BSUB -W 200:00
#BSUB -eo logs/%J.stderr

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pseudobulkQC

module load singularity

CLUSTER_CMD=("bsub -n {threads} -M {cluster.mem} -o {cluster.output} -J {cluster.name} -W {cluster.time}")

snakemake --jobs 80 \
  --cluster-config cluster_lsf.yaml \
  --cluster "${CLUSTER_CMD}" \
  --use-singularity \
  --singularity-args "--bind /home/william1 --bind /work" --keep-going
