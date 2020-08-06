#!/bin/sh
#$ -cwd            # Set the working directory for the job to the current directory
#$ -pe smp 1       # Request 1 core
#$ -l h_rt=120:0:0 # Request 24 hour runtime
#$ -l h_vmem=1G    # Request 1GB RAM
#$ -j y
#$ -o logs

module load singularity

CLUSTER_CMD="qsub -cwd -l h_rt={cluster.time} -l h_vmem={cluster.mem} -o {cluster.output} -j y -N {cluster.name} -pe smp {threads}"

snakemake --jobs 80 \
  --cluster-config cluster_sge.yaml \
  --cluster "${CLUSTER_CMD}" \
  --use-singularity \
  --singularity-args "--bind /home/william1 --bind /work" --keep-going
