#!/bin/bash -e
#SBATCH -t 96:00:00 -N 1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rabdill@uchicago.edu

source venv/bin/activate

snakemake --use-singularity --jobs 45 --slurm --default-resources slurm_account=blekhman slurm_partition=blekhman
