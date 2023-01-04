#!/bin/bash -e
#SBATCH -t 45:00:00 -N 1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=email@EMAIL_GOES_HERE.com

source venv/bin/activate

snakemake --use-singularity --jobs 45 --slurm --default-resources slurm_account=ACCOUNT_HERE slurm_partition=PARTITION_HERE
