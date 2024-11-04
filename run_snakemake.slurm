#!/bin/bash -e
#SBATCH -t 96:00:00 -N 1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rabdill@uchicago.edu

module load python

if [ ! -d "./venv" ]; then
    python -m venv venv
fi

source venv/bin/activate

if [ ! -f "./venv/bin/snakemake" ]; then
    pip install --upgrade pip
    pip install -r requirements.txt
fi

snakemake --use-singularity

# to submit:
# sbatch --job-name=next -o nextone.log run_snakemake.slurm
