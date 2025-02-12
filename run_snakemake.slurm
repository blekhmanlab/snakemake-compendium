#!/bin/bash -e
#SBATCH -t 96:00:00 -N 1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rabdill@uchicago.edu

module load python/3.11.9
module load apptainer

if [ ! -d "./venv" ]; then
    python -m venv venv
fi

source venv/bin/activate

if [ ! -f "./venv/bin/snakemake" ]; then
    pip install --upgrade pip
    pip install -r requirements.txt
fi

snakemake --default-resources slurm_account=pi-blekhman slurm_partition=blekhman --set-resources sra_prefetch:slurm_partition=build

# to submit:
# sbatch --job-name=next -o nextone.log run_snakemake.slurm
