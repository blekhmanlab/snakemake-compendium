#!/bin/bash -e
#SBATCH -t 1:00:00 -N 1
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH -p build
#SBATCH -A pi-blekhman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rabdill@uchicago.edu

if [ ! -d "./venv" ]; then
    python -m venv venv
fi

source venv/bin/activate

if [ ! -f "./venv/bin/snakemake" ]; then
    pip install --upgrade pip
    pip install -r requirements.txt
fi

sbatch --job-name=${SLURM_JOB_NAME}run -o ${SLURM_JOB_NAME}run.$(date +"%s").log ./run_snakemake.slurm
