#!/bin/bash
#SBATCH --job-name="clc_dram_setup"
#SBATCH -p mem
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 20 # number of logical cores/threads
#SBATCH --mem-per-cpu 60GB # Makign sure it has plenty of memory
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

dram_dir="/home/christopher.burgess/septoria/bin/dram"

module purge
module load miniconda

source activate DRAM

DRAM-setup.py prepare_databases \
    --output_dir $dram_dir \
    --threads $SLURM_CPUS_ON_NODE