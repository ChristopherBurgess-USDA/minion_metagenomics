#!/bin/bash
#SBATCH --job-name="clc_dram"
#SBATCH -p mem
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 20 # number of logical cores/threads
#SBATCH --mem-per-cpu 62GB # Makign sure it has plenty of memory
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module purge
module load miniconda

source activate DRAM

DRAM.py annotate -i ../03_assembly/combined_reads.fasta -o annotation_combined

DRAM.py distill -i annotation_combined/annotations.tsv \
    -o genome_summaries_combined \
    --trna_path annotation_combined/trnas.tsv \
    --rrna_path annotation_combined/rrnas.tsv