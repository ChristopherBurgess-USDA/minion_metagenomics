#!/bin/bash
#SBATCH --job-name="clc_maxbin"
#SBATCH -p mem
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 30 # number of logical cores/threads
#SBATCH --mem-per-cpu 30GB # Makign sure it has plenty of memory
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=END

module purge
module load maxbin
module load seqkit
contig_path="../03_assembly/combined_reads.fasta"
mkdir -p maxbin_out

seqkit fq2fa ../02_quality_filter/dorado_bc_nr_trim.fastq.gz \
    -j $SLURM_CPUS_ON_NODE \
    -o ../02_quality_filter/dorado_bc_nr_trim.fasta

run_MaxBin.pl -contig $contig_path \
    -out maxbin_out/maxbin_ \
    -reads ../02_quality_filter/dorado_bc_nr_trim.fasta \
    -thread $SLURM_CPUS_ON_NODE
