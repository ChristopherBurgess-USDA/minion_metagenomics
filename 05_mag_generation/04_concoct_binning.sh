#!/bin/bash
#SBATCH --job-name="clc_concoct"
#SBATCH -p mem
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 30 # number of logical cores/threads
#SBATCH --mem-per-cpu 30GB # Makign sure it has plenty of memory
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=BEGIN#SBATCH --mail-type=END

module purge
module load concoct
module load samtools
mkdir -p concoct_out concoct_out/fasta_bins

contig_path="../03_assembly/combined_reads.fasta"

cut_up_fasta.py $contig_path \
    -c 10000 \
    -o 0 \
    --merge_last \
    -b concoct_out/contigs_10K.bed \
    > concoct_out/contigs_10K.fa

samtools index coverage/minimap_mapping_sorted.bam \
    -@ $SLURM_CPUS_ON_NODE \
    > coverage/minimap_mapping_sort_index.bam

concoct_coverage_table.py concoct_out/contigs_10K.bed \
    coverage/minimap_mapping_sort_index.bam \
    > concoct_out/coverage_table.tsv


concoct --composition_file concoct_out/contigs_10K.fa \
    --coverage_file concoct_out/coverage_table.tsv \
    -b concoct_out/ \
    -l 2000 \
    --total_percentage_pca 90 \
    -t $SLURM_CPUS_ON_NODE

merge_cutup_clustering.py concoct_out/clustering_gt1000.csv \
    > concoct_out/clustering_merged.csv

extract_fasta_bins.py $contig_path \
    concoct_out/clustering_merged.csv \
    --output_path concoct_out/fasta_bins