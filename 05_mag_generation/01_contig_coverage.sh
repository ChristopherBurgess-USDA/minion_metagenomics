#!/bin/bash
#SBATCH --job-name="clc_contig_coverage"
#SBATCH -p mem-low
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 20 # number of logical cores/threads
#SBATCH --mem-per-cpu 64GB # Makign sure it has plenty of memory
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module purge
module load samtools metabat minimap2


temp_file_storage="/90daydata/septoria/clc_metagenome"
mkdir -p coverage ${temp_file_storage}/05_mag_generation

contig_path="../03_assembly/combined_reads.fasta"
reads_path="../02_quality_filter/dorado_bc_nr_trim.fastq.gz"

minimap2 -ax asm5 ${contig_path} \
    ${reads_path} \
    -t $SLURM_CPUS_ON_NODE \
    -I64g \
    > coverage/minimap_mapping.sam

samtools view -bS coverage/minimap_mapping.sam > coverage/minimap_mapping.bam

samtools sort coverage/minimap_mapping.bam \
    --threads $SLURM_CPUS_ON_NODE \
    > coverage/minimap_mapping_sorted.bam

samtools index coverage/minimap_mapping_sorted.bam \
    -@ $SLURM_CPUS_ON_NODE > \
    coverage/minimap_mapping_sorted.bam.bai

jgi_summarize_bam_contig_depths \
    --outputDepth coverage/minimap_coverage.tsv \
    coverage/minimap_mapping_sorted.bam

mv coverage/minimap_mapping.sam ${temp_file_storage}/05_mag_generation

printf "%0.s#" {1..100}; echo -e "\n"
printf "%0.s*" {1..100}; echo -e "\n"
printf "%0.s#" {1..100}; echo -e "\n"
echo "flagstat for minimap"