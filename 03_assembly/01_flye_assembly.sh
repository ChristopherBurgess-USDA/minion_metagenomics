#!/bin/bash
#SBATCH --job-name="clc_flye"
#SBATCH -p mem
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 56 # number of logical cores/threads
#SBATCH --mem-per-cpu 26GB # Makign sure it has plenty of memory
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=END


temp_file_storage="/90daydata/septoria/clc_metagenome"

mkdir -p ${temp_file_storage}/03_assembly flye_out flye_out_unmapped temp

module load flye samtools minimap2

###############################################################################
#*****************************************************************************#
#                            Assembly Round 1                                 #
#*****************************************************************************#
###############################################################################
flye --nano-hq ../02_quality_filter/dorado_bc_nr_trim.fastq.gz \
    -t $SLURM_CPUS_ON_NODE \
    -o flye_out_1 \
    --meta

# aligning unassembled reads to assembly.

minimap2 -ax asm5 flye_out/assembly.fasta \
    ../02_quality_filter/dorado_bc_nr_trim.fastq.gz \
    -t $SLURM_CPUS_ON_NODE \
    -I64g \
    > temp/minimap_mapping.sam

# Generating fastq of unmapped reads
samtools fastq -n -f 4 temp/minimap_mapping.sam \
    > temp/unassembled_reads.fastq

###############################################################################
#*****************************************************************************#
#                            Assembly Round 2                                 #
#*****************************************************************************#
###############################################################################

flye --nano-hq temp/unassembled_reads.fastq \
    -t $SLURM_CPUS_ON_NODE \
    -o flye_out_unmapped \
    --meta

cp flye_out/assembly.fasta flye_assembly.fasta

# Renaming contigs to avoid name conflicts

sed -i 's/contig_/contig_t1_/g' flye_assembly.fasta 

cp flye_out_unmapped/assembly.fasta temp.fasta

sed -i 's/contig_/contig_t2_/g' temp.fasta

cat temp.fasta >> flye_assembly.fasta

rm -rf temp.fasta

###############################################################################
#*****************************************************************************#
#                            Results from 2 assemblies                        #
#*****************************************************************************#
###############################################################################

minimap2 -ax asm5 flye_assembly.fasta \
    ../02_quality_filter/dorado_bc_nr_trim.fastq.gz \
    -t $SLURM_CPUS_ON_NODE \
    -I26g \
    > temp/minimap_mapping_take2.sam

printf "%0.s#" {1..100}; echo -e "\n"
printf "%0.s*" {1..100}; echo -e "\n"
printf "%0.s#" {1..100}; echo -e "\n"
echo "flagstat for minimap_before adding"

samtools flagstat temp/minimap_mapping.sam \
    --threads $SLURM_CPUS_ON_NODE

printf "%0.s#" {1..100}; echo -e "\n"
printf "%0.s*" {1..100}; echo -e "\n"
printf "%0.s#" {1..100}; echo -e "\n"
echo "flagstat for minimap after 2 assemblies adding"

samtools flagstat temp/minimap_mapping_take2.sam \
    --threads $SLURM_CPUS_ON_NODE