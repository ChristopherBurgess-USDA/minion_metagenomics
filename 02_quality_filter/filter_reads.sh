#!/bin/bash
#SBATCH --job-name="clc_metagenome"
#SBATCH -p mem768  # can switch to medium for all steps except dorado duplex
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 60
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

temp_file_storage="/90daydata/septoria/clc_metagenome"

module load samtools fastqc bbtools multiqc

mkdir -p temp ${temp_file_storage}/01_base_calls ${temp_file_storage}/02_quality_filter

## see what my average duplex rate across all channels: ~16% 
awk '/[[:digit:]]+\.[[:digit:]]+%/ {sub("%", ""); sum+=$NF count++} END { printf ("average percentage of duplex reads is %s", sum/count)}' ..//01_base_calls/jobs/slurm-*.out

## Merging channels by file

samtools merge --threads 60 -o dorado_bc.bam ../01_base_calls/temp/*.bam

samtools view -d dx:0 -d dx:1 -Ofastq --threads 60 -h -e 'length(seq)<=100000' dorado_bc.bam | gzip -9 > dorado_bc_nr.fastq.gz

## Quality Filtering

## Trimming
bbduk.sh in=dorado_bc_nr.fastq.gz \
    out=dorado_bc_nr_trim.fastq.gz \
    qtrim=rl \
    trimq=10 \
    minlen=200

## Removing without trimming
bbduk.sh in=dorado_bc_nr.fastq.gz \
    out=dorado_bc_nr_filter.fastq.gz \
    maq=10 \
    minlen=200