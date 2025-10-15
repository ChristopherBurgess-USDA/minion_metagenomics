#!/bin/bash
#SBATCH --job-name="clc_metagenome"
#SBATCH -p gpu-low   # can switch to medium for all steps except dorado duplex
#SBATCH --time=01:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 2
#SBATCH --array 1-651%2

## Due to GPU contraints on the cluster job was run as an array job

module load samtools miniconda
source activate /project/soil_micro_lab/conda/minion

mkdir -p temp

dorado_path="../../bin/dorado-0.4.3-linux-x64/bin"

file="$(awk -v awkvar=$SLURM_ARRAY_TASK_ID '{print $awkvar}' ../00_raw_data/pod5_files_to_run.txt)"

base="$(basename -s .pod5 $file)"
echo "My current job file: " $base


$dorado_path/dorado duplex \
    $dorado_path/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 \
    $file > temp/$base.bam
  
mv $file ../00_raw_data/pod5s_done/$base.pod5