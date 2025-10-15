---
tags:
  - metagenome
  - dorado
  - minion
  - slurm
  - ceres
  - gpu
---
## CLC metagenome basecall 

Chris Burgess
2023-11-06

**Updated 2023-11-07:** CHanged where the standard output files were getting stored modified the number of jobs being submitted at a time.


A lot of the steps here were already outlined in the lab_comp notebook; however, in order to run this on the ceres cluster(SLURM), I created an array submission job. This new job creates a new `.bam` file for each channel file. 
Here I use `dorado` to call bases on the channel split `.pod5` data. Since the samples were duplexed sequenced I can use the duplex option to basecall. The basecalls are stored as a `.bam` file. After duplex basecalls, I run the summary options by `dorado`.

dorado model used: `dna_r10.4.1_e8.2_400bps_sup@v4.2.0`

Sequence checking and QC was done in [[02_quality_filter]]


```bash
#!/bin/bash
#SBATCH --job-name="clc_metagenome"
#SBATCH -p gpu-low   # can switch to medium for all steps except dorado duplex
#SBATCH --time=01:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 2
#SBATCH --array 1-651%2

module load samtools
module load miniconda
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



```
