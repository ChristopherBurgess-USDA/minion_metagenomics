---
tags:
  - minion
  - fastqc
---
## Checking sequence quality and filtering

Chris Burgess

2023-11-20

### Goals

After doing the duplex runs I wanted to check the sequence quality out and so different filtering methods to see which would gain or loss information.

Job paramters are as listed:

```bash
#!/bin/bash
#SBATCH --job-name="clc_metagenome"
#SBATCH -p mem768  # can switch to medium for all steps except dorado duplex
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 60
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

temp_file_storage="/90daydata/septoria/clc_metagenome"

module load samtools
module load fastqc
module load bbtools
module load multiqc

mkdir -p temp ${temp_file_storage}/01_base_calls ${temp_file_storage}/02_quality_filter
```

I used `awk` to see what my average duplex rate across all channels: ~16% 

`awk '/[[:digit:]]+\.[[:digit:]]+%/ {sub("%", ""); sum+=$NF count++} END { printf ("average percentage of duplex reads is %s", sum/count)}' ..//01_base_calls/jobs/slurm-*.out`

### Merging channel files

Here I merge and filter out only the duplex and simplex with no duplex reads

the duplex and simplex reads are IDed by the tags `dx:1` and `dx:0`

```bash
samtools merge --threads 60 -o dorado_bc.bam ../01_base_calls/temp/*.bam

samtools view -d dx:0 -d dx:1 -Ofastq --threads 60 -h -e 'length(seq)<=100000' dorado_bc.bam | gzip -9 > dorado_bc_nr.fastq.gz
```

### Quality filtering

I used `bbduk` which is part of `bbtools` module to do two types of quality filtering:

1) Triming reads from either side till the read has a min Qscore of 10
2) Removing reads that have a Qscore less than 10
3) remove reads that are short (<200bp)

```bash
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
```

### Sequence Quality

Sequence quality was assessed using `fastq` and `multiqc` and resulting file is saved as `multiqc_report.html`

Looking over the results I think trimming reads is the way to go, median read length doesn't change and we get to keep slighly more reads trimming.

### Moving Forward

Next step is the do assembly of trimmed quality reads `dorado_bc_nr_trim.fastq.gz`

### Future work

It might be useful to explore the alternative filtering tools, why, and how to QC ONT sequencing data. 