---
tags:
  - metagenome
  - mags
---
## Annotating MAGS

Chris Burgess

2023-12-06

### Overview

Here instead of annotating the metagenome as a whole, we want to see if we can generate metagenome assembly genomes (MAGs). To do that I'm adapting a [kbase nature protocol for SCINet](https://www.nature.com/articles/s41596-022-00747-x)

Since we don't have access to the article, the authors shared a linked to it [here](https://www.nature.com/articles/s41596-022-00747-x.epdf?sharing_token=aYwPliWVzKmhGa9POJvputRgN0jAjWel9jnR3ZoTv0PTMzshsGuS6o-gH6M2KRreQjVgBUxXmKyxa7cIrbqb9yneiDuho6oMD1wKWKtyH2NMXjjojD9Q_05hBPjK8RlE3-fjD2VjCEjEVdttYZs_x3lckEs0oM8YKWmjGc5Cokw%3D)

The protocol breaks down into a few steps:

1) generate coverage of contigs (aka abundance)
2) Bin contigs with [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/), [Maxbin2](https://sourceforge.net/p/maxbin/code/ci/master/tree/), and [CONCOCT](https://github.com/BinPro/CONCOCT)
3) Clean up binned contigs with [DAS_Tool](https://github.com/cmks/DAS_Tool)
4) Use [CheckM](https://github.com/Ecogenomics/CheckM) to select good bins
5) Annotate with [DRAM](https://github.com/WrightonLabCSU/DRAM/wiki)

### 1 Generate coverage of contigs

Both `metabat2` and `maxbin2` can use a coverage file to improve binning. We also need the coverage information to give us the abundance of annotated genes aswell. To do this we use `minimap2` to map our reads to our assembled contigs. Then use a `metabat2` tool to fine the mean coverage of all contigs.

`flye` does output some coverage information in their `assembly_info.txt`; however, there are two problems with this coverage estimate. First, it is just the mean; however, more importantly, since `flye` was unable to assemble all the reads, there were many reads that were >= 2000bp that were not included see the update to `03_assembly.md` notes.

Also I tried to use `bbmap` to map reads to assembled contigs; however, the reads themselves are too long and I couldn't figure out how to get it to work right.

**File:** `01_contig_coverage.sh`

**Update 2024-01-05** I updated the script to covert the `.sam` output to `.bam` before sorting and indexing. Hopefully, this fixed the issues with `CONCOCT`. 

```bash
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
module load samtools
module load metabat
module load minimap2


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
```

**Note:** I ran into an error with `samtools`: `[E::sam_parse1] no SQ lines present in the header`. This happened because the reference (assembled contigs) was too large so I increased `minimap2` resources to fix the issue: [More details in minimap2 FAQ](https://github.com/lh3/minimap2/blob/master/FAQ.md#3-the-output-sam-doesnt-have-a-header)


### 2 Binning Contigs

Binning contigs using `metabat2`, `maxbin2`, and `concoct` is done in 3 different scripts

#### 02_metabat_binning.sh

**Update: took longer than 2 hours to run needed a non low queue**

```bash
#!/bin/bash
#SBATCH --job-name="clc_metabat"
#SBATCH -p mem
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 30 # number of logical cores/threads
#SBATCH --mem-per-cpu 30GB # Makign sure it has plenty of memory
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module purge
module load metabat
contig_path="../03_assembly/combined_reads.fasta"

mkdir -p metabat_out

metabat2 -i ${contig_path} \
    -a coverage/minimap_coverage.tsv \
    -o metabat_out/metabat_bin \
    -m 2000 \
    -v
```

#### 03_maxbin_binning.sh

**Update:** I got an error running this one (below), need to trouble shoot this.

```bash
Failed to get Abundance information for contig [a65bb967-c6e3-4525-9b14-c03d63ce70be;67da4d9d-7335-4133-a188-3df122f83bb5] in file [maxbin_out/maxbin_.contig.tmp.abund1]
```

Going to try it with feeding it the reads (in fasta format) and having it generate coverage using `Bowtie2`

```bash
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
```

#### 04_concoct_binning.sh

**Updated the script, needed to index the sorted `.bam` file.**

**Update 2: This fix didn't work, I think I [found the solution](https://github.com/BinPro/CONCOCT/issues/307#issuecomment-1579464556)** Basically, the problem is that you need to convert the orginal `.sam` into a `.bam` format and before sorting and indexing. Also you need the orignal presorted `.bam` in the same folder (maybe?). I added these steps to `01_contig_coverage.sh` 

Still didn't work, I think I need to try `bowtie2` mapping instead of `minimap2`

```bash
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
```

