---
tags:
  - minion
  - assembly
---
## Assemblying metagenomes

Chris Burgess

2023-11-27

### Overview

The next step in the metagenomic pipeline is to assembly the reads into longer contigs. I orginally, was going to compare results from a few different assembly pipelines; however, it seems like `flye` assembler is the only one who can use assemble metagenomes with just nanopore samples.

### Flye assembler

The Flye assembler was pretty straightforward to use. I did have an out of memory (OOM) error on the first attempt so I used a normal mem node with a lot of nodes.

**The output is in `03_assembly/flye_out/assembly.fasta`**

Here is the script: `01_flye_assembly.sh`

```bash
#!/bin/bash
#SBATCH --job-name="clc_flye"
#SBATCH -p mem # can switch to medium for all steps except dorado duplex
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 60
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END


temp_file_storage="/90daydata/septoria/clc_metagenome"

mkdir -p ${temp_file_storage}/03_assembly flye_out

module load flye


flye --nano-hq ../02_quality_filter/dorado_bc_nr_trim.fastq.gz \
    -t $SLURM_CPUS_ON_NODE \
    -o flye_out \
    --meta
```

### Assembly Results

The assembly stats look pretty promising:

```log
[2023-11-27 23:54:10] INFO: Assembly statistics:

	Total length:	731863414
	Fragments:	32700
	Fragments N50:	36445
	Largest frg:	1555424
	Scaffolds:	0
	Mean coverage:	13
```
**the N50 is defined as the sequence length of the shortest contig at 50% of the total assembly length.**

### Update 2023-12-28 and 2024-01-08

I realized that I had a long of long good quality reads that were not assembled into contigs. So I added a new script `02_alignment_compare.sh` script to add the unassembled reads to the contig file (code below). This should give us better results by allowing for sequences to be classified.

Since I ran into problems with using minimaps with `CONCOCT` I want to test out [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to see if my downstream problems are fixed and if my alignment results are as good as `minimap2`.

**Update: There was an error and the `samtools` was exporting a `fastq` sequences and not `fasta`. I fixed it and reran it.**


```bash
#!/bin/bash
#SBATCH --job-name="clc_contig_coverage"
#SBATCH -p mem
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 20 # number of logical cores/threads
#SBATCH --mem-per-cpu 64GB # Makign sure it has plenty of memory
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=END

module load minimap2
module load samtools
module load bowtie2

mkdir -p alignment_compare

###############################################################################
#*****************************************************************************#
#                            minimap2 alignment                               #
#*****************************************************************************#
###############################################################################

# # aligning unassembled reads to assembly.
# minimap2 -ax asm5 flye_out/assembly.fasta \
#     ../02_quality_filter/dorado_bc_nr_trim.fastq.gz \
#     -t $SLURM_CPUS_ON_NODE \
#     > alignment_compare/minimap_mapping.sam

# # adding unmapped/unassembled reads to assembled contigs
# cp flye_out/assembly.fasta alignment_compare/combined_reads_minimap.fasta

# samtools fasta -n -f 4 alignment_compare/minimap_mapping.sam \
#     >> alignment_compare/combined_reads_minimap.fasta

# # Checking to see how combined assembled and unassembled reads cover raw reads 
# # also needed for coverage information.
# minimap2 -ax asm5 alignment_compare/combined_reads_minimap.fasta \
#     ../02_quality_filter/dorado_bc_nr_trim.fastq.gz \
#     -t $SLURM_CPUS_ON_NODE \
#     -I64g \
#     > alignment_compare/minimap_mapped_all_reads.sam

###############################################################################
#*****************************************************************************#
#                            bowtie2 alignment                                #
#*****************************************************************************#
###############################################################################

# aligning unassembled reads to assembly.
bowtie2-build -f flye_out/assembly.fasta \
    --threads $SLURM_CPUS_ON_NODE \
    alignment_compare/bowtie_index

bowtie2 -x alignment_compare/bowtie_index \
    -U ../02_quality_filter/dorado_bc_nr_trim.fastq.gz \
    -p $SLURM_CPUS_ON_NODE \
    --very-sensitive \
    -S alignment_compare/bowtie_mapping.sam

# adding unmapped/unassembled reads to assembled contigs
cp flye_out/assembly.fasta alignment_compare/combined_reads_bowtie.fasta

samtools fasta -n -f 4 alignment_compare/bowtie_mapping.sam \
    --threads $SLURM_CPUS_ON_NODE \
    >> alignment_compare/combined_reads_bowtie.fasta


# Checking to see how combined assembled and unassembled reads cover raw reads 
# also needed for coverage information.
bowtie2-build -f alignment_compare/combined_reads_bowtie.fasta \
    --threads $SLURM_CPUS_ON_NODE \
    alignment_compare/bowtie_index_all_reads

bowtie2 -x alignment_compare/bowtie_index_all_reads\
    -U ../02_quality_filter/dorado_bc_nr_trim.fastq.gz \
    -p $SLURM_CPUS_ON_NODE \
    --very-sensitive \
    -S alignment_compare/bowtie_mapping_all_reads.sam

###############################################################################
#*****************************************************************************#
#                            alignment results                                #
#*****************************************************************************#
###############################################################################
printf "%0.s#" {1..100}; echo -e "\n"
printf "%0.s*" {1..100}; echo -e "\n"
printf "%0.s#" {1..100}; echo -e "\n"
echo "flagstat for minimap_before adding"

samtools flagstat alignment_compare/minimap_mapping.sam \
    --threads $SLURM_CPUS_ON_NODE

printf "%0.s#" {1..100}; echo -e "\n"
printf "%0.s*" {1..100}; echo -e "\n"
printf "%0.s#" {1..100}; echo -e "\n"
echo "flagstat for minimap all reads"

samtools flagstat alignment_compare/minimap_mapped_all_reads.sam \
    --threads $SLURM_CPUS_ON_NODE

printf "%0.s#" {1..100}; echo -e "\n"
printf "%0.s*" {1..100}; echo -e "\n"
printf "%0.s#" {1..100}; echo -e "\n"
echo "flagstat for bowtie2 before all reads"

samtools flagstat alignment_compare/bowtie_mapping.sam \
    --threads $SLURM_CPUS_ON_NODE

printf "%0.s#" {1..100}; echo -e "\n"
printf "%0.s*" {1..100}; echo -e "\n"
printf "%0.s#" {1..100}; echo -e "\n"
echo "flagstat for bowtie2 all reads"

samtools flagstat alignment_compare/bowtie_mapping_all_reads.sam \
    --threads $SLURM_CPUS_ON_NODE
```


### Update 2024-01-09

There was a hiccup when using `bowtie2` but I got a new verison of it install. I tried rerunning it but I was getting the same warnings of reads being too short. However, the data is quite large and all downstream analysis was taking a long time. I had an idea to try assembling all the reads that weren't able to assemble in the first round of `flye` assembly. I've updated the `01_flye_assembly.sh` to do two rounds of assembly and see how many reads were able to map to their combined assembly. 

```bash
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
```

#### Results for second assembly.

After running a second assembly, we only really saw a 3% increase in the number of mapped reads so I'm not entire sure if it was worth it at all.

```bash
####################################################################################################

****************************************************************************************************

####################################################################################################

flagstat for minimap_before adding
5456335 + 0 in total (QC-passed reads + QC-failed reads)
4783184 + 0 primary
50363 + 0 secondary
622788 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1586045 + 0 mapped (29.07% : N/A)
912894 + 0 primary mapped (19.09% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
####################################################################################################

****************************************************************************************************

####################################################################################################

flagstat for minimap after 2 assemblies adding
5512769 + 0 in total (QC-passed reads + QC-failed reads)
4783184 + 0 primary
52454 + 0 secondary
677131 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1744090 + 0 mapped (31.64% : N/A)
1014505 + 0 primary mapped (21.21% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

### Update 2024-01-16

The `bowtie2` alignment failed I think i need to do a `--large-index` option when building the index but I'm not sure it is worth it when `minimap2` is so much faster.

