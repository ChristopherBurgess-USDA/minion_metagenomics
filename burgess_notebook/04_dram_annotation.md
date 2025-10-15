---
tags:
  - dram
  - metagenome
---

## Annotating Metagenome

Chris Burgess

2023-12-06

### Overview

Here I annotate the metagenomes using the program [DRAM](https://github.com/WrightonLabCSU/DRAM) which is a progam which automates the annotation and distilation of the results. `DRAM` is designed to be used on MAGs; however, here I'm treating all sequences in the sample as 1 giant MAG.

After annotating with `DRAM` I still need to get read counts by getting a coverage number for each of my scaffods/contigs.

### DRAM setup

Before running DRAM I need to set it up, but this just needs to be done once. It took me a few attempts since I didn't anticipate the memory consumption when building the databases, aswell downloading the Uniref database timed out a few times and the pipeline crash.

It'd be really nice if this function had a resume function instead of making you redownload all the databases. After digging into the [DRAM github readme](https://github.com/WrightonLabCSU/DRAM?tab=readme-ov-file#getting-started-part-2-setup-databases) I found out it takes ALOT of ram so I made sure my memory node on SCINet had plenty of it.

```bash
#!/bin/bash
#SBATCH --job-name="clc_dram_setup"
#SBATCH -p mem
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 20 # number of logical cores/threads
#SBATCH --mem-per-cpu 60GB # Makign sure it has plenty of memory
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

dram_dir="/home/christopher.burgess/septoria/bin/dram"

module purge
module load miniconda

source activate DRAM

DRAM-setup.py prepare_databases \
    --output_dir $dram_dir \
    --threads $SLURM_CPUS_ON_NODE
```

### Annotating

Annotating the assembled contigs from `03_assembly/flye_out/assembly.fasta` was pretty straightforward straight from [their github readme](https://github.com/WrightonLabCSU/DRAM?tab=readme-ov-file#getting-started-part-3-usage).

```bash
#!/bin/bash
#SBATCH --job-name="clc_dram"
#SBATCH -p mem
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 20 # number of logical cores/threads
#SBATCH --mem-per-cpu 10GB # Makign sure it has plenty of memory
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module purge
module load miniconda

source activate DRAM

DRAM.py annotate -i ../03_assembly/flye_out/assembly.fasta \
    -o annotation

DRAM.py distill -i annotation/annotations.tsv \
    -o genome_summaries \
    --trna_path annotation/trnas.tsv \
    --rrna_path annotation/rrnas.tsv
```

It outputs a pretty neat `04_dram_annotation/genome_summaries/product.html` and a `04_dram_annotation/genome_summaries/metabolism_summary.xlsx`. They have a reall [nice readme on how to interprete the DRAM results](https://github.com/WrightonLabCSU/DRAM/wiki/4a.-Interpreting-the-Results-of-DRAM)

### Adding gene abundance

I still need to get a coverage of each of the scaffolds so I can do something like differential abundance testing on the genes. To do this I need to use `04_dram_annotation/annotation/scaffolds.fna` or the `04_dram_annotation/annotation/genes.fna` as a reference and map the reads back to it.

### Update 2023-12-28

I realized that I had a lot of long good quality reads that were not assembled into contigs. I'm rerunning `DRAM` with the updated fasta file to see if I get better results (code below).

```bash
## Same code just changed to output directory to have _combined
DRAM.py annotate -i ../03_assembly/combined_reads.fasta -o annotation_combined

DRAM.py distill -i annotation_combined/annotations.tsv \
    -o genome_summaries_combined \
    --trna_path annotation_combined/trnas.tsv \
    --rrna_path annotation_combined/rrnas.tsv
```

### Update: 2024-1-4

The previous DRAM run gave me a OOM error so I reran the job with updated resource request (below)

```bash
#!/bin/bash
#SBATCH --job-name="clc_dram"
#SBATCH -p mem
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 20 # number of logical cores/threads
#SBATCH --mem-per-cpu 62GB # Makign sure it has plenty of memory
#SBATCH --mail-user=christopher.burgess@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
```

### Update 2024-01-10

There was a new version of [DRAM](https://github.com/WrightonLabCSU/DRAM/releases/tag/v1.5.0) which has [CAMPER](https://github.com/WrightonLabCSU/CAMPER) (Curated Annotations for Microbial (Poly)phenol Enzymes and Reactions). However, when I updated teh DRAM package it broke with this error below.

```python
❯ DRAM-setup.py -h
  File "/home/christopher.burgess/.conda/envs/DRAM/bin/DRAM-setup.py", line 124
    set_db_locs_parser.add_argument('--camper_tar_gz_loc', default=None,
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
SyntaxError: invalid syntax. Perhaps you forgot a comma?
```

There was an error in the release `DRAM` which I was able to fix in the `DRAM-setup.py` script: [github issue here](https://github.com/WrightonLabCSU/DRAM/issues/326)

