---
tags:
  - minion
  - metagenome
  - dorado
---
# CLC metagenome basecall 

Chris Burgess

2023-11-06


**Update:** 11/20/2023 move the project to septoria project space.

Used the same steps I did for the lab_comp notebook except I excuted it on the ceres cluster

Details below did not change from `notebook_burgess/lab_comp/00_raw_data.md`


The data was pulled from ceres /project/soil_micro_lab/clc_metagenome.

This is the metagenomic sequencing of 1 soil sample. The sequencing data came in the old `.fast5` format so in the script `fast5_to_pod5_convert.sh`. The script converts the `.fast5` files to `.pod5` then I split the `.pod5` files by their channel.


The `.fast5` data is in file path is `00_raw_data/fast5/`

output of the pod5 conversion was parsed in [[01_base_calls]]

### Coverting the fast5 to pod5

```bash
mkdir pod5

pod5 convert fast5 fast5/*.fast5 --output pod5/ --one-to-one ./fast5

rm -rf fast5 ## fast5 are backed up so cleaning this up

```


### Spitting pod5 by channels

```bash
mkdir pod5_split_by_channel

pod5 view ./pod5 --include "read_id, channel" --output pod5_summary.tsv

pod5 subset ./pod5 --summary pod5_summary.tsv --columns channel --output ./pod5_split_by_channel

rm -rf pod5 ## cleaning up the unsplit data

```

### Created a list of files to run so I can reference it in my array job

`echo ../00_raw_data/pod5_split_by_channel/*.pod5 > pod5_files_to_run.txt`


