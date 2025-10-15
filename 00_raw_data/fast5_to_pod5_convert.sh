
### Coverting the fast5 to pod5

mkdir -p pod5

pod5 convert fast5 fast5/*.fast5 --output pod5/ --one-to-one ./fast5

rm -rf fast5 ## fast5 are backed up so cleaning this up

### Spitting pod5 by channels

mkdir -p pod5_split_by_channel

pod5 view ./pod5 --include "read_id, channel" --output pod5_summary.tsv

pod5 subset ./pod5 --summary pod5_summary.tsv --columns channel --output ./pod5_split_by_channel

rm -rf pod5 ## cleaning up the unsplit data

### Created a list of files to run so I can reference it in my array job

`echo ../00_raw_data/pod5_split_by_channel/*.pod5 > pod5_files_to_run.txt`

