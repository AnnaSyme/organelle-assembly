#!/usr/bin/env/ bash


#convert stats from seqkit stats into combined tables 

```
cut -f1,4,5,6,8 illumina_raw_fastp.stats > all_illumina_stats.tsv
cut -f1,4,5,6,8 illumina_read_stats.tsv >> all_illumina_stats.tsv

cut -f1,4,5,6,8 nanopore_raw.stats > all_nanopore_stats.tsv
cut -f1,4,5,6,8 nanopore_read_stats.tsv >> all_nanopore_stats.tsv

cut -f1,4,5,6,8  assembly_stats.tsv > all_assembly_stats.tsv
```




