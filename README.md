# organelle-assembly

A script to assemble a plastid or mitochondrial genome from long and short reads. 

## To run:

```
conda activate bio
bash assembler.sh -b baits.fasta -g 160000 -s 40000000 R1.fasta R2.fasta nano.fq.gz
```

## Inputs:

* `-b` baits file, e.g. gene sequences of related species
* `-g` expected genome size 
* `-s` target bases (for Filtlong - e.g. coverage (250) x genome size)
* R1 and R2 illumina reads, already trimmed and filtered (see fastp.sh script)
* Nanopore reads, raw

## How it works (abridged):

* Uses a baits file to extract nanopore organelle reads (e.g. mitochondrial or chloroplast) from **all** the sequencing reads (e.g. nuclear, mitochondrial, chloroplast)
* Assembles these reads (Flye); polishes assembly (Racon)
* Uses this assembly as the new baits file to re-extract nanopore organelle reads
* Assembles (Flye); polishes assembly (Racon)
* Uses this assembly as baits to extract **illumnina** organelle reads
* Uses these reads to polish assembly (Pilon)
* Further assemblies for comparison: Raven, Miniasm, Unicycler
* Reports read and assembly stats

## Where are the results?

* The `results` folder has assemblies, assembly graphs, extracted reads, read-mapping bam files, and read/assembly stats.
* The script run and screen output is saved as `logfile.txt`

## Tools:

Installed with conda.

```
minimap2
samtools
filtlong
flye
racon
raven
fastp
rasusa
bwa
pilon
unicycler
miniasm
minipolish
mummer
seqkit
```

See packagae-list.txt for full details. 







