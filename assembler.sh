#!/usr/bin/env bash

{

#set up logfile
now=$(date +%Y_%m_%d_%H_%M)
# print info in logfile
echo "log file for $0"
# $0 is the script name
echo "year_month_day_hour_minute $now"


#A script to assemble an organelle genome

#............................................................................
# How to use

# activate your conda environment with the tools needed
# bash assembler.sh -b baits -g 160000 R1 R2 nano

#to filter illumina reads for quality and trim adapters:
#--n_base_limit 3 - discard read/pair if > 3Ns
#--length_required 130
#--average_qual 35 - want very high qual for polishing

#fastp --in1 R1.fq.gz --out1 R1_fastp.fq.gz --in2 R2.fq.gz --out2 R2_fastp.fq.gz --verbose \
#--adapter_fasta adapters.fasta --n_base_limit 3 --length_required 130 \
#--average_qual 35 --threads 16

#............................................................................
# Defaults

set -e #exit if a command exits with non zero status
script=$(basename $0) #script name less file path
baits=""
genome_size=160000  #set as input arg
threads=16
target_bases=40000000 #set as input arg

#............................................................................
# Functions

function msg {
  echo -e "$*"
}
# $* is args passed in, -e means interpret backslashes

function err {
  echo "Error: $*" 1>&2
}
# redirects stdout to stderr

function banner {
  printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

function msg_banner {
  printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
  msg "$*"
  printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

function usage {
  banner
  msg "$script\n   Assemble an organelle genome with long and short reads."
  msg "Author:\n  Anna Syme, anna.syme@gmail.com"
  msg "Usage:\n   $script [options] R1 R2 nano"
  msg "Parameters:"
  msg "   Illumina R1 reads, trimmed, filtered   R1.fq.gz"
  msg "   Illumina R2 reads, trimmed, filtered   R2.fq.gz"
  msg "   Nanopore reads, raw                    nano.fq.gz"
  msg "Options:"
  msg "   -h                     Show this help"
  msg "   -t NUM                 Number of threads (default=16)"
  msg "   -g NUM                 genome size in bp (default=160000)"
  msg "   -s NUM                 target bases in bp (default=40000000)"
  msg "   -b FILE                bait sequences file name"
  msg "Example:"
  msg "   $script -t 8 R1.fq.gz R2.fq.gz minion.fq.gz"
  banner
  exit 1
  #exits script
}
#...........................................................................
# Parse the command line options

#loops through, sets variable or runs function
#have to add in a colon after the flag, if it's looking for an arg
#e.g. t: means flag t is set, then look for the arg (eg 32) to set $threads to
# : at the start disables default error handling for this bit
#instead it will label an odd flag as a ?
#the extra : at the end means it will label missing args as :

while getopts ':ht:g:s:b::' opt ; do
  case $opt in
    h)
      usage
      ;;
    t)
      threads=$OPTARG
      ;;
    g)
      genome_size=$OPTARG
      ;;
    s)
      target_bases=$OPTARG
      ;;
    b)
      baits=$OPTARG
      ;;
    \?)
      echo "Invalid option '=$OPTARG'"
      exit 1
      ;;
    :)
      echo "Option '-$OPTARG' requires an argument"
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))
#remove all options that has been parsed by getopts
#so that $1 refers to next arg passed to script
#e.g. a positional arg

if [ $# -ne 3 ]; then
  msg "\n **Please provide three input parameters** \n"
  usage
fi
#if number of pos params is not 3, print usage msg

#...........................................................................
#start
msg "\n"
msg_banner "now running $script"

#...........................................................................
#check inputs

#give variable names to the inputs
R1=$1
R2=$2
nano_raw=$3

msg "This script will use:"
msg "   Illumina reads R1:   $R1"
msg "   Illumina reads R2:   $R2"
msg "   Nanpore reads:       $nano_raw"
msg "   Bait sequences:      $baits"
msg "   Genome size:         $genome_size"
msg "   Target bases:        $target_bases"
msg "   Threads:             $threads"

#...........................................................................

conda env export --name assemblerenv > assemblerenv.yml
#saves conda env with tools and versions

#...........................................................................
#these read files will be created during run
R1_extracted=R1_extracted.fq.gz
R2_extracted=R2_extracted.fq.gz
R1_extracted_subset=R1_extracted_subset.fq.gz
R2_extracted_subset=R2_extracted_subset.fq.gz
nano_extracted=nano_extracted.fq.gz
nano_extracted_long=nano_extracted_long.fq.gz
nano_extracted2=nano_extracted2.fq.gz
nano_extracted_long2=nano_extracted_long2.fq.gz

#...........................................................................
#these assemblies will be created during run

#round 1 assembly
assembly_flye1=assembly_flye1.fasta
assembly_flye1_racon1=assembly_flye1_racon1.fasta
assembly_flye1_racon2=assembly_flye1_racon2.fasta

#round 2 assembly
assembly_flye2=assembly_flye2.fasta
assembly_flye2_racon1=assembly_flye2_racon1.fasta
assembly_flye2_racon2=assembly_flye2_racon2.fasta
assembly_flye2_racon_pilon1=assembly_flye2_racon_pilon1.fasta
assembly_flye2_racon_pilon2=assembly_flye2_racon_pilon2.fasta

#raven 
assembly_raven=assembly_raven.fasta
assembly_raven_pilon=assembly_raven_pilon.fasta

#unicycler
assembly_unicycler=assembly_unicycler.fasta

#miniasm
assembly_miniasm=assembly_miniasm.fasta
assembly_miniasm_minipolished=assembly_miniasm_minipolished.fasta
assembly_miniasm_minipolished_pilon1=assembly_miniasm_minipolished_pilon1.fasta

#...........................................................................
msg_banner "now extracting organelle nanopore reads from all reads - round 1"

minimap2 -a -x map-ont -t $threads $baits $nano_raw | 
samtools fastq -0 $nano_extracted -n -F 4 -
#map the reads to a reference set of organelle genes (baits)
#this pipes the output to samtools fastq,
#which extracts the fastq reads from the alignment
#the flag -F 4 means exclude unmapped reads (i.e., non organelle reads)

#...........................................................................
msg_banner "now keeping only the longest nanopore organelle reads - round 1"

filtlong --length_weight 1 --mean_q_weight 0 --window_q_weight 0 \
--target_bases $target_bases $nano_extracted | gzip > $nano_extracted_long

##keeps set of longest reads to required cov (eg X250) based on genome size (eg 160000)
##read score based on length only due to relative score weights
##alternative: filter by length and quality e.g.
##filtlong --target_bases 32000000 reads_nano_cp.fq.gz \
##| gzip > reads_nano_cp_filtered.fq.gz

#...........................................................................
msg_banner "now assembling nanopore organelle reads - round 1"

flye --nano-raw $nano_extracted_long --genome-size $genome_size \
--out-dir flye-out --threads $threads
#using uncorrected reads as read correction can lead to errors
#other flye options had little effect here - keep haplotpyes, meta, trestle

cp flye-out/assembly.fasta $assembly_flye1

#...........................................................................
msg_banner "now polishing flye assembly with long reads - round 1"

#round 1. make overlaps file, then run racon
minimap2 -x map-ont -t $threads $assembly_flye1 $nano_extracted_long \
| gzip > overlaps1.paf.gz

racon --threads $threads $nano_extracted_long overlaps1.paf.gz \
$assembly_flye1 > $assembly_flye1_racon1

#round 2
minimap2 -x map-ont -t $threads $assembly_flye1_racon1 $nano_extracted_long \
| gzip > overlaps2.paf.gz

racon --threads $threads $nano_extracted_long overlaps2.paf.gz \
$assembly_flye1_racon1 > $assembly_flye1_racon2

#here, further rounds of polishing made little difference
#option to add in medaka polishing here

#...........................................................................
msg_banner "now extracting organelle nanopore reads from all reads - round 2"

#baits file is now the first polished assembly
#need to increase minimum match value
#otherwise too many reads are extracted and they don't assemble

minimap2 -m 5000 -a -x map-ont -t $threads $assembly_flye1_racon2 $nano_raw | 
samtools fastq -0 $nano_extracted2 -n -F 4 -

#...........................................................................
msg_banner "now keeping only the longest nanopore organelle reads - round 2"

filtlong --length_weight 1 --mean_q_weight 0 --window_q_weight 0 \
--target_bases $target_bases $nano_extracted2 | gzip > $nano_extracted_long2

#...........................................................................
msg_banner "now assembling nanopore organelle reads - round 2"

flye --nano-raw $nano_extracted_long2 --genome-size $genome_size \
--out-dir flye-out2 --threads $threads

cp flye-out2/assembly.fasta $assembly_flye2

#...........................................................................
msg_banner "now polishing flye assembly with long reads - round 2"

#round 1. make overlaps file, then run racon
minimap2 -x map-ont -t $threads $assembly_flye2 $nano_extracted_long2 \
| gzip > overlaps1.paf.gz

racon --threads $threads $nano_extracted_long2 overlaps1.paf.gz \
$assembly_flye2 > $assembly_flye2_racon1

#round 2
minimap2 -x map-ont -t $threads $assembly_flye2_racon1 $nano_extracted_long2 \
| gzip > overlaps2.paf.gz

racon --threads $threads $nano_extracted_long2 overlaps2.paf.gz \
$assembly_flye2_racon1 > $assembly_flye2_racon2

#...........................................................................
msg_banner "now extracting illumina organelle reads from all reads"

minimap2 -a -x sr $assembly_flye2_racon2 $R1 $R2 \
| samtools fastq -1 $R1_extracted -2 $R2_extracted -F 0x4 -f 0x2 -

#extract organelle reads only by mapping to long-read assembly
#chose assembly made with longest nanopore reads
#needs end dash for stdin
#more info https://broadinstitute.github.io/picard/explain-flags.html
#samtools flags: -F4 exclude unmapped reads, -f2 include properly paired reads

#...........................................................................
msg_banner "now creating subset of illumina organelle reads"

rasusa -i $R1_extracted -i $R2_extracted --coverage 250 \
--genome-size $genome_size -o $R1_extracted_subset -o $R2_extracted_subset

#downsample to required coverage, x250

#...........................................................................
msg_banner "now polishing flye assembly with illumina"

#round 1pilon polish
bwa index $assembly_flye2_racon2

bwa mem -t $threads $assembly_flye2_racon2 $R1_extracted_subset $R2_extracted_subset \
| samtools sort > flye_aln1.bam

samtools index flye_aln1.bam

samtools faidx $assembly_flye2_racon2

pilon --genome $assembly_flye2_racon2 --frags flye_aln1.bam \
--output assembly_flye2_racon_pilon1 \
--fix bases --mindepth 0.5 --changes --threads $threads --verbose

#fix bases, not contig breaks in case that makes incorrect breaks

#round 2 pilon polish

bwa index $assembly_flye2_racon_pilon1

bwa mem -t $threads $assembly_flye2_racon_pilon1 $R1_extracted_subset $R2_extracted_subset \
| samtools sort > flye_aln1.bam

samtools index flye_aln1.bam

samtools faidx $assembly_flye2_racon_pilon1 

pilon --genome $assembly_flye2_racon_pilon1  --frags flye_aln1.bam \
--output assembly_flye2_racon_pilon2 \
--fix bases --mindepth 0.5 --changes --threads $threads --verbose

#...........................................................................
msg_banner "now running raven assembler"

#note: this uses the same long and short reads as flye+racon+pilon
#long and short reads are not re-extracted. 
#two rounds of racon polishing are included by default

raven --graphical-fragment-assembly raven.gfa \
-t $threads $nano_extracted_long2 > $assembly_raven

#...........................................................................
msg_banner "now polishing raven assembler with pilon"

#round 1pilon polish
bwa index $assembly_raven

bwa mem -t $threads $assembly_raven $R1_extracted_subset $R2_extracted_subset \
| samtools sort > raven_aln1.bam

samtools index raven_aln1.bam

samtools faidx $assembly_raven

pilon --genome $assembly_raven --frags raven_aln1.bam \
--output assembly_raven_pilon \
--fix bases --mindepth 0.5 --changes --threads $threads --verbose

#...........................................................................
msg_banner "now running unicycler assembler"

unicycler -1 $R1_extracted_subset -2 $R2_extracted_subset \
-l $nano_extracted_long2 -o unicycler \
--threads $threads --no_rotate --keep 2

#using fastp filtered illumina reads and longest nano organelle reads
#--keep 2 will keep final files but also SAM
#--no_rotate means don't rotate replicons to certain start pos

cp unicycler/assembly.fasta $assembly_unicycler

#note: this uses the same reads as the flye but not raven assembly

#...........................................................................
msg_banner "now running miniasm assembler"

#map reads to themselves - miniasm needs this as input
minimap2 -x ava-ont $nano_extracted_long2 $nano_extracted_long2 \
| gzip -1 > reads_overlaps.paf.gz

#assemble with miniasm - needs reads and overlaps (the paf file)
miniasm -f $nano_extracted_long2 reads_overlaps.paf.gz > miniasm.gfa

#convert gfa to fasta file of unitigs
awk '/^S/{print ">"$2"\n"$3}' miniasm.gfa > $assembly_miniasm

#...........................................................................
msg_banner "now polishing miniasm assembly with long reads"

#minipolish, uses racon and attempts to circularise

minipolish -t $threads $nano_extracted_long2 miniasm.gfa > minipolished.gfa

awk '/^S/{print ">"$2"\n"$3}' minipolished.gfa > $assembly_miniasm_minipolished

#alternative: use racon polish only
#racon polish round 1
#minimap2 -x map-ont $assembly_miniasm $nano_extracted_long | gzip > overlaps1.paf.gz
#racon --threads $threads $nano_extracted_long overlaps1.paf.gz \
#$assembly_miniasm > $assembly_miniasm_racon1

#racon polish round 2
#minimap2 -x map-ont $assembly_miniasm_racon1 $nano_extracted_long | gzip > overlaps2.paf.gz
#racon --threads $threads $nano_extracted_long overlaps2.paf.gz \
#$assembly_miniasm_racon1 > $assembly_miniasm_racon2

#...........................................................................
msg_banner "now polishing miniasm assembly with short reads"

#pilon polish round 1
bwa index $assembly_miniasm_minipolished

bwa mem -t $threads $assembly_miniasm_minipolished $R1_extracted_subset $R2_extracted_subset \
| samtools sort > mini_pilon_aln1.bam

samtools index mini_pilon_aln1.bam

samtools faidx $assembly_miniasm_minipolished

pilon --genome $assembly_miniasm_minipolished --frags mini_pilon_aln1.bam \
--output assembly_miniasm_minipolished_pilon1 \
--fix bases --mindepth 0.5 --changes --threads $threads --verbose

#fix bases, not contig breaks in case that makes incorrect breaks

#pilon polish round 2
#option. made no difference here
#bwa index $assembly_miniasm_racon_pilon1
#bwa mem -t $threads $assembly_miniasm_racon_pilon1 $R1_extracted_subset $R2_extracted_subset \
#| samtools sort > mini_pilon_aln2.bam
#samtools index mini_pilon_aln2.bam
#samtools faidx $assembly_miniasm_racon_pilon1
#pilon --genome $assembly_miniasm_racon_pilon1 --frags mini_pilon_aln2.bam \
#--output assembly_miniasm_racon_pilon2 \
#--fix bases --mindepth 0.5 --changes --threads $threads --verbose

#...........................................................................
msg_banner "now mapping long reads to final flye assembly"

minimap2 -ax map-ont $assembly_flye2_racon_pilon2 $nano_extracted_long2 \
| samtools sort -o longmapped.bam

samtools depth -a longmapped.bam > longdepths
awk '$3 == "0"' < longdepths > longawk
wc -l longawk > longmapped_positions_zero_reads.txt

#...........................................................................
msg_banner "now mapping short reads to final flye assembly"

bwa index $assembly_flye2_racon_pilon2

bwa mem -t $threads $assembly_flye2_racon_pilon2 $R1_extracted_subset $R2_extracted_subset \
| samtools sort -o shortmapped.bam

samtools depth -a shortmapped.bam > shortdepths
awk '$3 == "0"' < shortdepths > shortawk
wc -l shortawk > shortmapped_positions_zero_reads.txt

mkdir bams
mv longmapped.bam shortmapped.bam bams/
mv longmapped_positions_zero_reads.txt shortmapped_positions_zero_reads.txt bams/

#...........................................................................
msg_banner "organise assemblies and graphs"

#move all assemblies into one folder
mkdir assemblies
cp $assembly_flye1 $assembly_flye1_racon1 $assembly_flye1_racon2 \
$assembly_flye2 $assembly_flye2_racon1 $assembly_flye2_racon2 \
$assembly_flye2_racon_pilon1 $assembly_flye2_racon_pilon2 \
$assembly_raven $assembly_raven_pilon \
$assembly_unicycler \
$assembly_miniasm \
$assembly_miniasm_minipolished $assembly_miniasm_minipolished_pilon1 \
assemblies/

#move graphs into one folder
mkdir graphs
cp flye-out/assembly_graph.gfa graphs/flye1-assembly.gfa
cp flye-out2/assembly_graph.gfa graphs/flye2-assembly.gfa
cp unicycler/assembly.gfa graphs/unicycler.gfa
cp miniasm.gfa raven.gfa graphs/

#...........................................................................
msg_banner "organise extracted reads"


seqkit stats $nano_extracted $nano_extracted_long $nano_extracted2 $nano_extracted_long2 \
-Ta > nanopore_read_stats.tsv

mkdir organelle-reads-nanopore
mv $nano_extracted $nano_extracted_long \
$nano_extracted2 $nano_extracted_long2 organelle-reads-nanopore/

seqkit stats $R1_extracted $R2_extracted $R1_extracted_subset $R2_extracted_subset \
-Ta > illumina_read_stats.tsv

mkdir organelle-reads-illumina
mv $R1_extracted $R2_extracted \
$R1_extracted_subset $R2_extracted_subset organelle-reads-illumina/

#...........................................................................
msg_banner "now calculating stats"

seqkit stats assemblies/* -Ta > assembly_stats.tsv

mkdir stats
mv nanopore_read_stats.tsv illumina_read_stats.tsv assembly_stats.tsv stats/

#...........................................................................
msg_banner "get results"

mkdir results
mv organelle-reads-nanopore organelle-reads-illumina assemblies graphs bams stats results/

#...........................................................................
msg_banner "Script finished!"

} 2>&1 | tee logfile.txt 
