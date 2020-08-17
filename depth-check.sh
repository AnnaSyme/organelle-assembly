!#usr/bin/env bash

#find where short reads didn't map to an assembly

samtools depth -a shortmapped.bam > samout

#http://www.htslib.org/doc/samtools-depth.html
#out file is a text file with contig name, assembly position and number of reads that mapped
#each line is a position

awk '$3 == "0"' < samout > awkout

#if field 3 is equal to zero

wc -l awkout

#num lines is number of positions with no reads that mapped
#can see contiguous regions
#can view in bam file in a Galaxy-JBrowse file


