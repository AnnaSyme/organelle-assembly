#!/usr/bin/env bash

#A script to filter and trim illumina reads using fastp

#to filter illumina reads for quality and trim adapters:
#--n_base_limit 3 - discard read/pair if > 3Ns
#--length_required 130
#--average_qual 35 - want very high qual for polishing

fastp --in1 R1.fq.gz --out1 R1_fastp.fq.gz --in2 R2.fq.gz --out2 R2_fastp.fq.gz --verbose \
--adapter_fasta adapters.fasta --n_base_limit 3 --length_required 130 \
--average_qual 35 --threads 16
