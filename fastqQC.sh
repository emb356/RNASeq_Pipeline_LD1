#!/bin/bash -l
#$ -cwd
# file name:fastqQC.sh
# runs fastQC for LD1 pipeline
# called in post_demux_Novaseq6000.sh

f=$1 #fastq.gz file
outputDir=$2

~/tools/FastQC_0.11.8/fastqc $f -o $outputDir



