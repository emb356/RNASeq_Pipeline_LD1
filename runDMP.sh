#!/bin/bash -l
#$ -cwd
# filename: runDMP.sh
# wrapper for second_demultiplex.rb
# executed in plateseq_post_demux.sh

workDir=$1 # project dir in /home/local/ARCS/ngs/Pipeline_ld1/Projects/PLATESEQ
sampleSheet=$2 # barcode csv

cd $workDir
mkdir $workDir/demultiplex

ruby /home/local/ARCS/ngs/Pipeline_ld1/scripts/second_demultiplex.rb $workDir $workDir/demultiplex/ unknown $sampleSheet 4

echo "done"
