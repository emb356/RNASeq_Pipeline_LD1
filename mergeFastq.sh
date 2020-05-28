#!/bin/bash
#$ -cwd
# file name: mergeFastq.sh
# merges fastqs across lanes 
# executed in post_demux_Novaseq6000.sh
 
sampleID=$1 # provide sample name
dir=$2 # dir of fastqs

if [ ! -d  merged ]; then mkdir merged; fi

for f in $dir/$sampleID"_"*_R1_*.fastq.gz; do
    f3=`echo $f | sed 's/_R1_/_R2_/'` # check PE-seq or SE-seq
    
    if [ -e $f3 ];then 
	 f4=`basename $f`
	 f5=`basename $f3`
	 echo $f4
	 echo $f5
	
    fi # PE seq needs two fastq files

    if [ ! -e $f3 ]; then
	f4=`basename $f`
	echo $f4
    fi

done

cat $dir/$sampleID"_"*_R1_*.fastq.gz > merged/$f4
cat $dir/$sampleID"_"*_R2_*.fastq.gz > merged/$f5

