#!/bin/bash -l
#$ -cwd
# file name: trigRNA_pipeline_kallisto.sh 
# manually triggers kallisto pipeline for a project
# USAGE above fastq directory: sh /home/local/ARCS/ngs/Pipeline_ld1/scripts/trigRNA_pipeline_kallisto.sh [organism] fastq

RNACODE="/home/local/ARCS/ngs/Pipeline_ld1/scripts"

genome=$1 # provide genome name
dir=$2 # fastq directory

dir=$( readlink -f $dir )

for f in $dir/*_R1_*.fastq*; do
    f2=`echo $f | sed 's/_R1_/_R2_/'` # get second fastq
    echo $f
    if [ -e $f2 ];then 
	echo $f2
        echo "sh $RNACODE/RNA_pipeline_kallisto.sh $genome $f $f2" >> $dir/kall_command.list; 

    fi

# SE (deprecated)
    if [ ! -e $f2 ]; then
        echo "sh $RNACODE/RNA_pipeline_kallisto.sh $genome $f" >> $dir/kall_command.list; #SE seq

    fi

done

echo $genome

# run GNU parallel for the pipeline
parallel --progress --jobs 25% --joblog $dir/kallisto_joblog.txt < $dir/kall_command.list
