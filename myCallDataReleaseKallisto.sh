#!/bin/bash
# file name: myCallDataReleaseKallisto.sh
# execute in ProjDir/RunDir to release project to HPC

# gets kallisto dir
# NOTE: only grabs one directory if there are multiple
kallistoDir=$(ls | grep "_kallisto" | head -n 1 )
runDir=$(pwd -P )

# run the kallistoRelease code if it hasn't already been run
if [ ! -z $kallistoDir ]; then
    if [ ! -f $kallistoDir/summary.csv ]; then
        echo "running kallisto release"
	cd $kallistoDir
	sh /home/local/ARCS/ngs/Pipeline_ld1/scripts/kallistoRelease.sh
	cd $runDir
    fi
    kallistoFolder=$( ls -d $PWD/*_kallisto | head -n 1 )
else
    kallistoFolder=""
fi

projId=$( basename $( dirname $PWD ) )

fastqDir=$runDir"/fastq"
echo "sh /home/local/ARCS/ngs/Pipeline_ld1/scripts/dataReleaseKallisto.sh $projId $fastqDir $kallistoFolder"
sh /home/local/ARCS/ngs/Pipeline_ld1/scripts/dataReleaseKallisto.sh $projId $fastqDir $kallistoFolder 

