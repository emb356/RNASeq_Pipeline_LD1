#!/bin/bash -l
#$ -cwd
# file name: generateSummary.sh
# creates alignment percentages table
# executed in kallistoRelease.sh

outDir=$1 # [organism]_kallisto dir
cd $outDir 
echo "SampleID,Number_Reads,Pseudoaligned_Reads,Pseudoaligned_Ratio" > summary.csv
chmod 770 summary.csv

sample_dirs=`ls --hide="summary.csv"`

for i in $sample_dirs
do
    sampleID=`echo $i | cut -f1,3 -d '_'`
    numberReads=`grep "n_processed" $i/run_info.json | cut -f2 -d ":" | cut -c 2- | rev | cut -c 2- | rev`
    pseudoalignedReads=`grep "n_pseudoaligned" $i/run_info.json | cut -f2 -d ":" | cut -c 2- | rev | cut -c 2- | rev`
    pseudoalignedRatio=`grep "p_pseudoaligned" $i/run_info.json | cut -f2 -d ":" | cut -c 2- | rev | cut -c 2- | rev`
#    uniqueReads=`grep "n_unique" $i/run_info.json | cut -f2 -d ":" | cut -c 2- | rev | cut -c 2- | rev`
    echo $sampleID", "$numberReads", "$pseudoalignedReads", "$pseudoalignedRatio>>summary.csv

done

cd ..
