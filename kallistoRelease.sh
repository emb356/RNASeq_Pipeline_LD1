#!/bin/bash
#$ -cwd
# filename: kallistoRelease.sh
# creates all files for release
# executed in RNA_pipeline_kallisto.sh
# can execute manually in [organism]_kallisto directory

source $HOME/.bashrc

# set variables that link to the mounted HPC scripts
setting="/home/local/ARCS/ngs/Pipeline_ld1/scripts/global_setting.sh"
. $setting

echo "copy files"

LD1BASE="/home/local/ARCS/ngs/Pipeline_ld1/scripts"

# generate summary.csv
sh $LD1BASE/generateSummaryTable.sh ./

sample_dirs=`ls --hide="summary.csv"`

dirName="to_release"

if [ ! -d $dirName ]; then mkdir -p $dirName; fi

kallisto_out="to_release/kallisto_raw_output"
if [ ! -d $kallisto_out ]; then mkdir -p $kallisto_out; fi
chmod 775 $kallisto_out

# copy the tsvs to use in Rscripts 
for i in $sample_dirs; do
    cp $i/abundance.tsv $kallisto_out/$i.abundance.tsv 
    if [ -e $i/pseudoalignments.bam ]; then
        mv $i/pseudoalignments.bam $dirName/$i.pseudoalignments.bam ; mv $i/pseudoalignments.bam.bai $dirName/$i.pseudoalignments.bam.bai
    fi
    mv $i $kallisto_out
    chmod 775 $kallisto_out/$i
done 

# create concatenated estimated counts tables for transcripts and genes for all samples
cd $dirName
# set the .Rprofile for the Rscripts
source activate HPC_Pipeline
cp -p $LD1BASE/.Rprofile .
Rscript $LD1BASE/kallistoCountsTable.R
Rscript $LD1BASE/tximport_gene_counts.R

mkdir transcript_level
chmod 775 transcript_level 
mv countstable.txt transcript_level/est_counts_transcipts_kallisto.txt
mv tpmtable.txt transcript_level/tpm_values_transcripts_kallisto.txt
cp $LD1BASE/READMEKALLISTO.txt ./README.txt

cd ..

rm $kallisto_out/*.tsv

cp summary.csv $dirName
