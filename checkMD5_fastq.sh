#!/bin/bash -l
#$ -cwd
# file name: checkMD5_fastq.sh
# runs md5sum on fastqs and puts data live on web
# executed in dataReleaseKallisto.sh
 
infile=$1
outfile=$2
dir=$3
webdir="/home/local/ARCS/ngs/quicklinks_isilon/release"

md5sum $infile >> $outfile

echo "complete" > $infile"_md5Done.txt"

samples=`ls *"_md5Done.txt" | wc -l`

fastqs=`ls *fastq.gz | wc -l`


# check all md5 files are complete
# update release link when all files have been created
if [ "$samples" -eq "$fastqs" ]; then 
    cat *temp.txt > fastq_md5.txt
    rm *temp.txt
    rm *_md5Done.txt
    cd $webdir/
    # contains python 3 env
    # NOTE: still must change permissions on HPC (chmod -R 775 projectDir)
    source activate HPC_Pipeline
    /home/local/ARCS/ngs/Pipeline_ld1/scripts/index_gen.py -d $(pwd -P)/$dir
fi
