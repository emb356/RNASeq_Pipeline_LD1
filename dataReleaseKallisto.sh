#!/bin/bash
#$ -cwd
# file name: dataReleaseKallisto.sh
# sends project data to the web release directory and runs md5sum

dir=$1 # proj ID 
fastqDir=$2 # fastq directory
filesDir=$3 # kallisto directory 
webdir=/home/local/ARCS/ngs/quicklinks_isilon/release

cd $webdir/
# create a directory using ProjID name in Web/release if it doesnt already exist
if [[ ! -e $dir ]]; then
    mkdir $dir
else
    echo "EXITING: directory already exists in /home/local/ARCS/ngs/quicklinks_isilon/release"
    exit 0
fi

# copy all files in the to_release directory to Web/release
if [[ -d $filesDir ]]; then
    for f in $filesDir"/to_release/*"; do
	cp -pr $f "$webdir/"$dir
    done
fi

# move fastqs to Web/release if fastqDir exists
if [[ -d $fastqDir ]]; then
   for f in $fastqDir/*gz; do
      tf=$(readlink -f $f)
      mv $tf $webdir/${dir}

   done
   
   cd $webdir
   cp -p /home/local/ARCS/ngs/Pipeline_ld1/scripts/READMEKALLISTO.txt ./README.txt
   chmod -R 775 $dir

   # this block just calls the md5 code
   # script that creates html index is executed in checkMD5_fastq.sh
   cd $webdir/${dir}
   parallel "sh /home/local/ARCS/ngs/Pipeline_ld1/scripts/checkMD5_fastq.sh {} {}'fastq_md5_temp.txt' $dir" ::: *fastq.gz 
fi

