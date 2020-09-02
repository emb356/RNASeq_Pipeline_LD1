#!/bin/bash -lv
# file name: plateseq_post_demux.sh
# processes plateseq for second demultiplexing

projDir=$1
projDir=$( readlink -f $projDir )
pdir=$(basename $projDir)
# where the Project data will go
secondStep="/home/local/ARCS/ngs/Pipeline_ld1/Projects/PLATESEQ/${pdir}"

# make project directory in PLATESEQ
# make fastq folder based on Indexing#
indexdir=$(echo $projDir | cut -d/ -f9)
mkdir -p ${secondStep}/${indexdir}
webdir=/home/local/ARCS/ngs/quicklinks_isilon/release

cd $secondStep
cd $indexdir

for f in $projDir/*.gz; do
ln -s $f .
done

for samp in $(ls *.gz | cut -f 1 -d "_" | uniq); do 
    # make a dir for each plate name
    mkdir ${samp}
    cd ${samp}
    # symlink each plate fastq to its dir
    ln -s ../${samp}*gz .
    # copy barcode sample sheet
    rpi=$( echo $samp | cut -c 1-4 )
    cp /home/local/ARCS/ngs/Pipeline_ld1/References/${rpi}.csv .
 
    # set the lane correctly
    r1=$( ls PS*_R1_* )
    lane=$( echo $r1 | cut -f 3 -d "_" | sed "s:L00::" )
    sed -i "s:8,S:${lane},S:" ${rpi}.csv

    SampleSheet=${secondStep}/${indexdir}/${samp}/${rpi}.csv

    echo "sh /home/local/ARCS/ngs/Pipeline_ld1/scripts/runDMP.sh `pwd` $SampleSheet" >> ../plate_command.list
    cd ${secondStep}/${indexdir}
done

parallel --progress --joblog plate_joblog.txt < plate_command.list

# create a directory using ProjID name in Web/release if it doesnt already exist
if [[ ! -e $webdir/${pdir} ]]; then
    mkdir -p $webdir/${pdir}
fi

# move demultiplexed fastqs to release
mv PS*/PS*_{A,G,T,C}*gz $webdir/${pdir}

cd $webdir/${pdir}
fastqs=`ls *fastq.gz | wc -l`

# only intitiate md5sum creation if all 768 fastq.gzs are in release folder
if [ "$fastqs" -eq 768 ]; then
    echo "running md5sum for $webdir/${pdir}"
    parallel "sh /home/local/ARCS/ngs/Pipeline_ld1/scripts/checkMD5_fastq.sh {} {}'fastq_md5_temp.txt' $pdir" ::: *fastq.gz   
fi



