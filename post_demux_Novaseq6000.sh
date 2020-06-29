#!/bin/bash -l
#$ -S /bin/sh
#$ -cwd
# file name: post_demux_NovaSeq6000.sh
# starts data processing after demultiplexing
# executed in bclToFastqV2_Nova.sh or bclToFastqV2_Nextseq.sh

source $HOME/.bashrc

indexDir=$1
runid=$2
setting=$3

if [[ $indexDir == "" || $runid == ""  || $setting == "" ]]
    then
    echo "Error: Missing Argument :: $0 $1 $2 $3"
    echo "USAGE: $0 path/indexDir runID path/pipeline_global_settings_file "
    exit 1
fi

. $setting

cd $indexDir
fqRunDir=$( pwd -P )
for projectid in *_NOVASEQ; do
    projDir=${fqRunDir}/${projectid}
    cd $projDir
    echo "Starting post-demux analyses for $projectid in $projDir"

# if there are samples were sequenced in multiple lanes, merge them (while keeping the unmerged versions available, if necessary)
    echo "checking to merge fasts across lanes"
    numDupSamples=$(ls *_R1_*.gz | cut -f 1 -d "_" | sort | uniq -d | wc -l )
    if [ $numDupSamples -gt 0 ]; then
	mkdir unmergedFqs
	mv $( ls *fastq.gz | grep -iv undetermined ) unmergedFqs/
	cd unmergedFqs/
	ls *gz | cut -f1 -d "_" | uniq > $projDir/sampleID

	cd $projDir

        echo "merging..."
        cat $projDir/sampleID | parallel "sh /home/local/ARCS/ngs/Pipeline_ld1/scripts/mergeFastq.sh {} $PWD/unmergedFqs"

# move merged fastqs out of output dir
	mv $projDir/merged/*fastq.gz $projDir
	rmdir $projDir/merged
    fi
    echo "done merging"
    
    echo "starting fastQC"
    parallel "sh /home/local/ARCS/ngs/Pipeline_ld1/scripts/fastqQC.sh {} `pwd`" ::: *fastq.gz

    cd $fqRunDir    

    organism=`echo $projectid | cut -f5 -d '_' | tr '[:upper:]' '[:lower:]' `
    app=`echo $projectid | cut -f6 -d '_'`

    APP="RNA-seq/"
    echo "organism: $organism"
    echo "app: $app"
    
    if [[ $app =~ "RNA" || $app =~ "CDNA" || $app =~ "mirna-seq" ]]; then
        APP="RNA-seq/"
    elif [[ $app == "LIBRARY" ]]; then
	plrun=`echo $projectid | cut -f7 -d '_'`
	if [[ $plrun == "PLATESEQ" ]]; then
	    APP=$plrun
	fi
    else
        continue;
    fi

    # if this is a PLATEseq run, we need to run it's own pipeline script
    if [[ $APP == "PLATESEQ" ]]; then
	cmd="sh $LD1BASE/plateseq_post_demux.sh $projDir"
	echo $cmd
	$cmd
    fi


    # do the following, if this isn't a PLATESeq job
    if [[ $APP != "PLATESEQ" ]]; then

        echo "creating $Projects/$APP/$projectid/$runid/fastq"
        if [ ! -d $Projects/$APP/$projectid/$runid/fastq ]; then
    	mkdir -p $Projects/$APP/$projectid/$runid/fastq;
        fi
    
        ln -s $indexDir/$projectid/*fastq.gz $Projects/$APP/$projectid/$runid/fastq/
        echo "links done"
    
    
        cd $Projects/$APP/$projectid/$runid/
        
        for fastq in $indexDir/$projectid/*_R1_*gz; do
    	echo $fastq
            base_fq1=`basename $fastq`
            samplename=`echo $base_fq1| cut -f1 -d '_'`
            sampleID=`echo $base_fq1| cut -f2 -d '_'`
            lane=`echo $base_fq1| cut -f3 -d '_'`
            ending=`echo $base_fq1| cut -f5 -d '_'`
    
            base_fq2=$samplename"_"$sampleID"_"$lane"_R2_"$ending
    	    
    	if [ -e $Projects/$APP/$projectid/$runid/fastq/$base_fq2 ]; then
                ln_fq_1=$Projects/$APP/$projectid/$runid/fastq/$base_fq1
                ln_fq_2=$Projects/$APP/$projectid/$runid/fastq/$base_fq2
    	else
    	    ln_fq_1=$Projects/$APP/$projectid/$runid/fastq/$base_fq1
                ln_fq_2=""
    	fi
    	
    	echo $ln_fq_1
    	
    	if [[ $app =~ "RNA" || $app =~ "CDNA" || $app =~ "mirna-seq" ]]; then
    	    echo "sh $LD1BASE/RNA_pipeline_kallisto.sh $organism $ln_fq_1 $ln_fq_2" >> $Projects/$APP/$projectid/$runid/kall_command.list
    	else
    	    echo "ERROR APPLICATION: $app Unknown"
    	fi
        done 
        # run kallisto in parallel
        parallel --progress --jobs 25% --joblog $Projects/$APP/$projectid/$runid/kallisto_joblog.txt < $Projects/$APP/$projectid/$runid/kall_command.list

       # if all samples have finished processing with kallisto, create release files
       logCount=`ls -l logs/*txt | wc -l`
       sampCount=`ls -l $organism"_kallisto" | grep -c ^d`
       if [ "$logCount" -eq "$sampCount" ]; then
           cd $organism"_kallisto"
           cmd="sh $LD1BASE/kallistoRelease.sh"
           echo $cmd
           $cmd
       fi

    fi
done
