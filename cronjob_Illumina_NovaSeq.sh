#!/bin/bash 
# file name: cronjob_Illumina_NovaSeq.sh
# checks NovaSeq and NextSeq run dirs for completed and unprocessed runs
# triggers pipeline
# executed with cron

export PATH=$PATH:/usr/local/bin

if [[ $1 == "" ]]; then
    setting="/home/local/ARCS/ngs/Pipeline_ld1/scripts/global_setting.sh"
else
    setting=$1
fi

. $setting

# number of threads used for bcl2fastq demultiplex
nt=16

runs=`ls $NovaSeqRuns`
for f in $runs; do
    
    if [ ! -d $NovaSeqRuns/$f ]; then
	# skip if $f isn't a directory
	continue
    fi

    if [[ $f == "MyRun" ]]; then
	continue
    fi

    g=`grep -w $f $StatusDir/basecall.complete`
# if run has not been processed, we will not find it in basecall.complete file  
    if [[ $g == "" ]]; then
	# check if the run is complete
	sampleSheetScript="$LD1BASE/NovaSeqSampleSheet.rb"
	seqCompleteFname="SequenceComplete.txt"
        # check if sample sheet exists
        if [[ -s $SampleSheets/$f.csv && ! -e $NovaSeqRuns/$f/$seqCompleteFname ]]; then
            echo "SampleSheet ($SampleSheets/$f.csv) found. Waiting for run to finish." ;
            continue
        fi
# add run to basecall.complete if it is done, and run bcl2fastq
	if [[ -e $NovaSeqRuns/$f/$seqCompleteFname && -s $SampleSheets/$f.csv ]]; then
	    echo $NovaSeqRuns/$f/$seqCompleteFname
            cmd1="ruby $sampleSheetScript $SampleSheets/$f.csv $NovaSeqRuns/$f/SampleSheet.csv"
	    echo "$cmd1"
	    $cmd1
        fi
        if [[ -s $NovaSeqRuns/$f/SampleSheet.csv && -e $NovaSeqRuns/$f/$seqCompleteFname && -e $NovaSeqRuns/$f/CopyComplete.txt ]]; then
		    
	    echo "process $f"
	    echo -e "$f\t$NovaSeqRuns/$f/\t$FastqDir/$f" >> $StatusDir/basecall.complete
            cmd1="sh $LD1BASE/bclToFastqV2_Nova.sh -i $NovaSeqRuns/$f -o $FastqDir/$f -b ${bcl2fastqBin} -s $LD1BASE/global_setting.sh -n $nt -m 1"

            echo "$cmd1"
            $cmd1
	fi
    fi
done


runs=`ls $NextSeqRuns`
for f in $runs; do
    if [ ! -d $NextSeqRuns/$f ]; then
	# skip if $f isn't a directory
	continue
    fi

    if [[ $f == "to_be_deleted" ]] || [[ $f == "single_cell_10x" ]] || [[ $f == "single_cell_plate" ]]; then
	continue
    fi
    
    g=`grep -w $f $StatusDir/basecall.complete`
    if [[ $g == "" ]]; then
        # check if sample sheet exists
	if [[ -s $SampleSheets/$f.csv && ! -e $NextSeqRuns/$f/RTAComplete.txt ]]; then
            echo "SampleSheet ($SampleSheets/$f.csv) found. Waiting for run to finish." ;
	    continue
	fi

# add run to basecall.complete if it is done, and run bcl2fastq
	if [[ -e $NextSeqRuns/$f/RTAComplete.txt && -s $SampleSheets/$f.csv ]]; then
	    echo $NextSeqRuns/$f/RTAComplete.txt
            cmd1="ruby $LD1BASE/NextSeqSampleSheet.rb $SampleSheets/$f.csv $NextSeqRuns/$f/SampleSheet.csv"
            echo $cmd1
            $cmd1
        fi
	if [[ -s $NextSeqRuns/$f/SampleSheet.csv && -s  $NextSeqRuns/$f/RTAComplete.txt ]]; then
	    echo "process $f"
	    echo -e "$f\t$NextSeqRuns/$f/\t$FastqDir/$f" >> $StatusDir/basecall.complete
            cmd1="sh $LD1BASE/bclToFastqV2_Nextseq.sh -i $NextSeqRuns/$f -o $FastqDir/$f -b ${bcl2fastqBin} -s $LD1BASE/global_setting.sh -n $nt -m 1"
	    $cmd1
  	    echo $cmd1
        fi
    fi
done
exit 0
