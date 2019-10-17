#!/bin/bash 

export PATH=$PATH:/usr/local/bin

if [[ $1 == "" ]]; then
    setting="/home/local/ARCS/ngs/Pipeline_ld1/scripts/global_setting.sh"
else
    setting=$1
fi

. $setting

## number of threads used for bcl2fastq demultiplex
nt=16

runs=`ls $NovaSeqRuns`
for f in $runs; do
    
    if [ ! -d $NovaSeqRuns/$f ]; then
	## skip if $f isn't a directory
	continue
    fi

    if [[ $f == "MyRun" ]]; then
	continue
    fi

    g=`grep -w $f $StatusDir/basecall.complete`
#if run has not been processed, we will not find it in basecall.complete file  
    if [[ $g == "" ]]; then
	# check if the run is complete
	sampleSheetScript="$PIPEBASE/NovoSeqSampleSheet.rb"
	seqCompleteFname="SequenceComplete.txt"
#add run to basecall.complete if it is done, and run bcl2fastq
	echo $NovaSeqRuns/$f/$seqCompleteFname
	if [[ -e $NovaSeqRuns/$f/$seqCompleteFname ]]; then
	    cmd1="ruby $sampleSheetScript $SampleSheets/$f.csv $NovaSeqRuns/$f/SampleSheet.csv"
	    echo "$cmd1"
	    $cmd1
        fi
        if [[ -s $NovaSeqRuns/$f/SampleSheet.csv && -e $NovaSeqRuns/$f/$seqCompleteFname && -e $NovaSeqRuns/$f/CopyComplete.txt ]]; then
		    
	    echo "process $f"
	    echo -e "$f\t$NovaSeqRuns/$f/\t$FastqDir/LD1_$f" >> $StatusDir/basecall.complete




            cmd1="sh /home/local/ARCS/ngs/Pipeline_ld1/scripts/bclToFastqV2_Nova.sh -i $NovaSeqRuns/$f -o $FastqDir/"LD1_"$f -b ${bcl2fastqBin} -s /home/local/ARCS/ngs/Pipeline_ld1/scripts/global_setting.sh -n $nt -m 1"
		
	    echo "$cmd1"
            $cmd1
	fi
    fi
done


runs=`ls $NextSeqRuns`
for f in $runs; do
    if [ ! -d $NextSeqRuns/$f ]; then
	## skip if $f isn't a directory
	continue
    fi

    if [[ $f == "to_be_deleted" ]] || [[ $f == "single_cell_10x" ]] || [[ $f == "single_cell_plate" ]]; then
	continue
    fi
    
    g=`grep -w $f $StatusDir/basecall.complete`
    if [[ $g == "" ]]; then
	if [[ ! -s $SampleSheets/$f.csv ]]; then
	    echo "SampleSheet ($SampleSheets/$f.csv) is missing. Failed to start demultiplexing." ;
	    continue
	    ##exit ;
	fi
        cmd1="ruby $PIPEBASE/Hiseq2NextSeqSampleSheet.rb $SampleSheets/$f.csv $NextSeqRuns/$f/SampleSheet.csv"
        echo $cmd1
        $cmd1

	if [[ -s $NextSeqRuns/$f/SampleSheet.csv && -s  $NextSeqRuns/$f/RTAComplete.txt ]]; then
	    echo "process $f"
	    echo -e "$f\t$NextSeqRuns/$f/\t$FastqDir/LD1_$f" >> $StatusDir/basecall.complete
            cmd1="sh /home/local/ARCS/ngs/Pipeline_ld1/scripts/bclToFastqV2_Nova.sh -i $NextSeqRuns/$f -o $FastqDir/"LD1_"$f -b ${bcl2fastqBin} -s /home/local/ARCS/ngs/Pipeline_ld1/scripts/global_setting.sh -n $nt -m 1"
	    $cmd1
  	    echo $cmd1
        fi
    fi
done
exit 0
