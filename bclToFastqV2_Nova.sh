#!/bin/bash -l
#$ -cwd
uname -a

source $HOME/.bashrc
echo "Start $0:`date`"

setting="/home/local/ARCS/ngs/Pipeline_ld1/scripts/global_setting.sh"
USAGE="Usage: $0 -i RunDir  -o OutDir  -s setting [-n num_threads] [-m num of mismatches to be allowed in barcodes] [-c #reads to be placed in each fastq file] [-S flag to indicate demultiplexing of Read1 only to get reads stats] [-F flag to force repeat demultiplexing of run] "

nt=16 # default 16 threads
force="no"
SEforce="no"
barcodemismatch=1
chopFastqinMillions="250000000"

#doest work because lacks permissions
#bcl2fastqBin="/home/local/ARCS/ngs/mnt/data/bin/bcl2fastq/bcl2fastq2-v2.20.0.422/bin/bcl2fastq"
#bcl2fastqBin="/usr/local/bin/bcl2fastq"

while getopts i:o:b:s:n:m:c:hFS opt
  do
  case "$opt" in
      i) RunDir="$OPTARG";;
      o) OutDir="$OPTARG";;
      b) bcl2fastqBin="$OPTARG";;
      s) setting="$OPTARG";;
      n) nt="$OPTARG";;
      F) force="yes";;
      S) SEforce="yes";;
      m) barcodemismatch="$OPTARG";;
      c) chopFastqinMillions="$OPTARG";;
      h) echo $USAGE; exit 1
  esac
done

if [[ $RunDir == "" || $OutDir == ""  || $setting == "" ]]
    then
    echo $USAGE
    exit 1
fi

. $setting 
absIN=`readlink -f $RunDir`  ## not robust enough, need a better method
absOUT=`readlink -f $OutDir`
INTENSITY=$absIN/Data/Intensities/
BCALL=$absIN/Data/Intensities/BaseCalls/
runName=`basename $absIN`



if [[ $force == "yes" ]];then absOUT=$absOUT"FORCE" ; fi	##Incase bcl-fastq needs to be re-done for anyreason, use the -F flag to Force
if [[ $SEforce == "yes" ]];then absOUT=$absOUT"SE" ; fi		##Incase conversion needs to be on Read1, e.g. to get demultiplexing stats only.

echo "outdir: $OutDir     $absOUT "
if [[ ! -e $absOUT ]]; then mkdir $absOUT ; fi
chmod 770 $absOUT
pushd $absOUT

#### Demultiplex ##################################################################################

numlanes=`grep "<Lane>" $absIN/RunInfo.xml | wc -l`   
numReads=`grep "Read Num" $absIN/RunInfo.xml | wc -l`
lastread=`grep "Read Num" $absIN/RunInfo.xml | tail -1 `

sampleSheet=$SampleSheets/$runName.csv
echo "Using Sample Sheet : $sampleSheet"
cat $sampleSheet
if [[ ! -s $sampleSheet ]]; then  echo "SampleSheet is missing. Failed to start demultiplexing." ;  exit ; fi
demultiplexout=$absOUT/demultiplex
if [[ ! -e $demultiplexout ]]; then mkdir $demultiplexout ; fi
cp $sampleSheet $demultiplexout/$runName.csv
cd $demultiplexout
DIRNAME0="Indexing0"
DIRNAME6="Indexing6" # "DualIndex"	
DIRNAME7="Indexing7"
DIRNAME8="Indexing8"
DIRNAME8HALOPLEX="Indexing8HALOPLEX"
DIRNAME16="Indexing16"

lines_16_csv=0
lines_0_csv=0
lines_6_csv=0
lines_7_csv=0
lines_8_csv=0
lines_8_haloplex=0

function run_regular_casava(){
        local mask=$1
        local tsv=$2
        local indexing=$3
	local ends=$4 
	local csv=""
	local outdir=""

	if [[ $indexing == "NI" ]];then
		fastqcluster="40000000"
		head -1 $tsv > $tsv.0.csv
		grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ gsub(/-/,"",$5); print $0; } ' >> $tsv.0.csv
	else	
		head -1 $tsv > $tsv.0.csv
		head -1 $tsv > $tsv.6.csv
		head -1 $tsv > $tsv.7.csv
		head -1 $tsv > $tsv.8.csv
		head -1 $tsv > $tsv.8.haloplex.csv
		head -1 $tsv > $tsv.16.csv
		# grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ if( length($5)==0 || index($5,"-")==1 ) {gsub(/-/,"",$5); print $0; }}' >> $tsv.0.csv
		# grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ if(index($5,"-") > 1 || length($5)==16 ) { gsub(/-/,"",$5); print $0; } }' >> $tsv.16.csv
		grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ if( length($5)==0 || index($5,"-")==1 ) {print $0; }}' >> $tsv.0.csv
		grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ if(index($5,"-") > 1 || length($5)==16 ) {print $0; } }' >> $tsv.16.csv
		
		grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ if( length($5)==6) {print $0; }}' >> $tsv.6.csv
		grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ if( length($5)==7) {print $0; }}' >> $tsv.7.csv
		grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ if( length($5)==8) {print $0; }}' | grep -v "HALOPLEX" >> $tsv.8.csv
		grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ if( length($5)==8) {print $0; }}' | grep "HALOPLEX" >> $tsv.8.haloplex.csv

		fastqcluster=$chopFastqinMillions
	fi

	lines_0_csv=`wc -l $tsv.0.csv | cut -f1 -d " "`
	lines_6_csv=`wc -l $tsv.6.csv | cut -f1 -d " "`
	lines_7_csv=`wc -l $tsv.7.csv | cut -f1 -d " "`
	lines_8_csv=`wc -l $tsv.8.csv | cut -f1 -d " "`
	lines_8_haloplex=`wc -l $tsv.8.haloplex.csv | cut -f1 -d " "`
	lines_16_csv=`wc -l $tsv.16.csv | cut -f1 -d " "`

	
# Tip /ifs/data/c2b2/ngs_lab/ngs/usr/CASAVA/bin/configureBclToFastq.pl --help 
# Tip --use-bases-mask is important for different indexing schemes or SE or PE.
# Tip --adapter-sequence /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/PipelineCASAVA/new_scripts/truseq-adapters.fasta
# Tip --adapter-sequence /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/PipelineCASAVA/new_scripts/nextera-adapters.fasta

## Fire CASAVA command for the lanes with single indexing scheme 

	if [[ $SEforce == "yes" ]]; then R1="r1" ; fi
	let threads=$nt


	if [[ $lines_0_csv > 1 ]];
	then
		if [[ $ends == "SE" ]];then
			if [[ $indexing == "NI" ]];then mask="Y*" ;
			elif [[ $indexing == "MI" ]];then mask="Y*,Y*,Y*" ; 
			else
				mask="Y*,Y*";
			fi
		else	
			if [[ $indexing == "NI" ]];then mask="Y*,Y*" ;
			elif [[ $indexing == "MI" ]];then mask="Y*,Y*,Y*,Y*" ; 
			else
				mask="Y*,Y*,Y*";
			fi
		fi
		outdir=$demultiplexout/$DIRNAME0
		outdir=$OutDir/$DIRNAME0
		csv=$tsv.0.csv
		adapter=$PIPEBASE/truseq-adapters.fasta


		ruby $PIPEBASE/NovoSeqSampleSheet.rb $csv $csv.new.csv

		touch $csv.new.csv

		cmd="nohup $bcl2fastqBin --runfolder-dir $RunDir --output-dir $outdir --barcode-mismatches $barcodemismatch -r 4  -p $nt -w 4 --minimum-trimmed-read-length 10 --mask-short-adapter-reads 10 --sample-sheet $csv.new.csv --use-bases-mask $mask"
		echo $cmd
		$cmd

	fi
	
	if [[ $lines_6_csv > 1 ]];
	then
		if [[ $ends == "SE" ]];then
			if [[ $indexing == "MI" ]]
			then 
				mask="Y*,I6n*,n*"
			else
				mask="Y*,I6n*"
			fi
		else	
			if [[ $indexing == "MI" ]]
			then
				mask="Y*,I6n*,n*,Y*"
			else
				mask="Y*,I6n*,Y*"
			fi
		fi
		outdir=$demultiplexout/$DIRNAME6
		outdir=$OutDir/$DIRNAME6
		csv=$tsv.6.csv
		adapter=$PIPEBASE/truseq-adapters.fasta

		ruby $PIPEBASE/NovoSeqSampleSheet.rb $csv $csv.new.csv
		touch $csv.new.csv
                cmd="nohup $bcl2fastqBin --runfolder-dir $RunDir --output-dir $outdir --barcode-mismatches $barcodemismatch -r 4  -p $nt -w 4 --minimum-trimmed-read-length 10 --mask-short-adapter-reads 10 --sample-sheet $csv.new.csv --use-bases-mask $mask"
		echo $cmd
		$cmd

	fi
	if [[ $lines_7_csv > 1 ]];
	then
		if [[ $ends == "SE" ]];then
			if [[ $indexing == "MI" ]]
			then
				mask="Y*,I7n*,n*"
			else
				mask="Y*,I7n*"
			fi
		else	
			if [[ $indexing == "MI" ]]
			then
				mask="Y*,I7n*,n*,Y*"
			else
				mask="Y*,I7n*,Y*"
			fi
		fi
		outdir=$demultiplexout/$DIRNAME7
		outdir=$OutDir/$DIRNAME7
		csv=$tsv.7.csv
		adapter=$PIPEBASE/nextera-adapters.fasta

		ruby $PIPEBASE/NovoSeqSampleSheet.rb $csv $csv.new.csv

		touch $csv.new.csv
                cmd="nohup $bcl2fastqBin --runfolder-dir $RunDir --output-dir $outdir --barcode-mismatches $barcodemismatch -r 4  -p $nt -w 4 --minimum-trimmed-read-length 10 --mask-short-adapter-reads 10 --sample-sheet $csv.new.csv --use-bases-mask $mask"
		echo $cmd
		$cmd

	fi
	if [[ $lines_8_csv > 1 ]];
	then
		if [[ $ends == "SE" ]];then
			if [[ $indexing == "MI" ]]
			then
				mask="Y*,I8n*,n*"
			else
				mask="Y*,I8n*"
			fi
		else	
			if [[ $indexing == "MI" ]]
			then
				mask="Y*,I8n*,n*,Y*"
			else
				mask="Y*,I8n*,Y*"
			fi
		fi
		outdir=$demultiplexout/$DIRNAME8
		outdir=$OutDir/$DIRNAME8
		csv=$tsv.8.csv
		# adapter=$PIPEBASE/nextera-adapters.fasta
		adapter=$PIPEBASE/truseq-adapters.fasta
		
		ruby $PIPEBASE/NovoSeqSampleSheet.rb $csv $csv.new.csv

		touch $csv.new.csv
                cmd="nohup $bcl2fastqBin --runfolder-dir $RunDir --output-dir $outdir --barcode-mismatches $barcodemismatch -r 4  -p $nt -w 4 --minimum-trimmed-read-length 10 --mask-short-adapter-reads 10 --sample-sheet $csv.new.csv --use-bases-mask $mask"
		echo $cmd
		$cmd

	fi

	if [[ $lines_8_haloplex > 1 ]];
	then
		if [[ $ends == "SE" ]];then
			if [[ $indexing == "MI" ]]
			then
				mask="Y*,I8n*,n*"
			else
				mask="Y*,I8n*"
			fi
		else	
			if [[ $indexing == "MI" ]]
			then
				mask="Y*,I8n*,n*,Y*"
			else
				mask="Y*,I8n*,Y*"
			fi
		fi
		outdir=$demultiplexout/$DIRNAME8HALOPLEX
		outdir=$OutDir/$DIRNAME8HALOPLEX
		csv=$tsv.8.haloplex.csv
		adapter=$PIPEBASE/haloplex-adapters.fasta

		ruby $PIPEBASE/NovoSeqSampleSheet.rb $csv $csv.new.csv

		touch $csv.new.csv
                cmd="nohup $bcl2fastqBin --runfolder-dir $RunDir --output-dir $outdir --barcode-mismatches $barcodemismatch -r 4  -p $nt -w 4 --minimum-trimmed-read-length 10 --mask-short-adapter-reads 10 --sample-sheet $csv.new.ccsv --use-bases-mask $mask"
		echo $cmd
		$cmd

	fi

	if [[ $lines_16_csv > 1 ]];
	then
		if [[ $ends == "SE" ]];then
			mask="Y*,I8,I8"
		else	
			mask="Y*,I8n*,I8n*,Y*"
			
		fi
		outdir=$demultiplexout/$DIRNAME16
		outdir=$OutDir/$DIRNAME16
		csv=$tsv.16.csv
		adapter=$PIPEBASE/nextera-adapters.fasta

		ruby $PIPEBASE/NovoSeqSampleSheet.rb $csv $csv.new.csv
		touch $csv.new.csv
                cmd="nohup $bcl2fastqBin --runfolder-dir $RunDir --output-dir $outdir --barcode-mismatches $barcodemismatch -r 4  -p $nt -w 4 --minimum-trimmed-read-length 10 --mask-short-adapter-reads 10 --sample-sheet $csv.new.csv --use-bases-mask $mask"
		echo $cmd
		$cmd

	fi
}

mask="crap"
if [[ $numReads == 1 ]]; then
#mask="Y*"#SE - NI
    run_regular_casava $mask "$demultiplexout/$runName.csv" "NI" "SE"

elif [[ $numReads == 2 ]]; then
    if [[ "$lastread" == *"IsIndexedRead=\"Y\""* ]]; then
#mask="Y*,I6n"#SE - SI
	run_regular_casava $mask "$demultiplexout/$runName.csv" "SI" "SE"
	elif [[ "$lastread" == *"IsIndexedRead=\"N\""* ]];then
#mask="Y*,Y*"#PE - NI
	run_regular_casava $mask "$demultiplexout/$runName.csv" "NI" "PE"
	fi

elif [[ $numReads == 3 ]]; then
    if [[ "$lastread" == *"IsIndexedRead=\"Y\""* ]];then
	echo "run mixed- se"#SE - SI + DI
	run_regular_casava $mask "$demultiplexout/$runName.csv" "MI" "SE"
	elif [[ "$lastread" == *"IsIndexedRead=\"N\""* ]];then
#mask="Y*,I6n,Y*"#PE - SI
	run_regular_casava $mask "$demultiplexout/$runName.csv" "SI" "PE"
	fi

elif [[ $numReads == 4 ]]; then
    echo "run mixed-de" #PE - SI + DI
    run_regular_casava $mask "$demultiplexout/$runName.csv" "MI" "PE"
else
    echo "Something wrong with RunInfo.xml"
    exit
fi


cd $absOUT
touch $absOUT/BclToFastq_complete.txt
echo -e "bcl to fastq done"

###########################################################

echo -e "conversion done, Force=$force , SEforce=$SEforce" > $absOUT/BclToFastq_complete.txt
echo -e "$RunDir\t$outdir" >> $StatusDir/bcl2fastq.complete

###################################################################################


if [[ $SEforce == "yes" ]];then	exit ; fi ##Avoid trigering downstream if forcing SE on PE to get read numbers only
if [[ $force == "yes" ]]; then exit ; fi

popd

runID=`basename $OutDir`
chmod -R 770 $OutDir
echo -e "post demutilplex starts"

#for x in 0 6 7 8 16 8HALOPLEX    ## DIRNAME
#do
#    dir=$OutDir/Indexing$x/
#    if [ -e $dir ]
#    then
#	cd $dir
#	cmd="sh $PIPEBASE/post_demux_Hiseq4000.sh $dir $runID $setting"
#	echo $cmd
#	$cmd
#    fi
#done



echo "End $0:`date`"


