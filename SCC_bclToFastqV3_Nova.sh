#!/bin/bash -l
#$ -cwd
# file name: SCC_bclToFastqV3_Nova.sh
# demultiplexes run data from NovaSeq6000

uname -a

source $HOME/.bashrc
echo "Start $0:`date`"

# default global settings
setting="/home/local/ARCS/ngs/Pipeline_ld1/scripts/global_setting.sh"

USAGE="Usage: $0 -i RunDir  -o OutDir  -s SampleSheet [-n num_threads] [-m num of mismatches to be allowed in barcodes] [-c #reads to be placed in each fastq file] [-F flag to force demultiplexing of run w/o post-processing] "

# default params
nt=16
force="no"
barcodemismatch=1
chopFastqinMillions="250000000"

while getopts i:o:b:s:n:m:c:hFS opt
  do
  case "$opt" in
      i) RunDir="$OPTARG";;
      o) OutDir="$OPTARG";;
      b) bcl2fastqBin="$OPTARG";;
      s) SampleSheet="$OPTARG";;
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
absIN=`readlink -f $RunDir`  
absOUT=`readlink -f $OutDir`
INTENSITY=$absIN/Data/Intensities/
BCALL=$absIN/Data/Intensities/BaseCalls/
runName=`basename $absIN`



if [[ $force == "yes" ]];then absOUT=$absOUT"_FORCE" ; fi	# In case bcl2fastq needs to be re-done w/o post-processing, use the -F flag to Force

echo "outdir: $OutDir     $absOUT "
if [[ ! -e $absOUT ]]; then mkdir $absOUT ; fi
chmod 770 $absOUT
pushd $absOUT

#### Demultiplex ##################################################################################

numlanes=`grep "<Lane>" $absIN/RunInfo.xml | wc -l`   
numReads=`grep "Read Num" $absIN/RunInfo.xml | wc -l`
lastread=`grep "Read Num" $absIN/RunInfo.xml | tail -1 `

sampleSheet=$SampleSheet
echo "Using Sample Sheet : $sampleSheet"
cat $sampleSheet
if [[ ! -s $sampleSheet ]]; then  echo "SampleSheet is missing. Failed to start demultiplexing." ;  exit ; fi
demultiplexout=$absOUT/demultiplex
if [[ ! -e $demultiplexout ]]; then mkdir $demultiplexout ; fi

# remove SCC projects from samplesheet
awk '!/_SCC_/' $sampleSheet > $demultiplexout/$runName.csv
sampleSheet=$demultiplexout/$runName.csv
#cp $sampleSheet $demultiplexout/$runName.csv
cd $demultiplexout
DIRNAME0="Indexing0"
DIRNAME6="Indexing6" 	
DIRNAME7="Indexing7"
DIRNAME8="Indexing8"
DIRNAME8HALOPLEX="Indexing8HALOPLEX"
DIRNAME16="Indexing16"
DIRNAME20="Indexing20"

lines_20_csv=0
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
		grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ if( length($5)==0 || index($5,"-")==1 ) {print $0; }}' >> $tsv.0.csv
		grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ if(index($5,"-") > 1 && length($5)==17 ) {print $0; } }' >> $tsv.16.csv
		grep -v "FCID,Lane," $tsv | awk -F"," 'BEGIN{OFS=",";}{ if(index($5,"-") > 1 && length($5)==21 ) {print $0; } }' >> $tsv.20.csv
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
	lines_20_csv=`wc -l $tsv.20.csv | cut -f1 -d " "`

	
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
		outdir=$absOUT/$DIRNAME0
		csv=$tsv.0.csv

		ruby $LD1BASE/NovaSeqSampleSheet.rb $csv $csv.new.csv

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
		outdir=$absOUT/$DIRNAME6
		csv=$tsv.6.csv

		ruby $LD1BASE/NovaSeqSampleSheet.rb $csv $csv.new.csv
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
		outdir=$absOUT/$DIRNAME7
		csv=$tsv.7.csv

		ruby $LD1BASE/NovaSeqSampleSheet.rb $csv $csv.new.csv

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
		outdir=$absOUT/$DIRNAME8
		csv=$tsv.8.csv
		
		ruby $LD1BASE/NovaSeqSampleSheet.rb $csv $csv.new.csv

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
		outdir=$absOUT/$DIRNAME8HALOPLEX
		csv=$tsv.8.haloplex.csv

		ruby $LD1BASE/NovaSeqSampleSheet.rb $csv $csv.new.csv

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
		outdir=$absOUT/$DIRNAME16
		csv=$tsv.16.csv

		ruby $LD1BASE/NovaSeqSampleSheet.rb $csv $csv.new.csv
		touch $csv.new.csv
                cmd="nohup $bcl2fastqBin --runfolder-dir $RunDir --output-dir $outdir --barcode-mismatches $barcodemismatch -r 4  -p $nt -w 4 --minimum-trimmed-read-length 10 --mask-short-adapter-reads 10 --sample-sheet $csv.new.csv --use-bases-mask $mask"
		echo $cmd
		$cmd

	fi

	if [[ $lines_20_csv > 1 ]];
	then
		if [[ $ends == "SE" ]];then
			mask="Y*,I10,I10"
		else	
			mask="Y*,I10n*,I10n*,Y*"
			
		fi
		outdir=$absOUT/$DIRNAME20
		csv=$tsv.20.csv

		ruby $LD1BASE/NovaSeqSampleSheet.rb $csv $csv.new.csv
		touch $csv.new.csv
                cmd="nohup $bcl2fastqBin --runfolder-dir $RunDir --output-dir $outdir --barcode-mismatches $barcodemismatch -r 4  -p $nt -w 4 --minimum-trimmed-read-length 10 --mask-short-adapter-reads 10 --sample-sheet $csv.new.csv --use-bases-mask $mask"
		echo $cmd
		$cmd

	fi


}

mask="crap"
if [[ $numReads == 1 ]]; then
# mask="Y*" SE - NI
    run_regular_casava $mask "$demultiplexout/$runName.csv" "NI" "SE"

elif [[ $numReads == 2 ]]; then
    if [[ "$lastread" == *"IsIndexedRead=\"Y\""* ]]; then
# mask="Y*,I6n" SE - SI
	run_regular_casava $mask "$demultiplexout/$runName.csv" "SI" "SE"
	elif [[ "$lastread" == *"IsIndexedRead=\"N\""* ]];then
# mask="Y*,Y*" PE - NI
	run_regular_casava $mask "$demultiplexout/$runName.csv" "NI" "PE"
	fi

elif [[ $numReads == 3 ]]; then
    if [[ "$lastread" == *"IsIndexedRead=\"Y\""* ]];then
	echo "run mixed- se"#SE - SI + DI
	run_regular_casava $mask "$demultiplexout/$runName.csv" "MI" "SE"
	elif [[ "$lastread" == *"IsIndexedRead=\"N\""* ]];then
# mask="Y*,I6n,Y*" PE - SI
	run_regular_casava $mask "$demultiplexout/$runName.csv" "SI" "PE"
	fi

elif [[ $numReads == 4 ]]; then
    echo "run mixed-de" # PE - SI + DI
    run_regular_casava $mask "$demultiplexout/$runName.csv" "MI" "PE"
else
    echo "Something wrong with RunInfo.xml"
    exit
fi


cd $absOUT
touch $absOUT/BclToFastq_complete.txt
echo -e "bcl to fastq done"

###########################################################

echo -e "conversion done, Force=$force" > $absOUT/BclToFastq_complete.txt
echo -e "$RunDir\t$absOUT" >> $StatusDir/bcl2fastq.complete

###################################################################################


if [[ $force == "yes" ]]; then exit ; fi

popd

runID=`basename $absOUT`
chmod -R 770 $absOUT
echo -e "post demutilplex starts"

for x in 0 6 7 8 16 20 8HALOPLEX    # DIRNAME
do
    dir=$absOUT/Indexing$x/
    if [ -e $dir ]
    then
	cd $dir
	cmd="sh $LD1BASE/post_demux_Novaseq6000.sh $dir $runID $setting"
	echo $cmd
	$cmd
    fi
done



echo "End $0:`date`"


