#!/bin/bash
#$ -cwd
# file name: kallistoMapping.sh
# maps fastqs using kallisto
# executed in RNA_pipeline_kallisto.sh
 
USAGE="Usage: kallistoMapping.sh -s sample -g genome -l libtype -o output -f fastq1 -p fastq2"

while getopts "s:i:l:o:f:p:h" opt;
    do
    case "$opt" in
	s) sample="$OPTARG";;
        i) index="$OPTARG";;
        l) libtype="$OPTARG";;
        o) output="$OPTARG";;
        f) fastq1=("$OPTARG");;
        p) fastq2=("$OPTARG");;
	h) echo $USAGE; exit 1
    esac
done 
shift "$((OPTIND-1))"

# make output directory if it doesnt exist
if [[ ! -d $output ]]
then
    mkdir -p $output
fi

# check if SE (deprecated)
if [[ $fastq2 == "NA" ]]
   then
   if [[ $libtype == "N" ]]
      then 
      cmd="kallisto quant -i $index -o $output -t 4 -l 200 -s 20 --single $fastq1"
   else
      cmd="kallisto quant -i $index -o $output -t 4 -l 200 -s 20 --single $fastq1 --rf-stranded"
   fi
else   
   if [[ $libtype == "N" ]]
      then 
      cmd="kallisto quant -i $index -o $output -t 4 $fastq1 $fastq2"
   else
      cmd="kallisto quant -i $index -o $output -t 4 $fastq1 $fastq2 --rf-stranded"
   fi 
fi 

echo $fastq1
echo $fastq2
echo $cmd

$cmd

echo "Kallisto is done"
echo `date +"%H:%M:%S"`

# create a file indicating kallisto is done, used for executing kallistoRelease.sh from post_demux_Novaseq6000.sh
touch 'logs/'$sample'_kallisto_done.txt'
