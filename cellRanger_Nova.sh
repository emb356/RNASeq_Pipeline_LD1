#!/bin/bash -l
#$ -cwd
# file name: cellRanger_Nova.sh
# demultiplexes single cell data from Novaseq 6000

rundir=$1
outdir=$2
samplesheet=$3

USAGE="Usage: $0 -i [path to input run directory]  -o [path to output directory] -s [path to sample sheet]"

while getopts i:o:s:h opt
  do
  case "$opt" in
      i) RunDir="$OPTARG";;
      o) OutDir="$OPTARG";;
      s) SampSheet="$OPTARG";;
      h) echo $USAGE; exit 1
  esac
done

if [[ $RunDir == "" || $OutDir == ""  || $SampSheet == "" ]]
    then
    echo $USAGE
    exit 1
fi

# asssume only one index type per run for single cell
index8=`grep "SI-GA\|SI-NA" $SampSheet | wc -l`
ignore=""
# check if single indexes. If they are, ignore the dual indexes on the run (10bp)
if [ $index8 -gt 0 ]; then
#   mask="--use-bases-mask=Y101,I8N2,N10,Y101"
#   mask="--use-bases-mask=Y26n*,I10,I10,Y90n*"
   ignore="--filter-single-index"
else
#  mask="--use-bases-mask=Y28n*,I10,I10,Y90n*"
  ignore="--filter-dual-index"
fi
   
cmd="cellranger mkfastq --run=$RunDir --output-dir=$OutDir --samplesheet=$SampSheet --localcores=8 --localmem=34 $ignore --barcode-mismatches=0"
#cmd="cellranger mkfastq --run=$RunDir --output-dir=$OutDir --samplesheet=$SampSheet --localcores=8 --localmem=34 --barcode-mismatches=0"

#cmd="cellranger mkfastq --run=$RunDir --output-dir=$OutDir --samplesheet=$SampSheet --localcores=8 --localmem=34"
#cmd="cellranger mkfastq --run=$RunDir --output-dir=$OutDir --csv=$SampSheet --localcores=8 --localmem=34 --barcode-mismatches=0"

#cmd="cellranger mkfastq --run=$RunDir --output-dir=$OutDir --samplesheet=$SampSheet --localcores=8 --localmem=34 --barcode-mismatches=0 --use-bases-mask=Y27n*,I10,I10,Y91n* --ignore-dual-index"
echo $cmd
$cmd
