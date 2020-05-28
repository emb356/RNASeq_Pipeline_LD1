#!/bin/bash -l
#$ -cwd
# file name: RNA_pipeline_kallisto.sh
# executed in post_demux_Novaseq6000.sh 
# or can be called from trigRNA_pipeline_kallisto.sh for manually starting RNA-seq analysis
# three parameters are needed: genome, fastq1, fastq2 

genome=$1 

# fastq files
f1=$2 # R1
f2=$3 # R2


# path to code and reference transcripts
RNABASE="/home/local/ARCS/ngs/Pipeline_ld1/scripts" 
GENOMEBASE="/home/local/ARCS/ngs/Pipeline_ld1/scripts/genome_settings"

# if Ribozero, use transcriptome index with both exome and ncrna
RIBO=`echo $f1 | grep "RIBOZERO" | wc -c`
if [[ $RIBO -gt 0 ]]; then
    setting_genome="$GENOMEBASE/global_setting_"$genome"ALL.sh"
    . $setting_genome
else
    setting_genome="$GENOMEBASE/global_setting_"$genome".sh"
    . $setting_genome
fi

# check if genome in global settings
if [[ ! -s $setting_genome ]]; then
    echo "No such genome avaliable"
    exit 1
fi

# make log and kallisto directories
if [ ! -d  logs ]; then mkdir -p logs; fi
if [ ! -d  $genome"_kallisto" ]; then mkdir -p $genome"_kallisto"; fi

# finding library type
# libType = N for unstranded or S for stranded
libType="S"
findType1=`echo $f1 | grep "STRDPOLYA\|RIBOZERO\|PLATESEQ" | wc -c`
findType2=`echo $f1 | grep "CLONTECH\|TRUSEQ\|NEXTERA" | wc -c`

if [[ $findType1 -gt 0 ]]; then
    libType="S"
elif [[ $findType2 -gt 0 ]]; then
    libType="N"
else
    echo "unable to determine the library type of $f1 - exitting!"
    exit 12
fi
echo $libType
echo $genome
echo $f1
echo $f2
echo $INDEX

f1_base=`basename $f1`
long_sample=${f1_base%%.fastq.gz*}
short_sample=${f1_base%%_*}

echo $long_sample
echo $short_sample

# SE (single end) mapping (deprecated)
if [[ $f2 == "" ]]; then
     sh $RNABASE/kallistoMapping.sh -s "$long_sample" -i "$INDEX" -l "$libType" -o "$genome""_kallisto/""$short_sample" -f "$f1" -p "NA" 
else
     sh $RNABASE/kallistoMapping.sh -s "$long_sample" -i "$INDEX" -l "$libType" -o "$genome""_kallisto/""$short_sample" -f "$f1" -p "$f2" 
fi
