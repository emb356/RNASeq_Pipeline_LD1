# filename: global_setting.sh
# environment variables for ld1 bulk RNA-seq pipeline

export NovaSeqRuns=/home/local/ARCS/ngs/isilon/Illumina/
export NextSeqRuns=/home/local/ARCS/ngs/isilon/NextSeq

#test cronjob on ld1
#export StatusDir=/home/local/ARCS/ngs/Pipeline_ld1/Status
#export FastqDir=/home/local/ARCS/ngs/Pipeline_ld1/FastqCASAVA/
#export SampleSheets=/home/local/ARCS/ngs/bin/ngs_pipeline/restartOnNgsLd1/sampleSheets

export FastqDir=/home/local/ARCS/ngs/isilon/FastqCASAVA_ngsld1
export StatusDir=/home/local/ARCS/ngs/mnt/data/status/
export SampleSheets=/home/local/ARCS/ngs/mnt/data/status/sampleSheet
#export SampleSheets=/home/local/ARCS/ngs/mnt/data/status/dummySampleDir
export RUBY18=/home/local/ARCS/ngs/mnt/data/usr/local/bin/ruby

export NGSSHELL=/home/local/ARCS/ngs/mnt/data/code/shell/
export UTILS=/home/local/ARCS/ngs/mnt/data/code/NGS/utils/
export Projects=/home/local/ARCS/ngs/Pipeline_ld1/Projects
export PIPEBASE=/home/local/ARCS/ngs/mnt/data/code/NGS/Production/CASAVA-Pipeline/

export bcl2fastqBin="bcl2fastq"
export KALLISTO=/home/local/ARCS/ngs/mnt/data/code/NGS/RNA_seq/KALLISTO
export RNACode=/home/local/ARCS/ngs/mnt/data/code/NGS/RNA_seq/CASAVA
export LD1BASE=/home/local/ARCS/ngs/Pipeline_ld1/scripts
