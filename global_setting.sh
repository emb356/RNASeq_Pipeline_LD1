#environment variables for ld1 bulk RNA-seq pipeline
export NovaSeqRuns=/home/local/ARCS/ngs/isilon/Illumina/
export NextSeqRuns=/home/local/ARCS/ngs/isilon/NextSeq

#test cronjob on ld1
#export StatusDir=/home/local/ARCS/ngs/Pipeline_ld1/Status
#export FastqDir=/home/local/ARCS/ngs/Pipeline_ld1/FastqCASAVA/
#export SampleSheets=/home/local/ARCS/ngs/bin/ngs_pipeline/restartOnNgsLd1/sampleSheets

export FastqDir=/home/local/ARCS/ngs/isilon/FastqCASAVA_ngsld1
export StatusDir=/home/local/ARCS/ngs/mnt/data/status/
export SampleSheets=/home/local/ARCS/ngs/mnt/data/status/sampleSheet
export RUBY18=/home/local/ARCS/ngs/mnt/data/usr/local/bin/ruby

export RNABASE=/home/local/ARCS/ngs/Pipeline_ld1/RNASEQ
export NGSSHELL=/home/local/ARCS/ngs/mnt/data/code/shell/
export UTILS=/home/local/ARCS/ngs/mnt/data/code/NGS/utils/
export Projects=/home/local/ARCS/ngs/Pipeline_ld1/Projects
#export PIPEBASE=/home/local/ARCS/ngs/Pipeline_ld1/scripts
export PIPEBASE=/home/local/ARCS/ngs/mnt/data/code/NGS/Production/CASAVA-Pipeline/
#export NextSeqOut=/home/local/ARCS/ngs/isilon/FastqNextSeq/
export bcl2fastqBin="bcl2fastq"
