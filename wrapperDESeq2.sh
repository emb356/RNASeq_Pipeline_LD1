#!/bin/bash
#$ -cwd
# filename: wrapperDESeq2.sh
# just a wrapper that sets up the anaconda R environment

# set the .Rprofile for the Rscripts
source activate HPC_Pipeline
cp -p /home/local/ARCS/ngs/Pipeline_ld1/scripts/.Rprofile .

Rscript DESeq2_kallisto.R

