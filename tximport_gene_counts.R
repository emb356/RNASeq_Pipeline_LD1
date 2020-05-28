#!/bin/env Rscript
# file name: tximport_gene_counts.R
# creates gene counts table from kallisto files
# executed in kallistoRelease.sh

library(tximport)
WD=getwd()
dir = file.path(WD,"kallisto_raw_output")
dir
tsv_files <- list.files(path=dir, pattern = "*.tsv", full.names=TRUE)
names(tsv_files) <- sub('\\..*$', '', basename(tsv_files)) 
genome <- sub(".*/(.*)_.*/to_release$","\\1", WD) 
genome

# standard ensembl biomaRt attributes
transcript_id <- "ensembl_transcript_id"
gene_id <- "ensembl_gene_id"
gene_name <- "external_gene_name"

# load libraries
library(biomaRt,lib.loc="/home/local/ARCS/ngs/miniconda2/envs/HPC_Pipeline/lib/R/library")
library(dplyr,lib.loc="/home/local/ARCS/ngs/miniconda2/envs/HPC_Pipeline/lib/R/library")
library(curl,lib.loc="/home/local/ARCS/ngs/miniconda2/envs/HPC_Pipeline/lib/R/library")
library(rhdf5,lib.loc="/home/local/ARCS/ngs/miniconda2/envs/HPC_Pipeline/lib/R/library")
 
# check for valid genome name
if (genome == "human" | genome == "humanALL" | genome == "humanLNCRNA"){
  ensembl_dataset = "hsapiens_gene_ensembl"
} else if (genome == "mouse" | genome == "mouseALL" | genome == "mouseLNCRNA"){
  ensembl_dataset = "mmusculus_gene_ensembl"
} else if (genome == "drosophi"){
  ensembl_dataset = "dmelanogaster_gene_ensembl"
} else if (genome == "drosophi"){
  ensembl_dataset = "dmelanogaster_gene_ensembl"
} else if (genome == "horse"){
  ensembl_dataset = "ecaballus_gene_ensembl"
  transcript_id = "refseq_mrna"
} else if (genome == "monkey"){
  ensembl_dataset = "csabaeus_gene_ensembl"
} else if (genome == "rat"){
  ensembl_dataset = "rnorvegicus_gene_ensembl"
} else if (genome == "yeast"){
  ensembl_dataset = "scerevisiae_gene_ensembl"
  gene_name = "ensembl_gene_id"
} else {
  stop("Genome has been entered incorrectly. Check the script for valid genomes", call.=FALSE)
}


tx2g <- function(){
#  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = ensembl_dataset)
  mart <- biomaRt::useMart("ensembl", ensembl_dataset, host = "uswest.ensembl.org", ensemblRedirect = FALSE)
#  mart <- biomaRt::useMart("ensembl", ensembl_dataset, host = "ensembl.org", ensemblRedirect = FALSE)
#  mart <- biomaRt::useMart(biomart = "ensembl", host = "www.ensembl.org", ensemblRedirect = FALSE, dataset = ensembl_dataset)
  t2g <- biomaRt::getBM(attributes = c(transcript_id, gene_name), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = transcript_id, ext_gene = gene_name)

  return(t2g)
}
t2g <- tx2g()

txi.kallisto <- tximport(tsv_files, type = "kallisto", tx2gene = t2g[,c(1,2)], ignoreTxVersion = TRUE)

gene_counts = round(txi.kallisto$counts)
gene_tpm = txi.kallisto$abundance
write.table( gene_counts, file="est_counts_genes_kallisto.txt", sep="\t", row.names=T, quote=F)
write.table( gene_tpm, file="tpm_values_genes_kallisto.txt", sep="\t", row.names=T, quote=F)
