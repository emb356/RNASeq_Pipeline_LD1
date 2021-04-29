#!/bin/env Rscript
# file name: DESeq2_kallisto.R
# runs DESeq2 on kallisto using tximport
# copy this script to [organism]_kallisto
# execute using wrapperDESeq2.sh
# enter sample names at bottom

library(tximport)
library(biomaRt)
library(DESeq2)
library(ggplot2)
library(dplyr)


runDESeq2kallisto <- function( sampleList ) {
WD=getwd()

n <- length(sampleList[[1]])
m <- length(sampleList[[2]])

samples <- unlist(sampleList, recursive=FALSE)
tsv_files <- file.path("to_release/kallisto_raw_output", samples, "abundance.tsv")
names(tsv_files) <-samples
genome <- sub(".*/(.*)_kallisto$","\\1", WD) 
genome
tsv_files
names(tsv_files)

# check for valid genome name
if (genome == "human" | genome == "humanALL"){
  ensembl_dataset = "hsapiens_gene_ensembl"
} else if (genome == "mouse" | genome == "mouseALL"){
  ensembl_dataset = "mmusculus_gene_ensembl"
} else if (genome == "drosophi"){
  ensembl_dataset = "dmelanogaster_gene_ensembl"
} else if (genome == "drosophi"){
  ensembl_dataset = "dmelanogaster_gene_ensembl"
} else if (genome == "horse"){
  ensembl_dataset = "ecaballus_gene_ensembl"
  transcript_id = "refseq_mrna"
} else if (genome == "pig"){
  ensembl_dataset = "sscrofa_gene_ensembl"
  transcript_id = "refseq_mrna"
} else if (genome == "monkey"){
  ensembl_dataset = "csabaeus_gene_ensembl"
} else if (genome == "rat"){
  ensembl_dataset = "rnorvegicus_gene_ensembl"
} else {
  stop("Genome has been entered incorrectly. Check the script for valid genomes", call.=FALSE)
}
print(ensembl_dataset)
tx2g <- function(){
  #mart <- useMart("ensembl", ensembl_dataset, host = "useast.ensembl.org")
  mart <- biomaRt::useMart(biomart = "ensembl", host = "www.ensembl.org", ensemblRedirect = FALSE, dataset = ensembl_dataset)
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                         "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                         ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  return(t2g)
}
t2g <- tx2g()

txi.kallisto <- tximport(tsv_files, type = "kallisto", tx2gene = t2g[,c(1,3)], ignoreTxVersion = TRUE)

# run DESeq2
sampleTable <- data.frame(condition = factor(rep(c("ctrl", "exp"), c(n,m))))
# modify for batch effect acrosss projects, and change to dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~batch + condition)
# sampleTable$batch <- c("1","1","1","1","1","2","2","2","2","2","1","1","1","1","1","2","2","2","2","2","2")

rownames(sampleTable) <- colnames(txi.kallisto$counts)
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
res_df <- as.data.frame(res)
top_genes <- rownames(res_df)[1:50]
print(sampleTable)

# write results to to_release
release_dir = file.path(WD, "to_release")
setwd(release_dir)

fnameStem <- paste( "DESeq2_kallisto",
                    paste( sapply(sampleList[[1]], paste, collapse="-"),
                           collapse = "_" ), "vs",
                    paste( sapply(sampleList[[2]], paste, collapse="-"),
                           collapse = "_" ),
                    sep="_" )

# unix limits the length of filenames to 256 chars, so just to be careful, if the length is 200 chars, set it to something simple
if ( nchar( fnameStem ) > 200 )
     fnameStem <- paste( "DESeq2_kallisto",
                    paste( sapply(length(sampleList[[1]]), paste, collapse="-"),
                           collapse = "_" ), "samples_vs",
                    paste( sapply(length(sampleList[[2]]), paste, collapse="-"),
                           collapse = "_" ),
                    sep="_" )

fname <- paste( fnameStem,
                "DEG_list.csv",
                sep="-" )

write.csv(res , file=fname , quote=F)

fname <- paste( fnameStem,
                "DEG_Plots.pdf",
                sep="-" )

# vsd transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsdData <- assay(vsd)
vsd_df <- as.data.frame(vsdData)
vsd_df$genes <- rownames(vsd_df)
vsd_sub <- filter(vsd_df, vsd_df$genes %in% top_genes)
rownames(vsd_sub) <- vsd_sub$genes
vsd_sub <- within(vsd_sub, rm(genes))

fname <- paste( fnameStem,
                "vsd.csv",
               sep="_" )

write.csv(vsdData, file=fname, quote=F)

# normalized counts
normalTable <- counts( dds, normalized=TRUE )
fname <- paste( fnameStem,
                "counts_table_normalized.txt",
                sep="_" )
write.table( normalTable,
             file=fname,
             sep="\t",
             quote=FALSE,
             na="0")


# start PDF
fname <- paste( fnameStem,
                "DEG_Plots.pdf",
                sep="-" )
pdf( file=fname )

plotMA <- plotMA(res, ylim=c(-2,2), main = "MA-Plot")
print(plotMA)

res_df <- filter(as.data.frame(res), padj > 0)
res_df <- transform(res_df, significance = ifelse(padj<0.1, "True", "False"))
plotVolcano <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color=significance))+
         geom_point()+
         ggtitle("Volcano Plot")+
         xlab("log2FoldChange")+
         ylab("-log10(padj)")
print(plotVolcano)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plotPCA <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_text(aes(label=rownames(sampleTable) )) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ggtitle("PCA after VST")+
  scale_x_continuous(expand = c(.1, .1)) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
print(plotPCA)

library( "RColorBrewer" )
library( "gplots" )
hmcol <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
plotHeatmap <- heatmap.2( as.matrix(vsd_sub),
           main="Top 50 DE Genes after VST (by Significance)",
           col = hmcol,
           scale="row",
           trace="none")
print(plotHeatmap)
dev.off()
setwd(WD)
}

qw <- function(...) {
  sapply(match.call()[-1], deparse)
}

# enter sample names
A <- qw(CT097,CT098,CT100,CT101)
B <- qw(CX005,CX006,CX007)

runDESeq2kallisto( list( A, B ))


print("DESeq2 is done.")
