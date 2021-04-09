#!/bin/env Rscript
# file name: manualDESeq2.R
# runs DESeq2 using counts files (not kallisto output)
# counts data must be in the form of tab-delimited "to_release/[sample]_counts.txt" files, first row is genes second row is sample counts
# change end of file by sample name accordingly

library( DESeq2 )
library(dplyr)
library(ggplot2)
runDESeq2 <- function( sampleList ) {
WD=getwd()

n <- length(sampleList[[1]])
m <- length(sampleList[[2]])
print(sampleList)
samples <- unlist(sampleList, recursive=FALSE)
print(samples)
# getting counts files based on pattern
txt_files <- paste("to_release/",samples,"_counts.txt",sep="")
txt_files
names(txt_files) <-samples
names(txt_files)

# run DESeq2
sampleTable <- data.frame(condition = factor(rep(c("ctrl", "exp"), c(n,m))))
# modify for batch effect across projects, and change to dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~batch + condition)
# sampleTable$batch <- c("1","1","1","1","1","2","2","2","2","2","1","1","1","1","1","2","2","2","2","2","2")

cts <- Reduce(function(...) merge(..., by=1), lapply( txt_files, read.delim, stringsAsFactors=F ) )
head(cts)
rownames( cts ) <- cts[,1]

cts <- cts[,2:dim(cts)[2]]
cts <- round(cts)
head(cts)
colnames(cts) <- samples
rownames(sampleTable) <- colnames(cts)
sampleTable
dds <- DESeqDataSetFromMatrix( countData = cts,
                               colData = sampleTable,
                               design = ~ condition )
dds <- DESeq( dds )
res <- results(dds)
res <- res[order(res$padj),]
res_df <- as.data.frame(res)
top_genes <- rownames(res_df)[1:50]
print(sampleTable)
print(top_genes)
# write results to to_release
release_dir = file.path(WD, "to_release")
setwd(release_dir)


fnameStem <- paste( "DESeq2",
                     paste( sapply(samples, paste, collapse="-"),
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
print("test")
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsdData <- assay(vsd)
vsd_df <- as.data.frame(vsdData)
vsd_df$genes <- rownames(vsd_df)
vsd_sub <- filter(vsd_df, vsd_df$genes %in% top_genes)
rownames(vsd_sub) <- vsd_sub$genes
vsd_sub <- within(vsd_sub, rm(genes))
print("test")
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
hmcol <- colorRampPalette( brewer.pal( 9, "GnBu" ))( 200 )
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

A <- qw(KC005,KC006)
B <- qw(KC003,KC004)
C <- qw(KC011,KC012)
D <- qw(KC009,KC010)

runDESeq2( sampleList=list( A, B ))
runDESeq2( sampleList=list( C, D ))

print("DESeq2 is done.")


