Files included for each project (does not apply for user-mapped genomes and self-prepped libraries): 

1. summary.csv 
   - contains basic statistical information
	Sample name
	Number of reads 
	Number of pseudoaligned reads
	Ratio of pseudoaligned/total reads (equivalent to mapping ratio)


2. kallisto_raw_output directory
   contains abundance files (tsv and h5) 
   - sampleName/abundance.tsv: plaintext file of the abundance estimates. Includes estimated counts, TPM, effective length for each transcript.
   - sampleName/abundance.h5: HDF5 binary file containing run info, abundance estimates, bootstrap estimates, and transcript length information. This file can be read in by sleuth.


3. transcript_level directory
   - est_counts_transcripts_kallisto.txt: estimated counts for each transcript
   - tpm_values_transcripts_kallisto.txt: transcripts per million (tpm) values for each transcript 
  
 
4. est_counts_genes_kallisto.txt and tpm_values_genes_kallisto.txt
   - estimated counts and tpm values for each gene (computed from tximport R package)


5. Zipped raw fastq
   Paired-end sequence: sampleName_R1_001.fastq.gz (forward strand) and sampleName_R2_001.fastq.gz (reverse strand)


6. Differentially expressed gene analysis (Sleuth or DESeq2) - Upon Request
   Sleuth Output:
   qval      		false discovery rate adjusted p-value, using Benjamini-Hochberg
   pval      		p-value of the Wald test
   b         		'beta' value (effect size). Technically a biased estimator of the fold change
   se_b      		standard error of the beta
   mean_obs  		mean of natural log counts of observations
   var_obs   		variance of observation
   tech_var  		technical variance of observation from the bootstraps 
   sigma_sq  		raw estimator of the variance once the technical variance has been removed
   smooth_sigma_sq 	smooth regression fit for the shrinkage estimation
   final_sigma_sq 	max(sigma_sq, smooth_sigma_sq); used for covariance estimation of beta
   DESeq2 Output:
   baseMean             intermediate mean of normalized counts for all samples
   log2FoldChange       log2 fold change
   lfcSE                standard error
   stat                 Wald statistic
   pvalue               Wald test p-value
   padj                 BH adjusted p-values



* Some samples might be sequenced in multiple lanes, gene/isoform abundance calculation is based on the merged mapped reads.

* You could use wget to download all the files at once if you are familiar with Unix. Like this: wget -r webLink


Standard RNA-seq pipeline description: 

For STRDPOLYA: We use poly-A pull-down to enrich mRNAs from total RNA samples,
then proceed with library construction using Illumina TruSeq chemistry.
For RIBOZERO: We use ribosomal depletion to remove rRNAs from total RNA
samples, then proceed with library construction using Illumina TruSeq
chemistry.
For CLONTECH: We use the Clontech Ultra Low v4 kit for cDNA synthesis followed by NexteraXT.

Libraries are then sequenced using Illumina NovaSeq 6000 at Columbia Genome
Center. We multiplex samples in each lane, which yields targeted number of
paired-end 100bp reads for each sample.

We use RTA (Illumina) for base calling and bcl2fastq2 (version 2.19) for converting BCL to fastq format, coupled with adaptor trimming. We perform a pseudoalignment to a kallisto index created from transcriptomes (Human: GRCh38; Mouse: GRCm38) using kallisto (0.44.0). We test for differentially expressed genes under various conditions using Sleuth or DESeq2, R packages designed to test differential expression between two experimental groups from RNA-seq counts data.


GRANT ACKNOWLEDGEMENT
We are a Shared Resource of the HICCC, and we are required to report on publications that resulted from direct cost support from the P30 Cancer Center Support Grant. Please cite the P30 as follows: “This research was funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696 and used the Genomics and High Throughput Screening Shared Resource.” Please add the paper to your bibliography in My NCBI and link the paper with P30CA013696.
