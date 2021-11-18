#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(tximport)
library(DESeq2)
library(rtracklayer)
library(data.table)
library(dplyr)
library(tidyr)
dir.exists("output/011_deseq2")
#files and directories
salmon_dir <- snakemake@input[["salmon_dir"]]
samples_file <- snakemake@input[["samples_file"]]
gene_ids <- snakemake@input[["gene_ids"]]
emapper_annotations_file <- snakemake@input[["emapper_annotations_file"]]
emapper_insecta_annotations_file <- snakemake@input[["emapper_insecta_annotations_file"]]

dds_file <- snakemake@output[["dds_file"]]
deseq2_result_file <- snakemake@output[["deseq2_results_file"]]
deseq2_count_file <- snakemake@output[["deseq2_counts_file"]]
deseq2_result_with_annotations_file <- snakemake@output[["deseq2_results_with_annotations_file"]]

samples <- read.table(samples_file, header = TRUE)
files <- file.path(salmon_dir, "quant", samples$sample, "quant.sf")
gene_ids <- read.delim(gene_ids, header = FALSE)

# extract a data.frame of tx to gene
setnames(gene_ids, old = c('V1','V2'), new = c('GENEID', 'TXNAME'))
setcolorder(gene_ids, c('TXNAME', 'GENEID'))
tx2gene <- gene_ids

names(files) <- paste0(samples$sample)

#import quant files
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

#generate DESeq2 object
dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~group)
saveRDS(dds, dds_file)

# run DESeq2
dds <- DESeq(dds)
res <- results(dds)
res <- subset(res, padj < 0.01)
res_ord <- as.data.table(res[order(res$padj),])
counts <- as.data.frame(counts(dds, normalized=TRUE))

# write results table
fwrite(res_ord, deseq2_results_file)
fwrite(counts, deseq2_counts_file)

#join emapper annotations with DEseq2 results
# read results and add column names
my_results <- cbind(rownames(res), data.frame(res, row.names=NULL))
setnames(my_results, old = c('rownames(res)'), new = c('gene'))
my_results_dt <- as.data.table(my_results)
my_results_wtrans_ids <- merge(my_results_dt, tx2gene, by.x = "gene", by.y = "GENEID", all = TRUE)
setnames(my_results_wtrans_ids, old = c('TXNAME'), new = c('transcript'))

#read and prep. emapper insecta annotations
emapper_insecta_annot <- read.delim(emapper_insecta_annotations_file, header=FALSE, comment.char="#")
sep_insecta_emapper_annot <- separate(emapper_insecta_annot, "V1", c('transcript', 'V1'), sep = ".p")
prepd_insecta_annot <- subset(sep_insecta_emapper_annot, , -c(V1, V2, V3, V4, V5, V6, V7, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20))
prepd_insecta_annot <- prepd_insecta_annot[, c(1, 3, 2, 5, 4)]
colnames(prepd_insecta_annot) <- c('transcript', 'emapper_insecta_gene', 'emapper_insecta_description', 'emapper_insecta_PFAM', 'emapper_insecta_go')

# read and prep. emapper annotations
emapper_annot <- read.delim(emapper_annotations_file, header=FALSE, comment.char="#")
sep_emapper_annot <- separate(emapper_annot, "V1", c('transcript', 'V1'), sep = ".p")
prepd_annot <- subset(sep_emapper_annot, , -c(V1, V2, V3, V4, V5, V7, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20))
prepd_annot <- prepd_annot[, c(1, 2, 4, 3, 6, 5)]
colnames(prepd_annot) <- c('transcript', 'emapper_taxa', 'emapper_gene', 'emapper_description', 'emapper_PFAM', 'emapper_go')

# merge results with emapper annotations
merged_insecta <- merge(my_results_wtrans_ids, prepd_insecta_annot, by.x = "transcript", by.y = "transcript", all.x = TRUE)
merged <- merge(merged_insecta, prepd_annot, by.x = "transcript", by.y = "transcript", all.x = TRUE)
merged_ord <-as.data.table(merged[order(merged$padj),])

# write results with annotations
fwrite(merged_ord, deseq2_results_with_annotations_file)
