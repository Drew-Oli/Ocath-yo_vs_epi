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

#files and directories
salmon_dir <- snakemake@input[["salmon_dir"]]
samples_file <- snakemake@input[["samples_file"]]
gff <- snakemake@input[["gff"]]
emapper_annotations_file <- snakemake@input[["emapper_annotations_file"]]
emapper_insecta_annotations_file <- snakemake@input[["emapper_insecta_annotations_file"]]

dds_file <- snakemake@output[["dds_file"]]
deseq2_results_file <- snakemake@output[["deseq2_results_file"]]
deseq2_counts_file <- snakemake@output[["deseq2_counts_file"]]
deseq2_results_with_annotations_file <- snakemake@output[["deseq2_results_with_annotations_file"]]
my_results_wtrans_ids_file <- snakemake@output[["my_results_wtrans_ids"]]
tx2gene_file <- snakemake@output[["tx2gene_file"]]
my_results_dt_file <- snakemake@output[["my_results_dt_file"]]
samples <- read.table(samples_file, header = TRUE)
files <- file.path(salmon_dir, "quant", samples$sample, "quant.sf")

# read gff
gr <- import.gff3(gff, feature.type = c("exons", "CDS", "mRNA", "gene"))

# extract a data.frame of tx to gene
mrnas <- gr[gr$type == "mRNA",]
mrna_dt <- as.data.table(mcols(mrnas))
tx2gene <- data.frame(mrna_dt[, .(TXNAME = ID, GENEID = as.character(Parent))])
fwrite(tx2gene, tx2gene_file)

names(files) <- paste0(samples$sample)

#import quant files
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

#generate & save DESeq2 object
dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~group)
saveRDS(dds, dds_file)

# run DESeq2
dds_deseq <- DESeq(dds)
res <- results(dds_deseq)
res <- subset(res, padj < 0.01)
res_ord <- as.data.table(res[order(res$padj),])
cou <- as.data.frame(counts(dds_deseq, normalized=TRUE))

# write results table
fwrite(res_ord, deseq2_results_file)
fwrite(cou, deseq2_counts_file)

#join emapper annotations with DEseq2 results
# read results and add column names
my_results <- cbind(rownames(res), data.frame(res, row.names=NULL))
setnames(my_results, old = c('rownames(res)'), new = c('gene'))
my_results_dt <- as.data.table(my_results)
my_results_wtrans_ids <- merge(my_results_dt, tx2gene, by.x = 'gene', by.y = 'GENEID', all = TRUE)
setnames(my_results_wtrans_ids, old = c('TXNAME'), new = c('transcript'))
fwrite(my_results_dt, my_results_dt_file)
fwrite(my_results_wtrans_ids, my_results_wtrans_ids_file)
 
# read and prep. emapper annotations
emapper_annot <- read.delim(emapper_annotations_file, header=FALSE, comment.char="#")
#sep_emapper_annot <- separate(emapper_annot, "V1", c('gene', 'V1'), sep = "-T")
prepd_annot <- subset(emapper_annot, , -c(V2, V3, V4, V5, V7, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21))
prepd_annot <- prepd_annot[, c(1, 2, 4, 3, 5)]
colnames(prepd_annot) <- c('transcript', 'emapper_taxa', 'emapper_gene', 'emapper_description', 'emapper_GO')

# merge results with emapper annotations
merged <- merge(my_results_wtrans_ids, prepd_annot, by.x = "transcript", by.y = "transcript", all.x = TRUE, all.y = FALSE)
merged_ord <-as.data.table(merged[order(merged$padj),])

# write results with annotations
fwrite(merged_ord, deseq2_results_with_annotations_file)

# log
sessionInfo()
