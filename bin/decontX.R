#!/usr/bin/env Rscript

library(celda)
library(dplyr)
library(ggplot2)
library(glue)
library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(c("--matrix"), action = 'store', type = 'character', help = '[Required] Matrix'),
  make_option(c("--features"), action = 'store', type = 'character', help = '[Required] Features'),
  make_option(c("--barcodes"), action = 'store', type = 'character', help = '[Required] Barcodes'),
  make_option(c("--soup_matrix"), action = 'store', type = 'character', help = '[Required] Soup matrix'),
  make_option(c("--soup_features"), action = 'store', type = 'character', help = '[Required] Soup features'),
  make_option(c("--soup_barcodes"), action = 'store', type = 'character', help = '[Required] Soup barcodes'),
  make_option(c("--clusters"), action = 'store', type = 'character', help = '[Required] Initial cluster assignments'),
  make_option(c("--delta"), action = 'store', type = 'character', default='10,10', help = '[Optional] Default: 10,10'),
  make_option(c("--fix_delta"), action = 'store_true', type = 'logical', default=F, help = '[Optional] Default: estimate delta'),
  make_option(c("--prefix"), action = 'store', type = 'character', default = 'seurat.', help = '[Optional] Prefix of output files')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, description = '')
opts <- parse_args(option_parser)


args <- commandArgs(T)

MATRIX <- opts$matrix
FEATURES <- opts$features
BARCODES <- opts$barcodes
SOUP_MATRIX <- opts$soup_matrix
SOUP_FEATURES <- opts$soup_features
SOUP_BARCODES <- opts$soup_barcodes
PREFIX <- opts$prefix
CLUSTERS <- opts$clusters
DELTA <- as.numeric(strsplit(opts$delta, ',')[[1]])
ESTIMATE_DELTA <- !opts$fix_delta

read_mm <- function(matrix_file, features_file, barcodes_file) {
    tmp <- as(Matrix::readMM(matrix_file), 'dgCMatrix')
    features <- read.table(features_file, as.is=T, sep='\t', head=F)
    features <- apply(features, 1, function(x){paste(x, collapse='\t')})
    barcodes <- read.table(barcodes_file, as.is=T, sep='\t', head=F)
    barcodes <- apply(barcodes, 1, function(x){paste(x, collapse='\t')})
    dimnames(tmp) <- list(features, barcodes)
    return(tmp)
}


write_mm <- function(mat, prefix) {
    features <- sort(colnames(mat))
    barcodes <- sort(rownames(mat))
    feature_to_index <- seq(1, length(features))
    names(feature_to_index) <- features
    barcode_to_index <- seq(1, length(barcodes))
    names(barcode_to_index) <- barcodes
    stopifnot(!('barcode' %in% colnames(mat)))
    mat$barcode <- rownames(mat)
    long <- tidyr::gather(mat, key='feature', value='count', -barcode)
    long <- long[long$count>0,c('feature', 'barcode', 'count')]
    long$feature_idx <- feature_to_index[long$feature]
    long$barcode_idx <- barcode_to_index[long$barcode]

    features_file <- glue('{prefix}features.tsv')
    barcodes_file <- glue('{prefix}barcodes.tsv')
    matrix_file <- glue('{prefix}matrix.mtx')

    con <- file(matrix_file)
    writeLines(c('%%MatrixMarket matrix coordinate integer general', '%', paste(c(length(features), length(barcodes), nrow(long)), collapse=' ')), con)
    close(con)
    write.table(long[,c('feature_idx', 'barcode_idx', 'count')], file=matrix_file, sep='\t', quote=F, row.names=F, col.names=F, append = T)
    write.table(features, file=features_file, sep='\t', quote=F, row.names=F, col.names=F)
    write.table(barcodes, file=barcodes_file, sep='\t', quote=F, row.names=F, col.names=F)

    return(c(matrix_file, features_file, barcodes_file))
}


counts <- read_mm(MATRIX, FEATURES, BARCODES)
soup <- read_mm(SOUP_MATRIX, SOUP_FEATURES, SOUP_BARCODES)
clusters <- read.table(CLUSTERS, sep='\t', as.is=T, head=F, col.names=c('barcode', 'cluster'))
z <- clusters$cluster
names(z) <- clusters$barcode

#counts <- SingleCellExperiment(assays = list(counts = counts))
#soup <- SingleCellExperiment(assays = list(counts = soup))

decontaminate <- decontX(counts, background=soup, z=z, delta=DELTA, estimateDelta=ESTIMATE_DELTA)
saveRDS(decontaminate, "decontaminate.rds")

#decontaminate <- decontaminate@metadata$decontX
contamination <- as.data.frame(decontaminate$estimates$all_cells$contamination)
colnames(contamination) <- c('contamination')
contamination$nucleus <- rownames(contamination)
write.table(contamination[,c('nucleus', 'contamination')], file=glue('{PREFIX}contamination.txt'), append = F, quote=F, row.names=F, col.names=F, sep='\t')

MEDIAN <- median(contamination$contamination)
p <- ggplot(contamination) +
geom_histogram(aes(x=contamination)) +
geom_vline(xintercept = MEDIAN, color='red') +
xlab('DecontX-estimated contamination') +
ylab('Number of nuclei')
png(glue('{PREFIX}contamination-histogram.png'), height=5, width=6, units='in', res=300)
print(p)
dev.off()

new_counts <- as.data.frame(t(as.matrix(decontaminate$decontXcounts)))
new_counts_rounded <- round(new_counts)

write_mm(new_counts, prefix=glue('{PREFIX}decontaminated-counts.'))
write_mm(new_counts_rounded, prefix=glue('{PREFIX}decontaminated-counts-rounded.'))