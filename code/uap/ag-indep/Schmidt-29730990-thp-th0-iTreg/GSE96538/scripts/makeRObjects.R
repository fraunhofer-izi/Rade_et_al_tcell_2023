#!/usr/bin/env Rscript

inFile = commandArgs(trailingOnly = TRUE)[1]
type = sub("\\..*", "", basename(inFile))
folder = dirname(inFile)
rawCounts = read.delim(inFile, header = T, row.names = 1)
Counts = rawCounts
if (sum(grepl("\\|", rownames(Counts))) == nrow(Counts)) {
    # rows contain tarnscriptID|geneID|etc.
    message("Gene are extracted from the rownames ...")
    rownames(Counts) = gsub(".*(ENST[^|]*).*", "\\1", rownames(Counts))
    genes = gsub(".*(ENSG[^|]*).*", "\\1", rownames(rawCounts))
} else {
    # we need to map to the gene names
    message("The assembly gtf is loaded to get gene names ...")
    suppressPackageStartupMessages(library(rtracklayer))
    indexDir = dirname(gsub("/counts/", "/index/", inFile))
    if (!any(grepl("\\.gtf$", dir(indexDir)))) {
        indexDir = dirname(indexDir)
    }
    gtfFile = grep("\\.gtf$", dir(indexDir), value = T)[1]
    gtfFile = paste0(indexDir, "/", gtfFile)
    assemb = import.gff(gtfFile)
    map = match(rownames(Counts), assemb$transcript_id)
    genes = assemb$gene_id[map]
}
save(Counts, file = paste0(folder, "/transcripts_", type, ".RData"))

Counts = aggregate(rawCounts, list(gene = genes), sum, drop = F)
rownames(Counts) = Counts[, "gene"]
Counts = Counts[, -1]
save(Counts, file = paste0(folder, "/genes_", type, ".RData"))
