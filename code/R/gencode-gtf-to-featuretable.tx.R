## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## LIBRARIES
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("yaml")
.bioc_packages = c("rtracklayer","dplyr")
## library(ggplotify)

## Install CRAN packages (if not already installed)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Install bioconductor packages (if not already installed)
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}

list.of.packages = c(.cran_packages, .bioc_packages)

## Loading library
for (pack in list.of.packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")

gc.vXX = rtracklayer::import(manifest$gencode_v32)
gc.vXX.table = dplyr::select(data.frame(subset(gc.vXX, type == "transcript")),
                             transcript_id, gene_id, seqnames, gene_name,
                             transcript_name, transcript_type)

## merge subbiotypes according to ftp://ftp.sanger.ac.uk/pub/gencode/_README_stats.txt
gc.vXX.table$TRANSCRIPT_TYPE_CLUSTER = gc.vXX.table$transcript_type
Long_non_coding_RNA_genes = c("processed_transcript", "lincRNA",
                              "3prime_overlapping_ncRNA", "antisense",
                              "non_coding",  "sense_intronic", "sense_overlapping",
                              "TEC", "known_ncrna", "macro_lncRNA",
                              "bidirectional_promoter_lncRNA")
Small_non_coding_RNA_genes = c("snRNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA",
                               "misc_RNA", "miRNA", "ribozyme", "sRNA",
                               "scaRNA", "scRNA", "vaultRNA")
gc.vXX.table$TRANSCRIPT_TYPE_CLUSTER = gsub(paste(Long_non_coding_RNA_genes, collapse = '|'),
                                      'lncRNA', gc.vXX.table$TRANSCRIPT_TYPE_CLUSTER)
gc.vXX.table$TRANSCRIPT_TYPE_CLUSTER = gsub(paste(Small_non_coding_RNA_genes, collapse = '|'),
                                      'small_ncRNA', gc.vXX.table$TRANSCRIPT_TYPE_CLUSTER)
gc.vXX.table$TRANSCRIPT_TYPE_CLUSTER  = gsub("^TR_.*","Ig_TcR", gc.vXX.table$TRANSCRIPT_TYPE_CLUSTER )
gc.vXX.table$TRANSCRIPT_TYPE_CLUSTER  = gsub("^IG_.*","Ig_TcR", gc.vXX.table$TRANSCRIPT_TYPE_CLUSTER )
gc.vXX.table$TRANSCRIPT_TYPE_CLUSTER  = gsub(".*pseudogene","pseudogene", gc.vXX.table$TRANSCRIPT_TYPE_CLUSTER )

colnames(gc.vXX.table) = c("ENSEMBL_TX_ID", "ENSEMBL_ID", "CHR", "GENE_SYMBOL",
                           "GENE_SYMBOL_TX", "BIOTYPE", "BIOTYPE_CLUSTER")
gc.vXX.table$ENSEMBL_ID_ABBR = gsub("\\..*","", gc.vXX.table$ENSEMBL_ID)
rownames(gc.vXX.table) = gc.vXX.table$ENSEMBL_TX_ID

saveRDS(gc.vXX.table, file = paste0(manifest$workdata, manifest$gencode_v32_features_tx))
