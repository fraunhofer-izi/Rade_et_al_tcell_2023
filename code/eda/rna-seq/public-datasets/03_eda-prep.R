.cran_packages = c("dplyr", "stringr", "devtools", "yaml")
.bioc_packages = c()

## Install CRAN packages (if not already installed)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Install bioconductor packages (if not already installed)
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  library(BiocManager)
  BiocManager::install(.bioc_packages[!.inst], ask = F)
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
home = manifest$homes
work = manifest$workdata

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## LOAD PHENO DATA
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
phenoData = readRDS(paste0(work, manifest$phenodata_rnaseq))
lvls = c("0h", "0.5h", "1h", "2h", "4h", "6h", "12h", "24h", "48h", "72h")

for (i in 1:length(phenoData)) {
  pD.gr = phenoData[[i]]
  pD.gr.n = names(phenoData)[[i]]
  rownames(pD.gr) = pD.gr$SAMPLE_NAME
  pD.gr = subset(
    pD.gr,
    pD.gr$HOURS == "0h" | pD.gr$HOURS == "0.5h" | pD.gr$HOURS == "1h" |
    pD.gr$HOURS == "2h" | pD.gr$HOURS == "4h" | pD.gr$HOURS == "6h" | pD.gr$HOURS == "12h" |
    pD.gr$HOURS == "24h" | pD.gr$HOURS == "48h" | pD.gr$HOURS == "72h"
  )
  pD.gr = droplevels(pD.gr)
  if (pD.gr.n != "schmidt") {
    pD.gr$HOURS = factor(pD.gr$HOURS, lvls)
  } else {
    pD.gr$HOURS = factor(pD.gr$HOURS, levels = c("0h", "2h", "6h", "24h", "48h"))
  }
  phenoData[[i]] = pD.gr
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## HTSEQ COUNTS
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
htseq.l = readRDS(paste0(work, manifest$htseq_counts_thp_thx))

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## SUBSET SAMLPES AND SAME ORDER FOR METADATA
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for (i in names(phenoData)) {
  htseq.l[[i]] = htseq.l[[i]][, rownames(phenoData[[i]])]
  # print(all.equal(colnames(htseq.l[[i]]), rownames(phenoData[[i]])))
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## SAVE
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj = list(htseq.l = htseq.l, phenoData.l = phenoData)
saveRDS(obj, file = paste0(work, manifest$eda_prep_thp_thx))

