## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## LIBRARIES
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("dplyr")
.bioc_packages = c("GEOquery")
## library(ggplotify)

## Install CRAN packages (if not already installed)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Install bioconductor packages (if not already installed)
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}

list.of.packages = c(.cran_packages, .bioc_packages)

## Loading library
for(pack in list.of.packages) {
  suppressMessages(library(pack, quietly=TRUE, verbose=FALSE,
                           character.only=TRUE))
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## 
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sra.run.table = read.table("SraRunTable-original", sep = "\t", header = T)

gse = getGEO("GSE94396")
pd = gse$GSE94396_series_matrix.txt.gz
pd = pData(pd)

pd = subset(pd, `sample.group:ch1` == "G01" | `sample.group:ch1` == "G02" | `sample.group:ch1` == "G04")
sra.run.table = sra.run.table[sra.run.table$Sample_Name %in% rownames(pd), ]
sra.run.table$DATASTORE_provider = NULL
sra.run.table$DATASTORE_region = NULL

write.table(sra.run.table, file = "SraRunTable", sep = "\t", quote = F, row.names = F)

