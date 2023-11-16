## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## LIBRARIES
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("dplyr", "yaml")
.bioc_packages = c("edgeR")

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

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
home = manifest$homes
work = manifest$workdata
uap.out = manifest$uap

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Tuomela-26967054-GSE52260
htseq.tuomela.dir = paste0(uap.out,manifest$tuomela_main)
htseq.tuomela = data.frame(FILES = list.files(htseq.tuomela.dir,
                                              pattern = "GSM.*-htseq_counts.txt$",
                                              full=T, recursive = T),
                           stringsAsFactors = F)

## Ullah-29466736-GSE9056
htseq.ullah.dir = paste0(uap.out,manifest$ullah_main)
htseq.ullah = data.frame(FILES = list.files(htseq.ullah.dir,
                                            pattern = "GSM.*-htseq_counts.txt$",
                                            full=T, recursive = T),
                         stringsAsFactors = F)

## Schmidt-29730990-GSE94396
htseq.schmidt.dir = paste0(uap.out,manifest$schmidt_main)
htseq.schmidt = data.frame(FILES = list.files(htseq.schmidt.dir,
                                              pattern = "GSM.*-htseq_counts.txt$",
                                              full=T, recursive = T),
                           stringsAsFactors = F)

## Schmidt-29730990-GSE96538
htseq.schmidt.ind.dir = paste0(uap.out,manifest$schmidt_ind_main)
htseq.schmidt.ind = data.frame(FILES = list.files(htseq.schmidt.ind.dir,
                                                  pattern = "GSM.*-htseq_counts.txt$",
                                                  full=T, recursive = T),
                               stringsAsFactors = F)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
studies.counts = c("htseq.tuomela", "htseq.ullah", "htseq.schmidt", "htseq.schmidt.ind")
htseq.counts.l = list()

for (i in 1:length(studies.counts)) {

  study = eval(parse(text = studies.counts[i]))

  x = readDGE(study$FILES, columns=c(1,2), header=FALSE,
              labels = gsub("-htseq_counts.txt", "", basename(study$FILE)))
  raw.counts = as.data.frame( head(x$counts, -5) )

  htseq.counts.l[[ gsub("htseq.","", studies.counts[i]) ]] = raw.counts
}

saveRDS(htseq.counts.l, file = paste0(work, manifest$htseq_counts))

