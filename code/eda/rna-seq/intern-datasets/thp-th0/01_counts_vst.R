# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LIBRARIES
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("dplyr", "yaml", "naturalsort")
.bioc_packages = c("edgeR", "rtracklayer", "SummarizedExperiment", "DESeq2")

## library(ggplotify)

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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MANIFEST
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")

pdata = read.csv2("assets/izi-run-phenodata-thp-th0.csv", header = T)
htseq.dir = paste0(manifest$izi_thp_th0_main, "htseq_count")
gc29.ftrs = readRDS(paste0(manifest$workdata, manifest$gencode_v29_features))

pdata = pdata %>% dplyr::mutate(
  PAIRS = case_when(
    TIMEPOINT == "0h" ~ paste0(TIMEPOINT,"_",DONOR),
    TRUE ~ paste0(TIMEPOINT,"_",DONOR)
  )
)

pdata = pdata %>% dplyr::mutate(
  SAMPLE_ID = case_when(
    TIMEPOINT == "0h" ~ paste0(TIMEPOINT,"_",DONOR),
    CONDITION == "Medium" & TIMEPOINT != "0h" ~ paste0(TIMEPOINT,"_M_",DONOR),
    TRUE ~ paste0(TIMEPOINT,"_A_",DONOR)
  )
)

pdata = pdata %>% dplyr::mutate(
  GROUP_BASE = case_when(
    TIMEPOINT == "0h" ~ "naive",
    CONDITION == "Medium" & TIMEPOINT != "0h" ~ "control",
    TRUE ~ "activated"
  )
)

pdata$GROUP = paste0(pdata$GROUP_BASE, "_", pdata$TIMEPOINT)

pdata = pdata %>%
  mutate(
    CONDITION =case_when(
      GROUP_BASE == "naive" ~ "thp",
      GROUP_BASE == "control" ~ "ctr",
      TRUE ~ "th0"
  )
)

pdata$TISSUE = "PBMC"
pdata$STUDY = "IZI"
pdata = mutate_if(pdata, is.character, as.factor)
pdata$TIMEPOINT = factor(pdata$TIMEPOINT, levels = naturalsort(levels(pdata$TIMEPOINT)))
pdata$GROUP_BASE = factor(pdata$GROUP_BASE, levels = c("naive", "control", "activated"))
pdata$GROUP = factor(pdata$GROUP, levels = naturalsort(levels(pdata$GROUP)))
pdata$GROUP = relevel(pdata$GROUP, ref = "naive_0h")
pdata$SHORT_NAME = paste0(pdata$CONDITION, "_", pdata$TIMEPOINT)
pdata$SHORT_NAME_REP = paste0(pdata$SHORT_NAME, "_", pdata$DONOR)
pdata = mutate_if(pdata, is.character, as.factor)
pdata = droplevels(pdata)

pdata = pdata[naturalorder(pdata$TIMEPOINT), ]
rownames(pdata) = pdata$BARCODE

colnames(pdata)[colnames(pdata) == "TIMEPOINT"] = "HOURS"
# colnames(pdata)[colnames(pdata) == "BARCODE"] = "SAMPLE_NAME"
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# COUNT DF
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
htseq = data.frame(
  FILES = list.files(
    htseq.dir, pattern = "RIB.*-htseq_counts.txt$", full=T, recursive = T
  ),
  stringsAsFactors = F
)

x = readDGE(
  htseq$FILES, columns=c(1,2), header=FALSE, labels = gsub("-htseq_counts.txt", "", basename(htseq$FILE))
)
raw.counts = as.data.frame( head(x$counts, -5) )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SummarizedExperiment obj
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
gc29.ftrs = gc29.ftrs[rownames(raw.counts), ]
raw.counts = raw.counts[, rownames(pdata)]
identical(rownames(gc29.ftrs), rownames(raw.counts))
identical(rownames(pdata), colnames(raw.counts))

se = SummarizedExperiment(
  assays = list(raw = as.matrix(raw.counts)), rowData = gc29.ftrs, colData = pdata)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Save csv and SummarizedExperiment obj
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

saveRDS(se, file = paste0(manifest$workdata, manifest$htseq_counts_izi_thp_th0))
