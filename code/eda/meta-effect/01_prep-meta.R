# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LIBRARIES
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("ggplot2","devtools", "dplyr", "yaml", "naturalsort", "reshape2",
                   "tidyr", "forcats")
.bioc_packages = c("DESeq2", "limma", "edgeR", "preprocessCore")

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

source("code/R/Rmarkdown-style.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MANIFEST
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

gc.ftrs = readRDS(paste0(manifest$workdata, manifest$gencode_v29_features))

source("code/eda/studies.R")
source("code/R/Rmarkdown-style.R")
source("code/R/eda-plots.R")
source("code/R/FSQN.R")
source("code/eda/meta-effect/helper.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  add.izi = F
  add.arcelus = F
} else {
  add.izi = args[1]
  add.arcelus = args[2]
  print(paste0("Adding izi run: ", args[1]))
  print(paste0("Adding Arcelus run: ", args[2]))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
studies = list()
pd.col = c("SAMPLE_NAME", "TISSUE", "SHORT_NAME_REP", "SHORT_NAME", "STUDY", "HOURS", "PLATFORM")

# names(study.l)
# st = names(study.l)[6]
for (st in names(study.l)) {
  print(st)
  obj = readRDS(paste0(work, manifest[[ study.l[[st]]$dge.res ]]))
  se = NULL

  if (st == "Th0" | st == "iTreg" | st == "Th17") {
    se = obj$se
    se$PLATFORM = "RNA-SEQ"
  }

  if (st == "Th0Cd2") {
    se = obj$se[, obj$se$GROUP_BASE != "control"]
    se$PLATFORM = "RNA-SEQ"
    colnames(colData(se))[colnames(colData(se)) == "BARCODE"] = "SAMPLE_NAME"
  }

  if (st == "Th0Cd4Mem") {
    se = obj$se
    se$SAMPLE_NAME = gsub("_", "-", se$SAMPLE_NAME)
  }

  if (st == "Th1" | st == "Th2") {
    exprs.obj = obj$eset
    se = SummarizedExperiment(
      assays = list(rma = as.matrix(exprs(exprs.obj))),
      colData = pData(exprs.obj)
    )
    se$PLATFORM = "ARRAY"
  }

  colData(se) = colData(se)[pd.col]
  se$SUBTYPE = as.factor(gsub("\\_.+", "", se$SHORT_NAME))
  colnames(se) = paste0(se$SAMPLE_NAME, ".", st, "_", gsub("0\\.5h", "05h", se$HOURS))
  colData(se) = droplevels(colData(se))
  rownames(se) = gsub("\\..+", "", rownames(se))
  print(nrow(se))
  studies[[st]] = se
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SE: merged Seq & Arrays
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Same gene set accros all subtypes
ftrs.all.1 = rownames(studies$Th1) # without arcelus
ftrs.all.2 = ftrs.all.1[ftrs.all.1 %in% rownames(studies$Th0Cd4Mem)] # with arcelus

if(add.arcelus == FALSE) {
  ftrs.all = ftrs.all.1
} else {
  ftrs.all = ftrs.all.2
}

exprs.all.cpm = cbind(
  assays(studies$Th0[ftrs.all, ])$norm.cpm,
  assays(studies$Th17[ftrs.all, ])$norm.cpm,
  assays(studies$iTreg[ftrs.all, ])$norm.cpm,
  if(add.izi == TRUE) {assays(studies$Th0Cd2[ftrs.all, ])$cpm},
  if(add.arcelus == TRUE) {assays(studies$Th0Cd4Mem[ftrs.all, ])$norm.cpm},
  assay(studies$Th2[ftrs.all, ]),
  assay(studies$Th1[ftrs.all, ])
)

exprs.all.raw = cbind(
  assays(studies$Th0[ftrs.all, ])$raw,
  assays(studies$Th17[ftrs.all, ])$raw,
  assays(studies$iTreg[ftrs.all, ])$raw,
  if(add.izi == TRUE) {assays(studies$Th0Cd2[ftrs.all, ])$raw},
  if(add.arcelus == TRUE) {assays(studies$Th0Cd4Mem[ftrs.all, ])$raw},
  assay(studies$Th2[ftrs.all, ]),
  assay(studies$Th1[ftrs.all, ])
)

pd.all = rbind(
  colData(studies$Th0),
  colData(studies$Th17),
  colData(studies$iTreg),
  if(add.izi == TRUE) {colData(studies$Th0Cd2)},
  if(add.arcelus == TRUE) {colData(studies$Th0Cd4Mem)},
  colData(studies$Th2),
  colData(studies$Th1)
)

pd.all$ID = gsub(".+\\.", "", rownames(pd.all))

gc.ftrs.cp = gc.ftrs
gc.ftrs.cp = gc.ftrs.cp[gc.ftrs.cp$ENSEMBL_ID_ABBR %in% ftrs.all, ]
rownames(gc.ftrs.cp) = gc.ftrs.cp$ENSEMBL_ID_ABBR
gc.ftrs.cp = gc.ftrs.cp[ftrs.all, ]

se = SummarizedExperiment(
  assays = list(
    # vst2 = exprs.all.vst2,
    # vst2adj = exprs.all.vst2adj,
    cpm = exprs.all.cpm,
    raw = exprs.all.raw
  ),
  colData = pd.all,
  rowData = gc.ftrs.cp
)

print(nrow(se))

if(add.arcelus == TRUE) {
  lvls = c("0h","0.5h","1h","2h","4h","6h","8h","12h","24h","48h","72h")
  se$HOURS = factor(se$HOURS, levels = lvls)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Batch correction
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
array = se[, se$PLATFORM == "ARRAY"]
colData(array) = droplevels(colData(array))
assays(array)$adj = limma::removeBatchEffect(
  assay(array), batch = array$STUDY, design = model.matrix( ~ array$HOURS))

seq = se[, se$PLATFORM == "RNA-SEQ"]
colData(seq) = droplevels(colData(seq))
seq.dgelist = DGEList(counts = assays(seq)$raw, samples = data.frame(colData(seq)))
seq.dgelist = calcNormFactors(seq.dgelist, method = "TMM")
seq.cpm = cpm(seq.dgelist, log = T, prior.count = 0.5, normalized.lib.sizes = T)

assays(seq)$adj = limma::removeBatchEffect(
  seq.cpm, batch = seq$STUDY, design = model.matrix( ~ seq$HOURS))

assays(se)$adj = cbind(assays(seq)$adj, assays(array)$adj)

fsqn = quantileNormalizeByFeature(t(assays(array)$adj), t(assays(seq)$adj))
assays(se)$fsqn = cbind(assays(seq)$adj, t(fsqn))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(add.izi == TRUE & add.arcelus == FALSE) {
  saveRDS(se, file = paste0(work, manifest$meta_se_w_izi))
} else if (add.izi == TRUE & add.arcelus == TRUE) {
  saveRDS(se, file = paste0(work, manifest$meta_se_w_izi_arcelus))
} else if (add.izi == FALSE & add.arcelus == FALSE) {
  saveRDS(se, file = paste0(work, manifest$meta_se))
}

