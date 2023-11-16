.cran_packages = c("dplyr", "yaml", "naturalsort")
.bioc_packages = c("limma", "SummarizedExperiment")

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
for(pack in list.of.packages) {
  suppressMessages(library(pack, quietly=TRUE, verbose=FALSE, character.only=TRUE))
}

source("code/decomposition/nmf/nmf_helper.R")

randomize.x <- function(x){
  set.seed(12345)
  row.n = rownames(x)
  # resample the columns
  mat = apply(x, 2, function(c) sample(c, size=length(c)))
  rownames(mat) = row.n
  mat
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Arguments
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  threads = 'vp15'
  runs = 200
  ranks = 2:7
  random.run = FALSE
} else {
  threads = args[1]
  runs = as.integer(args[2])
  ranks = args[3]
  r1 = as.numeric(strsplit(ranks,":")[[1]])[1]
  r2 = as.numeric(strsplit(ranks,":")[[1]])[2]
  ranks = r1:r2
  random.run = as.logical(as.integer(args[4]))

  print(paste0("# threads: ", args[1]))
  print(paste0("# runs: ", runs))
  print(paste0("# ranks: ", args[3]))
  print(paste0("# random.run: ", random.run))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MANIFEST / OBJECTS TO WORK
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

output.evo = paste0(work, manifest$nmf_aijoe_thp_th1_evo)

# lfc = log2(1.5)
lfc = 0

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj = readRDS(paste0(work, manifest$meta_all))
meta.res = obj$meta.res

se = readRDS(paste0(work, manifest$meta_se))

se$HOURS = gsub("h", "", se$HOURS)
se.nmf = se[, grepl("^Th1_", se$ID)]
colData(se.nmf) = droplevels(colData(se.nmf))
se.nmf = se.nmf[, naturalorder(se.nmf$HOURS)]
se.nmf$HOURS = factor(se.nmf$HOURS, levels = unique(se.nmf$HOURS))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load GC.v28 annotation (GENE_NAME, GENE_TYPE, GENE_TYPE_CLUSTER)
# Clustered subbiotypes based on
# ftp://ftp.sanger.ac.uk/pub/gencode/_README_stats.txt
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
gc.v29.anno = readRDS(paste0(work, manifest$gencode_v29_features))

meta.res = lapply(meta.res, function(x){
  x$ENSEMBL_ID = gc.v29.anno$ENSEMBL_ID[match(x$GENE_SYMBOL, gc.v29.anno$GENE_SYMBOL_DUPL_MARKED)]
  x
})

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Final Feature set from Meta-analysis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Für diese Zeitpunkte habe ich 5 Subtypen
meta.res.1 = lapply(meta.res, function(x){
  res = x[x$ntimes == 5, ]
  res = subset(res, abs(confect) > lfc)$ENSEMBL_ID
  # res = res[!is.na(res$confect), ]$ENSEMBL_ID
})
meta.res.sign = unique(unlist(meta.res.1, use.names = F))
meta.res.sign = gsub("\\..+", "", meta.res.sign)
print(length(meta.res.sign))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PARSE
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.nmf.th1 = se.nmf[meta.res.sign, ]
rownames(se.nmf.th1) = rowData(se.nmf.th1)$GENE_SYMBOL_DUPL_MARKED
# assays(se.nmf.th1)$nmf = t(apply(2^(assays(se.nmf.th1)$raw), 1, function(x) x/sum(x)))
assays(se.nmf.th1)$nmf = 2^(assays(se.nmf.th1)$raw)

nmf.exprs.th1 = assays(se.nmf.th1)$nmf

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## NMF
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(NMF)

loss = list('brunet')

estim.r = nmf(nmf.exprs.th1, ranks, nrun = runs, method = loss, .opt = threads, seed = 123456)
if (random.run == TRUE) {
  X.random = randomize.x(nmf.exprs.th1)
  estim.r.random = nmf(X.random, ranks, nrun = runs, method = loss, .opt = threads, seed = 123456)
  assays(se.nmf.th1)$random = X.random
  saveRDS(list(estim.r = estim.r, estim.r.random = estim.r.random, se = se.nmf.th1), file = output.evo)
} else {
  saveRDS(list(estim.r = estim.r, se = se.nmf.th1), file = output.evo)
}

