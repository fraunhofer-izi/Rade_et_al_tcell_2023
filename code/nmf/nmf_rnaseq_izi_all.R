.cran_packages = c("dplyr", "yaml", "naturalsort")
.bioc_packages = c("DESeq2")
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

source("code/decomposition/nmf/nmf_helper.R")

randomize.x <- function(x){
  set.seed(12345)
  row.n = rownames(x)
  # resample the columns
  mat = apply(x, 2, function(c) sample(c, size=length(c)))
  rownames(mat) = row.n
  mat
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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

manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

meta.res = readRDS(paste0(work, manifest$meta_all))$meta.res

# lfc = log2(1.5)
lfc = 0

ftrs.flt = readRDS(paste0(work, manifest$signatures_discovery))$cd4.sgntrs
ftrs.flt$ENSEMBL_ID_ABBR = gsub("\\..+", "", ftrs.flt$ENSEMBL_ID)
# ftrs.flt = subset(ftrs.flt, CLUSTER != "0100")

nmf.th0 = readRDS(paste0(work, manifest$nmf_limma_th0_evo))

output.evo = paste0(work, manifest$nmf_limma_izi_thp_th0_all)
output.pre = paste0(work, manifest$nmf_limma_izi_thp_th0_all_pre)


limma.izi.cd2 = readRDS(paste0(work, manifest$limma_izi_thp_th0_act))
dge.sign = lapply(limma.izi.cd2$dge.res, function(x) {
  x = subset(x,  abs(confect) > lfc)
  print((nrow(x)))
  return(x)
})
dge.sign.union = unique(do.call(rbind.data.frame, dge.sign)$ENSEMBL_ID_ABBR)

se = limma.izi.cd2$se
se.nmf = se[, se$STUDY == "IZI"]
pd = data.frame(colData(se))
thp.act = pd[grepl("thp_|th0_", as.character(pd$SHORT_NAME)), ]
ctrl = pd[!grepl("thp_|th0_", as.character(pd$SHORT_NAME)), ]
se.nmf = se[, c(rownames(thp.act), rownames(ctrl))]
se$HOURS = gsub("h", "", se$HOURS)
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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Final Feature set from Meta-analysis (input for NMF)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Für diese Zeitpunkte habe ich 4 Subtypen
meta.res.1 = lapply(meta.res[c("6h")], function(x){
  res = x[x$ntimes == 4, ]
  res = subset(res, abs(confect) > lfc)$ENSEMBL_ID
})
meta.res.1 = unique(unlist(meta.res.1, use.names = F))

# Für diese Zeitpunkte habe ich 5 Subtypen
meta.res.2 = lapply(meta.res[c("12h","24h","48h","72h")], function(x){
  res = x[x$ntimes == 5, ]
  res = subset(res, abs(confect) > lfc)$ENSEMBL_ID
})
meta.res.2 = unique(unlist(meta.res.2, use.names = F))

meta.res.sign = union(meta.res.1, meta.res.2)
meta.res.sign = gsub("\\..+", "", meta.res.sign)
print(length(meta.res.sign))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
assays(se.nmf)$nmf  = quantile_norm(2^(assays(se.nmf)$cpm))
se.nmf.pre = se.nmf
rownames(se.nmf.pre) = rowData(se.nmf.pre)$GENE_SYMBOL_DUPL_MARKED
saveRDS(list(se = se.nmf.pre), file = output.pre)

# se.nmf = se.nmf[rownames(se.nmf) %in% ftrs.flt$ENSEMBL_ID_ABBR, ]
# rownames(se.nmf) = rowData(se.nmf)$GENE_SYMBOL_DUPL_MARKED

se.nmf = se.nmf[rownames(se.nmf) %in% rowData(nmf.th0$se)$ENSEMBL_ID_ABBR, ]
se.nmf = se.nmf[rownames(se.nmf) %in% dge.sign.union, ]
rownames(se.nmf) = rowData(se.nmf)$GENE_SYMBOL_DUPL_MARKED

# se.nmf = se.nmf[rownames(se.nmf) %in% dge.sign.union, ]
# se.nmf = se.nmf[rownames(se.nmf) %in% meta.res.sign, ]
# rownames(se.nmf) = rowData(se.nmf)$GENE_SYMBOL_DUPL_MARKED

nmf.exprs = assays(se.nmf)$nmf
print(nrow(nmf.exprs))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# NMF
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(NMF)

loss = list('brunet')

print("### Subtype: CD2 ###")
estim.r = nmf(nmf.exprs, ranks, nrun = runs, method = loss, .opt = threads, seed = 123456)
if (random.run == TRUE) {
  print("### Subtype: CD2 random ###")
  X.random = randomize.x(nmf.exprs)
  assays(se.nmf)$random = X.random
  estim.r.random = nmf(X.random, ranks, nrun = runs, method = loss, .opt = threads, seed = 123456)
  saveRDS(list(estim.r = estim.r, estim.r.random = estim.r.random, se = se.nmf), file = output.evo)
} else {
  saveRDS(list(estim.r = estim.r, se = se.nmf), file = output.evo)
}
