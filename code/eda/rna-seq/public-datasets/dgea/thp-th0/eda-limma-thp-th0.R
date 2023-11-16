.cran_packages = c("dplyr", "stringr", "devtools", "yaml", "naturalsort")
.bioc_packages = c("DESeq2", "limma", "edgeR", "topconfects")

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

# devtools::install_github("zhangyuqing/sva-devel")

list.of.packages = c(.cran_packages, .bioc_packages, c("sva"))

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
source("code/R/eda-plots.R")
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

eda.prep = readRDS(paste0(work, manifest$eda_prep_thp_thx))

lfc = log2(1.5)
alpha = 0.05
out = paste0(work, manifest$limma_thp_th0)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Load GC.v29 annotation (GENE_NAME, GENE_TYPE, GENE_TYPE_CLUSTER, etc)
## Clustered subbiotypes based on
## ftp://ftp.ebi.ac.uk/pub/databases/gencode/_README_stats.txt
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
gc.v29.anno = readRDS(paste0(work, manifest$gencode_v29_features))
# Remove PAR locus genes (45 genes). This avoids EnemblID duplications after Suffix truncation
gc.v29.anno = gc.v29.anno[!grepl("PAR_Y$", gc.v29.anno$ENSEMBL_ID), ]
# Remove Ensembl suffix
rownames(gc.v29.anno) = gsub("\\..+", "", rownames(gc.v29.anno))

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for (i in names(eda.prep$phenoData.l)) {
  eda.prep$phenoData.l[[i]] = subset(eda.prep$phenoData.l[[i]], CONDITION == "thp" | CONDITION == "th0")
  eda.prep$phenoData.l[[i]] = droplevels(eda.prep$phenoData.l[[i]])
  eda.prep$htseq.l[[i]] = eda.prep$htseq.l[[i]][, rownames(eda.prep$phenoData.l[[i]])]
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## EdgeR obj
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dds.l = list()
for (i in names(eda.prep$htseq.l)) {
  dds.l[[i]] = DGEList(counts = eda.prep$htseq.l[[i]], samples = eda.prep$phenoData.l[[i]])
}
dds = cbind(dds.l[[1]], dds.l[[2]], dds.l[[3]], dds.l[[4]])

# Remove PAR locus genes (45 genes). This avoids EnemblID duplications after Suffix truncation
dds$counts = dds$counts[!grepl("PAR_Y$", rownames(dds$counts)), ]
# Remove Ensembl suffix
rownames(dds$counts) = gsub("\\..+", "", rownames(dds$counts))

dds$samples$BIOLOGICAL_REPLICATE = paste0("R", dds$samples$BIOLOGICAL_REPLICATE, "_", dds$samples$STUDY)
lvls = naturalsort(unique(dds$samples$BIOLOGICAL_REPLICATE))
dds$samples$BIOLOGICAL_REPLICATE = factor(dds$samples$BIOLOGICAL_REPLICATE, levels = lvls)

# CPM
# dds.norm = calcNormFactors(dds, method = "TMM")
norm.cpm = cpm(dds, log = T, prior.count = 0.5)
norm.cpm = limma::removeBatchEffect(
  norm.cpm,
  batch = dds$samples$STUDY,
  design = model.matrix( ~ dds$samples$HOURS)
)

gc.v29.anno = gc.v29.anno[rownames(dds), ]
identical(rownames(gc.v29.anno), rownames(dds))
se = SummarizedExperiment(
  assays=list(raw = dds$counts, norm.cpm = norm.cpm),
  rowData = gc.v29.anno,  colData = dds$samples
)

# Removing genes that are lowly expressed
group = dds$samples$HOURS
keep.exprs = filterByExpr(dds, group=group)
dds.fltrd = dds[keep.exprs,, keep.lib.sizes=FALSE]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LIMMA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TMM norm factors
dds.fltrd = calcNormFactors(dds.fltrd, method = "TMM")

# Make design matrix
design = model.matrix(~ 0 + STUDY + HOURS,  data = dds.fltrd$samples)
colnames(design) = gsub("HOURS", "h_", colnames(design))

contrast = makeContrasts(
  HOURS_05h_vs_0h = h_0.5h,
  HOURS_1h_vs_0h = h_1h,
  HOURS_2h_vs_0h = h_2h,
  HOURS_4h_vs_0h = h_4h,
  HOURS_6h_vs_0h = h_6h,
  HOURS_12h_vs_0h = h_12h,
  HOURS_24h_vs_0h = h_24h,
  HOURS_48h_vs_0h = h_48h,
  HOURS_72h_vs_0h = h_72h,
  levels = colnames(design)
)

v = voom(dds.fltrd, design, plot = F)

# Fit model
vfit = lmFit(v, design)
# Compute cofficients for constrasts
vfit2 = contrasts.fit(vfit, contrast)
# Bayes shrinkage
efit = treat(vfit2, lfc = lfc)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Topconfects
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
topconf = list()
for (ctrst in colnames(vfit2$coefficients)) {
  topconf[[ctrst]] = limma_confects(vfit2, coef = ctrst, fdr=alpha, full = T)
}

dge.res = list()
for (i in 1:ncol(contrast)) {
  cond = colnames(contrast)[[i]]
  tt = topTreat(topconf[[1]]$limma_fit, coef=i, n = "Inf", adjust.method="BH", p.value=1, lfc=0, confint = T)
  tt =  data.frame(tt)

  tt$ENSEMBL_ID_ABBR = rownames(tt)
  confect.ctrst = topconf[[cond]]$table
  confect.ctrst = confect.ctrst %>% dplyr::select(rank, confect, effect, fdr_zero, name)
  tt = dplyr::left_join(tt, confect.ctrst, by = c("ENSEMBL_ID_ABBR" = "name"))
  rownames(tt) = tt$ENSEMBL_ID_ABBR
  head(tt)
  dge.res[[cond]] = tt
}

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## bayes sd for each contrast
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bayes.sd = list()
i = "72h"
for (i in c("0.5h", "1h", "2h", "4h", "6h", "12h", "24h", "48h", "72h")) {
  obj = dds.fltrd
  obj = obj[, grepl(paste0("^0h$|^", i, "$"), obj$samples$HOURS)]
  obj$samples = droplevels(obj$samples)
  # obj = calcNormFactors(obj, method = "TMM")

  design.sub = model.matrix(~ 0 + STUDY + HOURS,  data = obj$samples)
  colnames(design.sub) = gsub("HOURS", "h_", colnames(design.sub))

  v.sub = voom(obj, design.sub, plot = F)
  fit.sub = lmFit(v.sub, design.sub)
  cf = paste0("h_", i)
  ebayes.sub = eBayes(contrasts.fit(fit.sub, makeContrasts("CTRST" = cf, levels = colnames(design.sub))))
  ctrst = paste0("HOURS_", i , "_vs_0h")
  sd.post = sqrt(ebayes.sub$s2.post)
  names(sd.post) = rownames(ebayes.sub$coefficients)
  bayes.sd[[ctrst]] = sd.post
}
bayes.sd = do.call("cbind", bayes.sd)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## SUMMARY, FILTER SIGN. GENES, PRINT SIGN. REG. GENES
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Annotation
for (i in 1:length(dge.res)) {
  dge.res.contr = dge.res[[i]]
  dge.res.contr = dplyr::left_join(dge.res.contr, gc.v29.anno, by = "ENSEMBL_ID_ABBR")
  rownames(dge.res.contr) = dge.res.contr$ENSEMBL_ID_ABBR
  dge.res[[i]] = dge.res.contr
}

dge.res.sign.005 = lapply(dge.res, function(x) {
  # x = subset(x,  adj.P.Val < alpha)
  x = subset(x,  abs(confect) > lfc)
  print((nrow(x)))
  return(x)
})

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SAVE
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj = list(
  se = se,
  dge.res = dge.res,
  dge.res.sign.005 = dge.res.sign.005,
  bayes.sd = bayes.sd
)
saveRDS(obj, file = out)

