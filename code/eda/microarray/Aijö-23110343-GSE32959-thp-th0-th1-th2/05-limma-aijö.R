## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## LIBRARIES
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("yaml", "dplyr", "naturalsort")
.bioc_packages = c("limma", "affy", "topconfects")

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
## MANIFEST
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
workdata     = manifest$workdata

eda = readRDS(paste0(workdata, manifest$eda_prep_aijoe_thp_th1))
eset = eda$esets$rma

# Remove PAR locus genes (45 genes). This avoids EnemblID duplications after Suffix truncation
eset = eset[!grepl("PAR_Y$", rownames(eset)), ]
# Remove Ensembl suffix
rownames(eset) = gsub("\\..+", "", rownames(eset))

lfc = log2(1.5)
alpha = 0.05
out = paste0(workdata, manifest$limma_aijoe_thp_th1)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Remove AFFX and ERCC probesets
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
affxCtrls = grepl("^AFFX|^ERCC", rownames( exprs(eset)))
eset = eset[!affxCtrls, ]

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## ESPRESSION FILTER (PERCENTILE)
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
percentile = 0.25
cutoff = round(quantile(exprs(eset), percentile), 2)

eset.core = eset[rowSums(exprs(eset) > cutoff) >= 3, ]
pass.percentile.filter = nrow(exprs(eset.core))

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## ESPRESSION FILTER (IQR)
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
qn_iqr_auto = function(iqr_vals, binwidth = 0.05){
  istart = round(min(iqr_vals), digits = 1)
  iend = round(max(iqr_vals), digits = 1) + binwidth
  iqr_breaks = seq(from = istart, to=iend, by = binwidth)
  nbreaks = length(iqr_breaks)
  iqr_labels = iqr_breaks[2:nbreaks]
  iqr_tb = table(cut(iqr_vals, breaks=iqr_breaks, labels=iqr_labels))
  iqr_df = data.frame(iqr = iqr_labels, freq = as.vector(iqr_tb))
  iqr_df$csum = cumsum(iqr_df$freq)
  iqr_df$csum2 <- abs(iqr_df$csum - max(cumsum(iqr_df$freq)))
  iqr_df
}

M = exprs(eset.core)
Iqrs = apply(M, 1, IQR)

iqr.df = qn_iqr_auto(Iqrs)

iqr = iqr.df$iqr[iqr.df$freq == max(iqr.df$freq)]
if(length(iqr) > 1){
  iqr = iqr[length(iqr)]
}
Midx = as.vector(which(Iqrs > iqr, arr.ind=T))

eset.core = eset.core[Midx, ]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LIMMA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pData(eset.core)$BIOLOGICAL_REPLICATE = paste0("R",pData(eset.core)$BIOLOGICAL_REPLICATE)
lvls = naturalsort(unique(pData(eset.core)$BIOLOGICAL_REPLICATE))
pData(eset.core)$BIOLOGICAL_REPLICATE = factor(pData(eset.core)$BIOLOGICAL_REPLICATE, levels = lvls)

# Make design matrix
design = model.matrix(~ 0 + BIOLOGICAL_REPLICATE + HOURS,  data = pData(eset.core))
colnames(design) = gsub("HOURS", "h_", colnames(design))

contrast = makeContrasts(
  HOURS_12h_vs_0h = h_12h,
  HOURS_24h_vs_0h = h_24h,
  HOURS_48h_vs_0h = h_48h,
  HOURS_72h_vs_0h = h_72h,
  levels = colnames(design)
)

# Fit model
fit = lmFit(eset.core, design)
# Compute cofficients for constrasts
fit2 = contrasts.fit(fit, contrast)
# Bayes shrinkage
efit = treat(fit2, lfc = lfc)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Topconfects
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
topconf = list()
for (ctrst in colnames(fit2$coefficients)) {
  topconf[[ctrst]] = limma_confects(fit2, coef = ctrst, fdr=alpha, full = T)
}

dge.res = list()
i = 1
for (i in 1:ncol(contrast)) {
  cond = colnames(contrast)[[i]]
  tt = topTreat(topconf[[1]]$limma_fit, coef=i, n = "Inf", adjust.method="BH", p.value=1, lfc=0, confint = T)
  tt =  data.frame(tt)

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
for (i in c("12h", "24h", "48h", "72h")) {
  obj = eset.core
  obj = obj[, grepl(paste0("^0h$|^", i, "$"), obj$HOURS)]
  pData(obj) = droplevels(pData(obj))

  design.sub = model.matrix(~ 0 + BIOLOGICAL_REPLICATE + HOURS,  data=pData(obj))
  colnames(design.sub) = gsub("HOURS", "h_", colnames(design.sub))

  fit.sub = lmFit(obj, design.sub)
  cf = paste0("h_", i)
  ebayes.sub = eBayes(contrasts.fit(fit.sub, makeContrasts("CTRST" = cf, levels = colnames(design.sub))))
  ctrst = paste0("HOURS_", i , "_vs_0h")
  bayes.sd[[ctrst]] = sqrt(ebayes.sub$s2.post)
}
bayes.sd = do.call("cbind", bayes.sd)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## SUMMARY, FILTER SIGN. GENES, PRINT SIGN. REG. GENES
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# dge.res = list()
# for (i in 1:ncol(contrast)) {
#   ctrst = colnames(contrast)[[i]]
#   df = topTreat(efit, coef=ctrst, n = "Inf", adjust.method="BH", p.value=1, lfc=0, confint = T)
#   df$ENSEMBL_ID_ABBR = rownames(df)
#   dge.res[[ctrst]] = df
# }

dge.res.sign.005 = lapply(dge.res, function(x) {
  # x = subset(x,  adj.P.Val < alpha)
  x = subset(x,  abs(confect) > lfc)
  print((nrow(x)))
  return(x)
})

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## SAVE OBJECT
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj = list(
  eset = eset,
  eset.core = eset.core,
  dge.res = dge.res,
  dge.res.sign.005 = dge.res.sign.005,
  bayes.sd = bayes.sd
)
saveRDS(obj, file = out)
