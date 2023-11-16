## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## LIBRARIES
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("yaml", "dplyr", "reshape2", "naturalsort")
.bioc_packages = c("SummarizedExperiment", "DESeq2", "limma", "edgeR", "topconfects")

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
# Manifest
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

se = readRDS(paste0(work, manifest$htseq_counts_izi_thp_th0))
# Remove PAR locus genes (45 genes). This avoids EnemblID duplications after Suffix truncation
se = se[!grepl("PAR_Y$", rownames(se)), ]
# Remove Ensembl suffix
rownames(se) = gsub("\\..+", "", rownames(se))

lfc = log2(1.5)
lfc = 0
alpha = 0.05

out.act = paste0(work, manifest$limma_izi_thp_th0_act)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load GC.v29 annotation (GENE_NAME, GENE_TYPE, GENE_TYPE_CLUSTER, etc)
# Clustered subbiotypes based on
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/_README_stats.txt
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
gc.v29.anno = readRDS(paste0(manifest$workdata, manifest$gencode_v29_features))
# Remove PAR locus genes (45 genes). This avoids EnemblID duplications after Suffix truncation
gc.v29.anno = gc.v29.anno[!grepl("PAR_Y$", gc.v29.anno$ENSEMBL_ID), ]
# Remove Ensembl suffix
rownames(gc.v29.anno) = gsub("\\..+", "", rownames(gc.v29.anno))

# Removing genes that are lowly expressed
dds = DGEList(counts = assay(se), samples = colData(se))

# For post logFC calculation
norm.cpm = cpm(dds, log = T, prior.count = 0.5)

# Removing genes that are lowly expressed
group = dds$samples$HOURS
keep.exprs = filterByExpr(dds, group=group)
dds = dds[keep.exprs,, keep.lib.sizes=FALSE]

assays(se)$cpm = norm.cpm

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TMM norm factors
dds = calcNormFactors(dds, method = "TMM")

# Make design matrix
design = model.matrix(~ 0 + DONOR + GROUP,  data = dds$samples)
colnames(design) = gsub("GROUP", "", colnames(design))

contrast <- makeContrasts(
  ACT_6h_vs_0h = activated_6h,
  ACT_12h_vs_0h = activated_12h,
  ACT_24h_vs_0h = activated_24h,
  ACT_48h_vs_0h = activated_48h,
  ACT_72h_vs_0h = activated_72h,
  levels = colnames(design)
)

v = voom(dds, design, plot = F)

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
for (i in c("6h", "12h", "24h", "48h", "72h")) {
  obj = dds
  obj = obj[, grepl(paste0("^0h$|^", i, "$"), obj$samples$HOURS)]
  obj = obj[, !grepl("control", obj$samples$GROUP_BASE)]
  obj$samples = droplevels(obj$samples)
  # obj = calcNormFactors(obj, method = "TMM")

  design.sub = model.matrix(~ 0 + DONOR + GROUP,  data = obj$samples)
  colnames(design.sub) = gsub("GROUPactivated_", "h_", colnames(design.sub))

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
# dge.res = list()
# for (i in 1:ncol(contrast)) {
#   ctrst = colnames(contrast)[[i]]
#   df = topTreat(efit, coef=ctrst, n = "Inf", adjust.method="BH", p.value=1, lfc=0, confint = T)
#   df$ENSEMBL_ID_ABBR = rownames(df)
#   dge.res[[ctrst]] = df
# }

# Annotation
for (i in 1:length(dge.res)) {
  dge.res.contr = dge.res[[i]]
  dge.res.contr = dplyr::left_join(dge.res.contr, gc.v29.anno, by = "ENSEMBL_ID_ABBR")
  rownames(dge.res.contr) = dge.res.contr$ENSEMBL_ID_ABBR
  dge.res[[i]] = dge.res.contr
}

dge.res.sign.005 = lapply(dge.res, function(x) {
  x = subset(x,  abs(confect) > lfc)
  print((nrow(x)))
  return(x)
})

# subsetting obj (columns)
geneID = "ENSEMBL_ID_ABBR"
obj.s = lapply(dge.res, function(x) {
  dplyr::select(x, all_of(geneID), all_of("confect"), all_of("logFC"))
})

# merging obj
for (st in names(obj.s)) {
  colnames(obj.s[[st]]) = paste(colnames(obj.s[[st]]), st, sep = "_")
  colnames(obj.s[[st]])[grep(geneID, colnames(obj.s[[st]]))] = geneID
}

obj.merged = Reduce(function(x, y)
  merge(x = x, y = y, by = geneID, all = TRUE), obj.s)
rownames(obj.merged) = obj.merged$ENSEMBL_ID

dge.res.confect = obj.merged[, grepl("confect", colnames(obj.merged))]
dge.res.logFC = obj.merged[, grepl("logFC", colnames(obj.merged))]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SAVE
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dge.res.act = dge.res[grepl("^ACT", names(dge.res))]
names(dge.res.act) = gsub("ACT", "HOURS", names(dge.res.act))
dge.res.sign.005.act = dge.res.sign.005[grepl("^ACT", names(dge.res.sign.005))]
names(dge.res.sign.005.act) = gsub("ACT", "HOURS", names(dge.res.sign.005.act))
obj = list(
  se = se,
  dge.res = dge.res.act,
  dge.res.sign.005 = dge.res.sign.005.act,
  dge.res.confect = dge.res.confect,
  dge.res.logFC = dge.res.logFC,
  bayes.sd = bayes.sd
)
saveRDS(obj, file = out.act)
