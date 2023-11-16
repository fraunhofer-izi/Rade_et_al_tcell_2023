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

counts = read.table(paste0(work, "studies/data/Arcelus_GSE140244/GSE140244_rnaseq_gene_counts.txt"), row.names = 1, header = T)
pd =  read.table(paste0(work, "studies/data/Arcelus_GSE140244/GSE140244_rnaseq_meta_data.txt"), header = T)

gc.v19.anno = readRDS(paste0(work, manifest$gencode_v19_features))
# Remove PAR locus genes (45 genes). This avoids EnemblID duplications after Suffix truncation
gc.v19.anno = gc.v19.anno[!grepl("PAR_Y$", gc.v19.anno$ENSEMBL_ID), ]
# Remove Ensembl suffix
rownames(gc.v19.anno) = gsub("\\..+", "", rownames(gc.v19.anno))

lfc = 0
alpha = 0.05

out = paste0(work, manifest$limma_arcelus_cd4_th0)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rownames(pd) = pd$ExpressionMatrix_SampleID
pd = pd[colnames(counts), ]
pd$SAMPLE_NAME = pd$ExpressionMatrix_SampleID
pd$TISSUE = "PBMC"
pd$STUDY = "Arcelus"
pd$PLATFORM = "RNA-SEQ"
pd$HOURS = paste0(pd$Time_point, "h")

pd = data.frame(pd) %>% dplyr::mutate(
  SUBTYPE = case_when(
    HOURS == "0h" ~ "thp",
    TRUE ~ "th0"
  )
)
pd$SHORT_NAME = paste0(pd$SUBTYPE, "_", pd$HOURS)
pd$SHORT_NAME_REP = paste0(pd$SHORT_NAME, "_", pd$Replicate)
pd$ID = paste0("Th0Mem_", pd$HOURS)

pd = pd %>% dplyr::select("SAMPLE_NAME", "TISSUE", "SHORT_NAME_REP", "SHORT_NAME", "STUDY",
                          "HOURS", "PLATFORM", "RNAseq_plate", "Total_fragments", "Percent_mappedOK")


pd = pd[!grepl("_B$", pd$SAMPLE_NAME), ] # Removing replicate B !!!
counts = counts[, colnames(counts) %in% pd$SAMPLE_NAME]
rownames(pd) = pd$SAMPLE_NAME
pd = pd[colnames(counts), ]
pd$RNAseq_plate = paste0("Plate_", pd$RNAseq_plate)

pd = mutate_if(pd, is.character, as.factor)
lvls = naturalsort(as.character(unique(pd$HOURS)))
pd$HOURS = factor(pd$HOURS, levels = lvls)
pd$BIOLOGICAL_REPLICATE = factor(paste0("R_", gsub("\\_.+", "", pd$SAMPLE_NAME)))
pd = droplevels(pd)
rownames(pd) = pd$SAMPLE_NAME

identical(rownames(pd), colnames(counts))

dds = DGEList(counts = counts, samples = pd)
dds.deseq = DESeqDataSetFromMatrix(counts, pd, ~ 1)

# CPM
# dds.norm = calcNormFactors(dds, method = "TMM")
norm.cpm = cpm(dds, log = T, prior.count = 0.5)
# norm.cpm = limma::removeBatchEffect(
#   norm.cpm,
#   batch = dds$samples$RNAseq_plate,
#   design = model.matrix( ~ dds$samples$HOURS)
# )

gc.v19.anno = gc.v19.anno[rownames(dds.deseq), ]
identical(rownames(gc.v19.anno), rownames(dds.deseq))

se = SummarizedExperiment(
  assays=list(raw = as.matrix(counts(dds.deseq))),
  rowData = gc.v19.anno, colData = colData(dds.deseq)
)

assays(se)$norm.cpm = norm.cpm

# Removing genes that are lowly expressed
group = dds$samples$HOURS
keep.exprs = filterByExpr(dds, group=group)
dds = dds[keep.exprs,, keep.lib.sizes=FALSE]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TMM norm factors
dds = calcNormFactors(dds, method = "TMM")

# Make design matrix
design = model.matrix(~ 0 + BIOLOGICAL_REPLICATE + HOURS,  data = dds$samples)
colnames(design) = gsub("HOURS", "h_", colnames(design))
colnames(design) = gsub("BIOLOGICAL_REPLICATER", "R", colnames(design))

contrast = makeContrasts(
  HOURS_2h_vs_0h = h_2h,
  HOURS_4h_vs_0h = h_4h,
  HOURS_8h_vs_0h = h_8h,
  HOURS_12h_vs_0h = h_12h,
  HOURS_24h_vs_0h = h_24h,
  HOURS_48h_vs_0h = h_48h,
  HOURS_72h_vs_0h = h_72h,
  levels = colnames(design)
)

v = voom(dds, design, plot = T)

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
for (i in c("2h", "4h", "8h", "12h", "24h", "48h", "72h")) {
  obj = dds
  obj = obj[, grepl(paste0("^0h$|^", i, "$"), obj$samples$HOURS)]
  obj$samples = droplevels(obj$samples)
  # obj = calcNormFactors(obj, method = "TMM")

  design.sub = model.matrix(~ 0 + BIOLOGICAL_REPLICATE + HOURS,  data = obj$samples)
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
  dge.res.contr$ENSEMBL_ID_ABBR = rownames(dge.res.contr)
  dge.res.contr = dplyr::left_join(dge.res.contr, gc.v19.anno, by = "ENSEMBL_ID_ABBR")
  rownames(dge.res.contr) = dge.res.contr$ENSEMBL_ID_ABBR
  dge.res[[i]] = dge.res.contr
}

dge.res.sign.005 = lapply(dge.res, function(x) {
  x = subset(x,  abs(confect) > lfc)
  print((nrow(x)))
  return(x)
})

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SAVE
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj1 = list(
  se = se,
  dge.res = dge.res,
  dge.res.sign.005 = dge.res.sign.005,
  bayes.sd = bayes.sd
)
saveRDS(obj1, file = out)

