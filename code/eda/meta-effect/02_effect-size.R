# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LIBRARIES
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("ggplot2","devtools", "dplyr", "yaml", "naturalsort", "esc",
                   "meta", "metafor", "parallel", "reshape2", "tidyr", "forcats")
.bioc_packages = c("SummarizedExperiment", "limma", "ComplexHeatmap", "topconfects")

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

source("code/eda/studies.R")
source("code/eda/meta-effect//helper.R")
source("code/R/Rmarkdown-style.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  effect.size.all = F
} else {
  effect.size.all = T
}

grp_mean_sd = function(obj) {
  group.mean.l = list()
  grp.sd.l = list()
  for (i in paste0("^", unique(colnames(obj)), "$")) {

    grp = obj[, grepl(i , colnames(obj))]

    grp.sd = data.frame(SD = matrixStats::rowSds(grp), row.names =  rownames(grp))
    colnames(grp.sd) = unique(colnames(grp))

    grp.mean = data.frame(rowMedians(grp), row.names = rownames(grp))
    colnames(grp.mean) = unique(colnames(grp))

    group.mean.l[[unique(colnames(grp))]] = grp.mean
    grp.sd.l[[unique(colnames(grp))]] = grp.sd
  }
  group.mean.df = do.call("cbind", group.mean.l)
  group.sd.df = do.call("cbind", grp.sd.l)
  return(list(median = group.mean.df, sd = group.sd.df))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MANIFEST
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

if(effect.size.all == T) {
  se = readRDS(paste0(work, manifest$meta_se_w_izi_arcelus))
} else {
  se = readRDS(paste0(work, manifest$meta_se))
  # Meta-Analysis based on "Discovery Set"
  study.l$Th0Cd2 = NULL
  study.l$Th0Cd4Mem = NULL
}

alpha = 0.05
lfc = 0

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sources
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# how to compute cohen's d and se:
# http://supportupgrade.bioconductor.org/p/70175/
# http://supportupgrade.bioconductor.org/p/71747/

# Toro-Domínguez et al,
# "A survey of gene expression meta-analysis: methods and applications",
# Briefings in Bioinformatics, https://doi.org/10.1093/bib/bbaa019

# Brown et al.,
# "Meta-analysis of transcriptomic datasets identifies genes enriched in the mammalian circadian pacemaker.",
# Nucleic Acids Res. 2017, doi: 10.1093/nar/gkx714. PMID: 28973476; PMCID: PMC5737434.
# see Supplement

# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/random.html

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# studies[Cellsubtype][contrast]
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# st = names(study.l)[4]

studies = list()
studies.all.genes = list()
print("Contrast analyzed in the meta-analysis")
for (st in names(study.l)) {
  print("-------")
  print(st)
  obj = readRDS(paste0(work, manifest[[ study.l[[st]]$dge.res ]]))
  res = obj$dge.res
  bayes.sd = obj$bayes.sd
  bayes.sd = bayes.sd[!duplicated(gsub("\\..+", "", rownames(bayes.sd))), ]
  rownames(bayes.sd) = gsub("\\..+", "", rownames(bayes.sd))
  colnames(bayes.sd) = gsub("_0.5h", "_05h", colnames(bayes.sd))

  if(st == "Th1" | st == "Th2") {
    exprs.mat = exprs(obj$eset.core)
    colnames(exprs.mat) = as.character(pData(obj$eset.core)$SAMPLE_NAME)
    colnames(exprs.mat) = rownames(colData(se))[match(colnames(exprs.mat), se$SAMPLE_NAME)]
  } else {
    exprs.mat = assays(se)$cpm
  }

  # colnames(exprs.mat) = gsub(".+\\.", "", colnames(exprs.mat))
  # exprs.sd = grp_mean_sd(exprs.mat)$sd

  samples = study.l[[st]]$samples

  res.ctrst = list()
  res.ctrst.all = list()
  for (ctrst in names(res)) {
    res.sub = res[[ctrst]]

    # d_value based on bayes sd START
    bayes.sd.ctrst = bayes.sd[rownames(bayes.sd) %in% res.sub$ENSEMBL_ID_ABBR, ]
    bayes.sd.ctrst = bayes.sd.ctrst[res.sub$ENSEMBL_ID_ABBR, ]
    bayes.sd.ctrst = data.frame(bayes.sd.ctrst[, ctrst, drop = F])
    res.sub$sd_bayes = bayes.sd.ctrst[[1]]
    res.sub$d_value = res.sub$logFC /res.sub$sd_bayes
    # d_value based on bayes sd END

    # Correction factor
    n1 = samples[[ctrst]][1]
    n2 = samples[[ctrst]][2]
    j = (1 - ((3) / (4 * (n1 + n2) - 9)))
    res.sub$g_value = res.sub$d_value * j

    res.sub$var_d = ( (n1 + n2) / (n1 * n2) ) + ((res.sub$d_value^2) / (2*(n1 + n2)) )
    res.sub$var_g = ( (j^2) * (res.sub$var_d) )

    res.sub = res.sub %>%
      dplyr::select(
        ENSEMBL_ID, ENSEMBL_ID_ABBR, GENE_SYMBOL = GENE_SYMBOL_DUPL_MARKED, logFC,
        P.Value, adj.P.Val = fdr_zero, sd_bayes = sd_bayes, CI.L, CI.R, AveExpr, d_value, g_value,
        var_g, confect, rank)

    res.ctrst.all[[ctrst]] = res.sub
    # res.sub = res.sub[res.sub$adj.P.Val < alpha, ]
    res.sub = subset(res.sub,  abs(confect) > lfc)
    print(paste0(ctrst, " -> ", nrow(res.sub)))
    res.ctrst[[ctrst]] = res.sub

    stopifnot(rownames(res.sub) == res.sub$ENSEMBL_ID_ABBR)
  }
  studies[[st]] = res.ctrst
  studies.all.genes[[st]] = res.ctrst.all
}

if(effect.size.all == T) {
  saveRDS(list(studies = studies, studies.all.genes = studies.all.genes), file = paste0(work, manifest$meta_effect_sizes_all))
  stop("This is the end")
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# meta.in[contrast][Cellsubtype]
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ctrsts = c(
  "_05h_vs_0h",
  "_1h_vs_0h",
  "_2h_vs_0h",
  "_4h_vs_0h",
  "_6h_vs_0h",
  "_12h_vs_0h",
  "_24h_vs_0h",
  "_48h_vs_0h",
  "_72h_vs_0h")

meta.in = list()
for (ctrst in ctrsts) {
  hr = gsub("_", "", gsub("_vs.+", "", ctrst))
  st.c  = list()
  for (st in names(studies)) {
    if (any(grepl(ctrst, names(studies[[st]])))) {
      df = studies[[st]][[ grep(ctrst, names(studies[[st]])) ]]
      df$GROUP = hr
      st.c[[st]] = df
    }
  }
  meta.in[[hr]] = st.c
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# meta.in[contrast][Cellsubtype]
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
meta.in.fltrd = list()
print("-------")
print("Nbr. of genes analyzed")
for (ctrst in names(meta.in)) {
  ftrs = lapply(meta.in[[ctrst]], function(x){x$ENSEMBL_ID_ABBR})
  m = make_comb_mat(ftrs)
  m = m[comb_degree(m) > 1]
  ftrs.top = lapply(comb_name(m), function(nm) extract_comb(m, nm))
  ftrs.top = unlist(ftrs.top)

  print(paste0(ctrst, " -> ", length(ftrs.top)))
  st.l = list()
  for (st in names(meta.in[[ctrst]])) {
    df = meta.in[[ctrst]][[st]]
    rownames(df) = NULL
    df = df[df$ENSEMBL_ID_ABBR %in% ftrs.top, ]
    rownames(df) = df$ENSEMBL_ID_ABBR
    st.l[[st]] = df
  }
  meta.in.fltrd[[ctrst]] = st.l
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Run
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
meta.res = list()
meta.res.obj = list()
print("-------")
print("metafor:")
for (ctrst in names(meta.in.fltrd)) {
  print(ctrst)
  cl.res = list()
  cl.res.obj = list()
  query = meta.in.fltrd[[ctrst]]
  res = comb_effect(
    .obj = query, .geneID = "GENE_SYMBOL", .effectSize = "g_value", .se = "var_g",
    .ncores = 20, .confectsFDR = 0.05
  )

  print(table(!is.na(res$df$confect))[2])

  res$df$GROUP = ctrst
  meta.res[[ctrst]] = res$df
  meta.res.obj[[ctrst]] = res$obj
}

# Multiples testing
meta.res = lapply(meta.res, function(x) {
  x$est_pvalue_adj = p.adjust(x$est_pvalue, method = "BH")
  x
})

saveRDS(list(
  studies = studies,
  studies.all.genes = studies.all.genes,
  meta.in = meta.in,
  meta.in.fltrd = meta.in.fltrd,
  meta.res = meta.res),
  file = paste0(work, manifest$meta_all))

saveRDS(meta.res.obj, file = paste0(work, manifest$meta_objects))

