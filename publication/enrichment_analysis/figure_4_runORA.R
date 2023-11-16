.cran_packages = c("miscTools", "yaml", "NMF", "ggrepel", "ggplot2", "dplyr",
                   "viridis", "cowplot", "naturalsort", "scico", "scales", "inlmisc")
.bioc_packages = c("SummarizedExperiment", "pheatmap", "ComplexHeatmap")

# backports, vctrs

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

library(clusterProfiler)
library(ReactomePA)

source("code/R/eda-plots.R")
source("paper-helper.R")
theme_set(mytheme())

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

gc.anno = readRDS(paste0(work, manifest$gencode_v29_features))
signtrs = readRDS(paste0(work, manifest$signatures_all))

d.signtrs = subset(signtrs, Housekeeping == TRUE)
# d.signtrs = subset(signtrs)
v1.signtrs = subset(signtrs, v1_pass == TRUE & Housekeeping == TRUE)
v1v2.signtrs = subset(signtrs, v1_pass == TRUE & v2_pass == TRUE & v2_ctr_pass == TRUE & Housekeeping == TRUE)

entrez_id = function(ftrs) {
  entrez.sets = list()
  for (i in unique(ftrs$CLUSTER)) {
    ftrs.top = subset(ftrs, CLUSTER == i)
    ftrs.top.entrez = gc.anno$ENTREZ[match(ftrs.top$GENE_SYMBOL, gc.anno$GENE_SYMBOL_DUPL_MARKED)]
    ftrs.top.entrez = ftrs.top.entrez[!is.na(ftrs.top.entrez)]
    entrez.sets[[i]] = ftrs.top.entrez
  }
  entrez.sets = entrez.sets[lengths(entrez.sets) > 0]
  entrez.sets
}

d.entrez.sets = entrez_id(d.signtrs)
v1.entrez.sets = entrez_id(v1.signtrs)
v1v2.entrez.sets = entrez_id(v1v2.signtrs)

params = list(p_value = 0.05, minGSSize = 10, maxGSSize = 500)

ora_run = function(query, out) {
  res.cP.cc.go = compareCluster(
    geneCluster   = query,
    fun           ="enrichGO",
    minGSSize     = params$minGSSize,
    maxGSSize     = params$maxGSSize,
    OrgDb         ="org.Hs.eg.db",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = params$p_value,
    readable      = T
  )
  res.cP.cc.go.simp = simplify(res.cP.cc.go)
  saveRDS(res.cP.cc.go.simp, file = out)
}

ora_run(d.entrez.sets, "enrichment_analysis/figure_4_ora_go_d.rds")
ora_run(v1.entrez.sets, "enrichment_analysis/figure_4_ora_go_v1.rds")
ora_run(v1v2.entrez.sets, "enrichment_analysis/figure_4_ora_go_v1v2.rds")

ora_run = function(query, out) {
  res.cP.reactome = compareCluster(
    geneCluster   = query,
    fun           = "enrichPathway",
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = params$p_value
  )
  saveRDS(res.cP.reactome, file = out)
}

ora_run(d.entrez.sets, "enrichment_analysis/figure_4_ora_reactome_d.rds")
ora_run(v1.entrez.sets, "enrichment_analysis/figure_4_ora_reactome_v1.rds")
ora_run(v1v2.entrez.sets, "enrichment_analysis/figure_4_ora_reactome_v1v2.rds")
