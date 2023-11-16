# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LIBRARIES
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("ggplot2","devtools", "yaml", "naturalsort", "dplyr",
                   "reshape2", "tidyr", "ggrepel", "cowplot", "pathfindR")
.bioc_packages = c("SummarizedExperiment")

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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata
obj = readRDS(paste0(work, manifest$meta_all))
meta.res = obj$meta.res

se = readRDS(paste0(work, manifest$meta_se))

m.res = lapply(meta.res, function(ctrst) {
  ctrst %>% dplyr::select(GENE_SYMBOL, est_pvalue, est_pvalue_adj, est_ci_l, est_ci_r,
                          est_effect, est_se, het_I2, rank, confect, confect_fdr_zero,
                          GROUP, signcon, ntimes)
})
m.res.df = do.call("rbind", m.res)

gc.anno = readRDS(paste0(work, manifest$gencode_v29_features))


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Clusterprofiler
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tops = subset(m.res.df, abs(confect) > lfc)
tops$ID = paste0(tops$GENE_SYMBOL, "_", tops$GROUP)
tops.1 = tops[tops$GROUP == "05h" |tops$GROUP == "1h" |tops$GROUP == "2h" | tops$GROUP == "4h" |tops$GROUP == "6h", ]
tops.1 = subset(tops.1, ntimes >= 4)
tops.2 = tops[tops$GROUP == "12h" |tops$GROUP == "24h" |tops$GROUP == "48h" |tops$GROUP == "72h", ]
tops.2 = subset(tops.2, ntimes >= 5)
tops = rbind(tops.1, tops.2)
tops$GROUP = gsub("05h", "0.5h", tops$GROUP)

meta.res = split(tops, tops$GROUP)
meta.res = meta.res[naturalsort(names(meta.res))]

entrez.sets = lapply(meta.res, function(ctrst) {
  symbol = ctrst[!is.na(ctrst$confect), ]$GENE_SYMBOL
  entrez = gc.anno$ENTREZ[match(symbol, gc.anno$GENE_SYMBOL_DUPL_MARKED)]
  as.character(entrez[!is.na(entrez)])
})
names(entrez.sets) = paste0(1:length(names(entrez.sets)), "_", names(entrez.sets))
entrez.sets = entrez.sets[lengths(entrez.sets) > 0]

universe = gc.anno$ENTREZ
universe = as.character(universe[!is.na(universe)])

params = list(p_value = 0.05, minGSSize = 10, maxGSSize = 500)

res.cP.cc.reactome = compareCluster(
  geneCluster   = entrez.sets,
  universe      = universe,
  fun           = "enrichPathway",
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = params$p_value
)

saveRDS(res.cP.cc.reactome, file = "enrichment_analysis/meta_analyis_ora_reactome.rds")

res.cP.cc.go = compareCluster(
  geneCluster   = entrez.sets,
  universe      = universe,
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

saveRDS(res.cP.cc.go.simp, file = "enrichment_analysis/meta_analyis_ora_go.rds")

# dotplot(res.cP.cc.reactome, showCategory = 5)

