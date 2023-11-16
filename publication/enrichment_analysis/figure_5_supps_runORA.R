# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries & Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("devtools","reshape2", "ggrepel", "stringr", "dplyr", "data.table", "yaml")
.bioc_packages = c("SummarizedExperiment", "edgeR")

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

source("code/R/eda-plots.R")
source("paper-helper.R")
theme_set(mytheme())

library(clusterProfiler)
library(ReactomePA)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

dge.act = readRDS(paste0(work, manifest$limma_izi_thp_th0_act))
dge.ctr = readRDS(paste0(work, manifest$limma_izi_ctrl))

dge.act.res = dge.act$dge.res
dge.ctr.res = dge.ctr$dge.res
rem = names(dge.ctr.res)[!grepl("^ACT_vs_CTR", names(dge.ctr.res))]
dge.ctr.res = dge.ctr.res[rem]

gc.anno = readRDS(paste0(work, manifest$gencode_v29_features))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dge.act.res.cp = dge.act.res
dge.ctr.res.cp = dge.ctr.res
names(dge.act.res.cp) = gsub("HOURS_|CTR_", "", names(dge.act.res.cp))
names(dge.ctr.res.cp) = gsub("HOURS_|CTR_", "", names(dge.ctr.res.cp))

l = list()
for (i in names(dge.act.res.cp)) {

  act = dge.act.res.cp[[i]]
  act = act[!is.na(act$confect), ]
  ctrl = dge.ctr.res.cp[[i]]
  ctrl = ctrl[!is.na(ctrl$confect), ]

  a = act[setdiff(act$ENSEMBL_ID_ABBR, ctrl$ENSEMBL_ID_ABBR), ]
  a.up = a[a$logFC > 0, ]
  a.up = a.up %>%
    dplyr::select(ENTREZ, ENSEMBL_ID_ABBR, logFC, BIOTYPE_CLUSTER) %>%
    dplyr::mutate(CTRST = gsub("_vs_.+", "", i), GROUP = "ACT_UP", GROUP_BASE = "ACT")
  a.dwn = a[a$logFC < 0, ]
  a.dwn = a.dwn %>%
    dplyr::select(ENTREZ, ENSEMBL_ID_ABBR, logFC, BIOTYPE_CLUSTER) %>%
    dplyr::mutate(CTRST = gsub("_vs_.+", "", i), GROUP = "ACT_DWN", GROUP_BASE = "ACT")

  c = ctrl[setdiff(ctrl$ENSEMBL_ID_ABBR, act$ENSEMBL_ID_ABBR), ]
  c.up = c[c$logFC > 0, ]
  c.up = c.up %>%
    dplyr::select(ENTREZ, ENSEMBL_ID_ABBR, logFC, BIOTYPE_CLUSTER) %>%
    dplyr::mutate(CTRST = gsub("_vs_.+", "", i), GROUP = "CTRL_UP", GROUP_BASE = "CTRL")
  c.dwn = c[c$logFC < 0, ]
  c.dwn = c.dwn %>%
    dplyr::select(ENTREZ, ENSEMBL_ID_ABBR, logFC, BIOTYPE_CLUSTER) %>%
    dplyr::mutate(CTRST = gsub("_vs_.+", "", i), GROUP = "CTRL_DWN", GROUP_BASE = "CTRL")

  act.up = act[act$logFC > 0, ]$ENSEMBL_ID_ABBR
  act.dwn = act[act$logFC < 0, ]$ENSEMBL_ID_ABBR
  ctrl.up = ctrl[ctrl$logFC > 0, ]$ENSEMBL_ID_ABBR
  ctrl.dwn = ctrl[ctrl$logFC < 0, ]$ENSEMBL_ID_ABBR

  i.up = intersect(act.up, ctrl.up)
  i.up.act = act[i.up, ] %>%
    dplyr::select(ENTREZ, ENSEMBL_ID_ABBR, logFC, BIOTYPE_CLUSTER) %>%
    dplyr::mutate(CTRST = gsub("_vs_.+", "", i), GROUP = "INT_UP", GROUP_BASE = "INT")
  # i.up.ctrl = ctrl[i.up, ] %>%
  #   dplyr::select(ENTREZ, ENSEMBL_ID_ABBR, logFC, BIOTYPE_CLUSTER) %>%
  #   dplyr::mutate(CTRST = gsub("_vs_.+", "", i), GROUP = "INT_UP_CTRL")

  i.diff = c(intersect(act.dwn, ctrl.up), intersect(act.up, ctrl.dwn))
  i.diff.act = act[i.diff, ] %>%
    dplyr::select(ENTREZ, ENSEMBL_ID_ABBR, logFC, BIOTYPE_CLUSTER) %>%
    dplyr::mutate(CTRST = gsub("_vs_.+", "", i), GROUP = "INT_DIFF", GROUP_BASE = "INT")
  # i.diff.ctrl = ctrl[i.diff, ] %>%
  #   dplyr::select(ENTREZ, ENSEMBL_ID_ABBR, logFC, BIOTYPE_CLUSTER) %>%
  #   dplyr::mutate(CTRST = gsub("_vs_.+", "", i), GROUP = "INT_DIFF_CTRL")


  i.dwn = intersect(act.dwn, ctrl.dwn)
  i.dwn.act = act[i.dwn, ] %>%
    dplyr::select(ENTREZ, ENSEMBL_ID_ABBR, logFC, BIOTYPE_CLUSTER) %>%
    dplyr::mutate(CTRST = gsub("_vs_.+", "", i), GROUP = "INT_DWN", GROUP_BASE = "INT")
  # i.dwn.ctrl = ctrl[i.dwn, ] %>%
  #   dplyr::select(ENTREZ, ENSEMBL_ID_ABBR, logFC, BIOTYPE_CLUSTER) %>%
  #   dplyr::mutate(CTRST = gsub("_vs_.+", "", i), GROUP = "INT_DWN_CTRL")


  c.all = ctrl[!ctrl$ENSEMBL_ID %in% a, ] %>%
    dplyr::select(ENTREZ, ENSEMBL_ID_ABBR, logFC, BIOTYPE_CLUSTER) %>%
    dplyr::mutate(CTRST = gsub("_vs_.+", "", i), GROUP = "CTRL_ALL", GROUP_BASE = "CTRL_ALL")
  c.all = c.all[!c.all$ENSEMBL_ID %in% i.diff, ]
  c.all.up = c.all[c.all$logFC > 0, ]
  c.all.up$GROUP = "CTRL_ALL_UP"
  c.all.dwn = c.all[c.all$logFC < 0, ]
  c.all.dwn$GROUP = "CTRL_ALL_DWN"

  l[[i]] = rbind(
    a.up,
    a.dwn,
    c.up,
    c.dwn,
    i.up.act,
    ## i.up.ctrl,
    i.diff.act,
    ## i.diff.ctrl,
    i.dwn.act,
    ## i.dwn.ctrl
    c.all = c.all,
    c.all.up = c.all.up,
    c.all.dwn = c.all.dwn

  )
}

df = do.call("rbind", l)
df = df[!is.na(df$ENTREZ), ]
df$ENTREZ = as.character(df$ENTREZ)
table(df$CTRST, df$GROUP_BASE)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
universe = gc.anno$ENTREZ
universe = as.character(universe[!is.na(universe)])

params = list(p_value = 0.05, minGSSize = 10, maxGSSize = 500)

# ora_run = function(query, out) {
#   res.cP.cc.go = compareCluster(
#     geneCluster   = ENTREZ ~ CTRST + GROUP,
#     data          = query,
#     universe      = universe,
#     fun           ="enrichGO",
#     minGSSize     = params$minGSSize,
#     maxGSSize     = params$maxGSSize,
#     OrgDb         ="org.Hs.eg.db",
#     ont           = "BP",
#     pAdjustMethod = "BH",
#     pvalueCutoff  = params$p_value,
#     readable      = T
#   )
#   res.cP.cc.go.simp = simplify(res.cP.cc.go)
#   saveRDS(res.cP.cc.go.simp, file = out)
# }
#
# ora_run(df, "enrichment_analysis/figure_5_supps_ora_go_act_neg_ctrl.rds")
#
# ora_run = function(query, out) {
#   res.cP.reactome = compareCluster(
#     geneCluster   = ENTREZ ~ CTRST + GROUP,
#     data          = query,
#     universe      = universe,
#     fun           = "enrichPathway",
#     organism      = "human",
#     pAdjustMethod = "BH",
#     pvalueCutoff  = params$p_value
#   )
#   saveRDS(res.cP.reactome, file = out)
# }
#
# ora_run(df, "enrichment_analysis/figure_5_supps_ora_reactome_act_neg_ctrl.rds")

################################################################################

# ora.in = subset(df, (GROUP_BASE == "ACT") | (df$GROUP == "CTRL_ALL_UP") | (df$GROUP == "CTRL_ALL_DWN"))
ora.in =subset(df, (GROUP_BASE == "ACT") | (df$GROUP == "CTRL_ALL_UP") | (df$GROUP == "CTRL_ALL_DWN") | (df$GROUP == "INT_DIFF"))
# ora.in = ora.in[abs(ora.in$logFC) > log2(2), ]
# table(ora.in$CTRST, ora.in$GROUP)


ora_run = function(query, out) {
  res.cP.cc.go = compareCluster(
    geneCluster   = ENTREZ ~ CTRST + GROUP,
    data          = query,
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
  saveRDS(res.cP.cc.go.simp, file = out)
}

ora_run(
  ora.in,
  "enrichment_analysis/figure_5_supps_ora_go_act_neg_ctrl_comb_2.rds"
)

ora_run = function(query, out) {
  res.cP.reactome = compareCluster(
    geneCluster   = ENTREZ ~ CTRST + GROUP,
    data          = query,
    universe      = universe,
    fun           = "enrichPathway",
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = params$p_value
  )
  saveRDS(res.cP.reactome, file = out)
}

ora_run(
  ora.in,
  "enrichment_analysis/figure_5_supps_ora_reactome_act_neg_ctrl_comb_2.rds"
)

# tmp = subset(df, (GROUP_BASE == "ACT") | (df$GROUP == "CTRL_ALL_UP") | (df$GROUP == "CTRL_ALL_DWN"))
# table(tmp$CTRST, tmp$GROUP)

