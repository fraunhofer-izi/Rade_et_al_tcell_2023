# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries & Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("plyr", "devtools","reshape2", "ggrepel", "stringr", "ggpmisc",
                   "dplyr", "scales", "yaml", "data.table",
                   "egg", "tibble", "NMF", "viridis", "cowplot")
.bioc_packages = c("SummarizedExperiment", "ComplexHeatmap")

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
source("code/R/eda-plots.R")
source("code/R/Rmarkdown-style.R")
theme_set(mytheme())
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load Data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

load("data/Housekeeping_GenesHuman.RData")
hk = as.character(Housekeeping_Genes$Gene.name)

# Load GC.v29 annotation (GENE_NAME, GENE_TYPE, GENE_TYPE_CLUSTER)
# Clustered subbiotypes based on
# ftp://ftp.sanger.ac.uk/pub/gencode/_README_stats.txt
gc.anno = readRDS(paste0(work, manifest$gencode_v29_features))

meta.prep = readRDS(paste0(work, manifest$meta_se))
# meta.prep.w.izi = readRDS(paste0(work, manifest$meta_se_w_izi))
meta.prep.w.izi.arcelus = readRDS(paste0(work, manifest$meta_se_w_izi_arcelus))
rownames(meta.prep.w.izi.arcelus) = rowData(meta.prep.w.izi.arcelus)$GENE_SYMBOL_DUPL_MARKED

signatures.discovery = readRDS(paste0(work, manifest$signatures_discovery))
cd4.sgntrs = signatures.discovery$cd4.sgntrs

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Arcelus run (Verification Set 1)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nmf.v1 = readRDS(paste0(work, manifest$nmf_limma_arcelus_evo))
nmf.v1$estim.r = nmf.v1$estim.r$brunet
nmf.v1$estim.r.random = nmf.v1$estim.r.random$brunet
rank.v1 = "5"

nmf.v1 = anno_metagenes(obj = nmf.v1, rank = rank.v1, group = "HOURS", reorder = TRUE)

# BASIS FEATURES
ftrs.v1.wrk = anno_w(nmf.v1, rank.v1)
# ftrs.v1.wrk = anno_w(nmf.v1, rank.v1, .rank_filter = c(1,2,3,5))
ftrs.v1.wrk.all = anno_w(nmf.v1, rank.v1)


cd4.sgntrs$CLUSTER = unlist(add_pseudo_mg_4to5(cd4.sgntrs$CLUSTER))
ftrs.v1.wrk$cluster = unlist(ch_mg_5(ftrs.v1.wrk$cluster))
# table(cd4.sgntrs$CLUSTER)
# table(ftrs.v1.wrk.all$cluster)

# Matching with clusters from the Discovery Set
ftrs.v1.wrk$cluster_discovery = cd4.sgntrs$CLUSTER[match(rownames(ftrs.v1.wrk), cd4.sgntrs$GENE_SYMBOL)]
ftrs.v1.wrk = ftrs.v1.wrk[!is.na(ftrs.v1.wrk$cluster_discovery), ] # Filtering out genes not in the discovery set
ftrs.v1.wrk$v1_pass = ftrs.v1.wrk$cluster == ftrs.v1.wrk$cluster_discovery

cd4.sgntrs$v1_pass = ftrs.v1.wrk$v1_pass[match(cd4.sgntrs$GENE_SYMBOL, rownames(ftrs.v1.wrk))]
cd4.sgntrs$v1_cluster = ftrs.v1.wrk$cluster[match(cd4.sgntrs$GENE_SYMBOL, rownames(ftrs.v1.wrk))]
cd4.sgntrs$v1_cluster[is.na(cd4.sgntrs$v1_pass)] = FALSE # entweder nicht signifikant oder nicht im Input
cd4.sgntrs$v1_pass[is.na(cd4.sgntrs$v1_pass)] = FALSE # entweder nicht signifikant oder nicht im Input

v1.sgntrs = cd4.sgntrs[cd4.sgntrs$v1_pass == T, ]


# tmp = v1.sgntrs
cd4.sgntrs[cd4.sgntrs$GENE_SYMBOL %in% c("TNFRSF9"), ]

cd4.sgntrs$CLUSTER = paste0(substring(cd4.sgntrs$CLUSTER, 1, 3), substring(cd4.sgntrs$CLUSTER, 5, 5))
cd4.sgntrs$v1_cluster = paste0(substring(cd4.sgntrs$v1_cluster, 1, 3), substring(cd4.sgntrs$v1_cluster, 5, 5))


saveRDS(list(nmf.v1 = nmf.v1, rank.v1 = rank.v1, ftrs.v1.wrk = ftrs.v1.wrk, ftrs.v1.wrk.all = ftrs.v1.wrk.all),
  file = paste0(work, manifest$pipeline_verification_v1)
)

stopifnot(identical(sort(subset(cd4.sgntrs, v1_pass == TRUE)$GENE_SYMBOL), sort(v1.sgntrs$GENE_SYMBOL)))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# izi run (Verification Set 2)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nmf.v2 = readRDS(paste0(work, manifest$nmf_limma_izi_thp_th0_evo))
nmf.v2$estim.r = nmf.v2$estim.r$brunet
nmf.v2$estim.r.random = nmf.v2$estim.r.random$brunet
rank.v2 = "3"

# plot(nmf.v2$estim.r)
# nmf_stats(obj = nmf.v2, anno = c("HOURS"))

nmf.v2 = anno_metagenes(obj = nmf.v2, rank = rank.v2, group = "HOURS", reorder = TRUE)
# kinetic_mg_profile(obj = nmf.v2, rank = rank.v2, tol10qualitative)

# BASIS FEATURES
# ftrs.v2.wrk = anno_w(nmf.v2, rank.v2, .rank_filter = c(1,2,4))
ftrs.v2.wrk = anno_w(nmf.v2, rank.v2)
ftrs.v2.wrk.all = ftrs.v2.wrk

# Matching with clusters from the Discovery Set
ftrs.v2.wrk$cluster_discovery = cd4.sgntrs$CLUSTER[match(rownames(ftrs.v2.wrk), cd4.sgntrs$GENE_SYMBOL)]
ftrs.v2.wrk = ftrs.v2.wrk[!is.na(ftrs.v2.wrk$cluster_discovery), ] # Filtering out genes not in the discovery set

ftrs.v2.wrk$cluster_discovery = paste0(
  substring(ftrs.v2.wrk$cluster_discovery, 1, 1),
  substring(ftrs.v2.wrk$cluster_discovery, 3, 4)
)
ftrs.v2.wrk$v2_pass = ftrs.v2.wrk$cluster == ftrs.v2.wrk$cluster_discovery
ftrs.v2.wrk = subset(ftrs.v2.wrk, cluster != "0100")

cd4.sgntrs$v2_pass = ftrs.v2.wrk$v2_pass[match(cd4.sgntrs$GENE_SYMBOL, rownames(ftrs.v2.wrk))]
cd4.sgntrs$v2_cluster = ftrs.v2.wrk$cluster[match(cd4.sgntrs$GENE_SYMBOL, rownames(ftrs.v2.wrk))]
cd4.sgntrs$v2_cluster[is.na(cd4.sgntrs$v2_pass)] = FALSE # entweder nicht signifikant oder nicht im Input
cd4.sgntrs$v2_pass[is.na(cd4.sgntrs$v2_pass)] = FALSE # entweder nicht signifikant oder nicht im Input


v1v2.sgntrs = subset(cd4.sgntrs, v1_pass == TRUE & v2_pass == TRUE)
v1v2.sgntrs = rbind(v1v2.sgntrs, subset(cd4.sgntrs, CLUSTER == "0100" & v1_pass == TRUE))
v1v2.sgntrs$v2_pass[v1v2.sgntrs$CLUSTER == "0100"] = TRUE
rownames(v1v2.sgntrs) = v1v2.sgntrs$GENE_SYMBOL

table(cd4.sgntrs$CLUSTER)
table(v1.sgntrs$CLUSTER)
table(v1v2.sgntrs$CLUSTER)
length(v1v2.sgntrs$CLUSTER)

saveRDS(list(
  nmf.v2 = nmf.v2, rank.v2 = rank.v2, ftrs.v2.wrk = ftrs.v2.wrk, ftrs.v2.wrk.all = ftrs.v2.wrk.all),
  file = paste0(work, manifest$pipeline_verification_v2)
)

cd4.sgntrs$v2_pass[cd4.sgntrs$CLUSTER == "0100"] = TRUE
stopifnot(identical(sort(subset(cd4.sgntrs, v1_pass == TRUE & v2_pass == TRUE)$GENE_SYMBOL), sort(v1v2.sgntrs$GENE_SYMBOL)))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# izi run (Verification Set Ctrl)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
v2.sgntrs = subset(cd4.sgntrs, v2_pass == TRUE)

dge.act = readRDS(paste0(work, manifest$limma_izi_thp_th0_act))
dge.ctr = readRDS(paste0(work, manifest$limma_izi_ctrl))

lfc.act.thp = dge.act$dge.res.logFC
rownames(lfc.act.thp) = gc.anno$GENE_SYMBOL_DUPL_MARKED[match(rownames(lfc.act.thp), gc.anno$ENSEMBL_ID_ABBR)]
lfc.act.thp = lfc.act.thp[v2.sgntrs$GENE_SYMBOL, ]

conf.act.thp = dge.act$dge.res.confect
rownames(conf.act.thp) = gc.anno$GENE_SYMBOL_DUPL_MARKED[match(rownames(conf.act.thp), gc.anno$ENSEMBL_ID_ABBR)]
conf.act.thp = conf.act.thp[v2.sgntrs$GENE_SYMBOL, ]

lfc.ctr = dge.ctr$dge.res.logFC
rownames(lfc.ctr) = gc.anno$GENE_SYMBOL_DUPL_MARKED[match(rownames(lfc.ctr), gc.anno$ENSEMBL_ID_ABBR)]
lfc.ctr.thp = lfc.ctr[, grepl("_vs_0h$", colnames(lfc.ctr))]
lfc.ctr.thp = lfc.ctr.thp[v2.sgntrs$GENE_SYMBOL, ]
lfc.act.ctr = lfc.ctr[, !grepl("_vs_0h$", colnames(lfc.ctr))]
lfc.act.ctr = lfc.act.ctr[v2.sgntrs$GENE_SYMBOL, ]

conf.ctr = dge.ctr$dge.res.confect
rownames(conf.ctr) = gc.anno$GENE_SYMBOL_DUPL_MARKED[match(rownames(conf.ctr), gc.anno$ENSEMBL_ID_ABBR)]
conf.ctr.thp = conf.ctr[, grepl("_vs_0h$", colnames(conf.ctr))]
conf.ctr.thp = conf.ctr.thp[v2.sgntrs$GENE_SYMBOL, ]
conf.act.ctr = conf.ctr[, !grepl("_vs_0h$", colnames(conf.ctr))]
conf.act.ctr = conf.act.ctr[v2.sgntrs$GENE_SYMBOL, ]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Ctrl vs Thp: Genes cannot filter out if in all contrasts FDR < 0.05.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
not.fail.1 = rowSums(is.na(conf.ctr.thp)) == 5
not.fail.1 = names(not.fail.1[not.fail.1])
identical(rownames(lfc.act.thp), rownames(lfc.ctr.thp))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Centroids of the signatures from the discovery set
# For each signature -> max distance comparing 0h to all other time points
# Use of the same time points as for izi sequencing.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# exprs.cpm = assays(meta.prep.w.izi.arcelus)$fsqn
exprs.cpm = assays(meta.prep)$fsqn
rownames(exprs.cpm) = rowData(meta.prep)$GENE_SYMBOL_DUPL_MARKED

colnames(exprs.cpm) = gsub(".+\\.", "", colnames(exprs.cpm))
exprs.med = aveByGrp(exprs.cpm, .ave = "median")

df = exprs_med_melt(.exprs.mat = exprs.med, .ftrs = cd4.sgntrs, .exclude = "Th1_")$df # also standardize

# Base centroid accross all phenotypes
cl.cores.base = df %>%
  dplyr::group_by(HOURS, CLUSTER) %>%
  dplyr::summarise(CENTROIDS = median(EXPRS)) %>%
  data.frame()

max.t = sapply(unique(cl.cores.base$CLUSTER), function(i){
  cl.core = subset(cl.cores.base, CLUSTER == i)
  cl.core$HOURS = as.character(cl.core$HOURS)
  cl.core$CENTROIDS_DIFF = cl.core$CENTROIDS - cl.core[cl.core$HOURS == "0h", ]$CENTROIDS
  cl.core = cl.core[cl.core$HOURS %in% c("0h", "6h", "12h", "24h", "48h", "72h"), ]

  cl.core.max = cl.core[which.max((abs(cl.core$CENTROIDS_DIFF))), ]
  tmp = cl.core %>% .[!cl.core$HOURS %in% cl.core.max$HOURS, ]
  cl.core.max.2 = tmp[which.max((abs(tmp$CENTROIDS_DIFF))), ]
  c(cl.core.max$HOURS, cl.core.max.2$HOURS)
}, simplify = FALSE,USE.NAMES = TRUE)

saveRDS(max.t, file = paste0(work, manifest$timepoints_cutoff_v2))

lfc.cut = 1
act.ctr.p = list()
for (cl in unique(v2.sgntrs$CLUSTER)) {
  print(cl)
  print(max.t[[cl]][1])
  act.ctr = abs(conf.act.ctr[v2.sgntrs[v2.sgntrs$CLUSTER == cl, ]$GENE_SYMBOL, ])
  b =  grepl(paste0("CTR_", max.t[[cl]][1]), colnames(act.ctr))
  act.ctr = act.ctr[!is.na(act.ctr[, b]), ]
  act.ctr = act.ctr[act.ctr[, b, drop = F] >= lfc.cut, ]
  act.ctr.p[[cl]] = act.ctr
}
act.ctr.p[["0100"]] = NULL

man.fltr.p = unlist(lapply(act.ctr.p, function(x) {rownames(x)}))

v2.sgntrs.p =  subset(v2.sgntrs, CLUSTER != "0100")
v2.sgntrs.p = v2.sgntrs.p[v2.sgntrs.p$GENE_SYMBOL %in% union(man.fltr.p, c(not.fail.1)), ]
v2.sgntrs.p = rbind(v2.sgntrs.p, subset(v2.sgntrs, CLUSTER == "0100"))

v1v2.sgntrs.p = v1v2.sgntrs[v1v2.sgntrs$GENE_SYMBOL %in% v2.sgntrs.p$GENE_SYMBOL, ]

table(cd4.sgntrs$CLUSTER)
table(v1.sgntrs$CLUSTER)
table(v1v2.sgntrs$CLUSTER)
table(v1v2.sgntrs.p$CLUSTER)

cd4.sgntrs$v2_ctr_pass = cd4.sgntrs$GENE_SYMBOL %in% v2.sgntrs.p$GENE_SYMBOL
stopifnot(identical(sort(subset(cd4.sgntrs, v1_pass == TRUE & v2_pass == TRUE & v2_ctr_pass == TRUE)$GENE_SYMBOL), sort(v1v2.sgntrs.p$GENE_SYMBOL)))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cd4.sgntrs$Housekeeping = !cd4.sgntrs$GENE_SYMBOL %in% hk

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
saveRDS(cd4.sgntrs, file = paste0(work, manifest$signatures_all))
