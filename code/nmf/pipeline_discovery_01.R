# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries & Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("plyr", "devtools","reshape2", "ggrepel", "stringr", "ggpmisc",
                   "dplyr", "scales", "yaml", "data.table", "cowplot",
                   "egg", "tibble", "NMF", "viridis")
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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load Data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
work = manifest$workdata

# Load GC.v29 annotation (GENE_NAME, GENE_TYPE, GENE_TYPE_CLUSTER)
# Clustered subbiotypes based on
# ftp://ftp.sanger.ac.uk/pub/gencode/_README_stats.txt
gc.anno = readRDS(paste0(work, manifest$gencode_v29_features))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# NMF results
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nmf.th0 = readRDS(paste0(work, manifest$nmf_limma_th0_evo))
nmf.th17 = readRDS(paste0(work, manifest$nmf_limma_th17_evo))
nmf.itreg = readRDS(paste0(work, manifest$nmf_limma_iTreg_evo))
nmf.th2 = readRDS(paste0(work, manifest$nmf_elo_thp_th2_evo))
nmf.th1 = readRDS(paste0(work, manifest$nmf_aijoe_thp_th1_evo))

nmf.th0$estim.r = nmf.th0$estim.r$brunet
nmf.th0$estim.r.random = nmf.th0$estim.r.random$brunet

nmf.th17$estim.r = nmf.th17$estim.r$brunet
nmf.th17$estim.r.random = nmf.th17$estim.r.random$brunet

nmf.itreg$estim.r = nmf.itreg$estim.r$brunet
nmf.itreg$estim.r.random = nmf.itreg$estim.r.random$brunet

nmf.th2$estim.r = nmf.th2$estim.r$brunet
nmf.th2$estim.r.random = nmf.th2$estim.r.random$brunet

nmf.th1$estim.r = nmf.th1$estim.r$brunet
nmf.th1$estim.r.random = nmf.th1$estim.r.random$brunet

# plot(nmf.th0$estim.r)
# plot(nmf.th17$estim.r)
# plot(nmf.itreg$estim.r)
# plot(nmf.th2$estim.r)
# plot(nmf.th1$estim.r)
# nmf_stats(obj = nmf.th0, anno = c("HOURS"), hm.fonsize = 10)

rank.th0 = "4"
rank.th17 = "4"
rank.itreg = "4"
rank.th2 = "5"
rank.th1 = "3"

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Ammotate metagenes by median "expression" peak and reorder H and W
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nmf.th0 = anno_metagenes(obj = nmf.th0, rank = rank.th0, group = "HOURS", reorder = TRUE)
nmf.itreg = anno_metagenes(obj = nmf.itreg, rank = rank.itreg, group = "HOURS", reorder = TRUE)
nmf.th17 = anno_metagenes(obj = nmf.th17, rank = rank.th17, group = "HOURS", reorder = TRUE)
nmf.th2 = anno_metagenes(obj = nmf.th2, rank = rank.th2, group = "HOURS", reorder = TRUE)
nmf.th1 = anno_metagenes(obj = nmf.th1, rank = rank.th1, group = "HOURS", reorder = TRUE)

test.ftrs.th0 = anno_w(nmf.th0, rank.th0)
test.ftrs.th17 = anno_w(nmf.th17, rank.th17)
test.ftrs.itreg = anno_w(nmf.itreg, rank.itreg)
test.ftrs.th2 = anno_w(nmf.th2, rank.th2)
# test.ftrs.th2 = anno_w(nmf.th2, rank.th2, .rank_filter = c(1,2,3,5))
# test.ftrs.th2.all = anno_w(nmf.th2, rank.th2)
test.ftrs.th1 = anno_w(nmf.th1, rank.th1)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# BASIS FEATURES: Basis -> sigmoid transformation -> hclust -> pattern stats
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ftrs.th0.wrk = test.ftrs.th0
ftrs.th17.wrk = test.ftrs.th17
ftrs.itreg.wrk = test.ftrs.itreg
ftrs.th2.wrk = test.ftrs.th2
# ftrs.th2.wrk.all = test.ftrs.th2.all
ftrs.th1.wrk = test.ftrs.th1

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Clustering filter; genes must pass in all celltypes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ftrs.th0.fltrd = test.ftrs.th0
ftrs.th17.fltrd = test.ftrs.th17
ftrs.itreg.fltrd = test.ftrs.itreg
# ftrs.itreg.fltrd.all = test.ftrs.itreg.all
ftrs.th2.fltrd = test.ftrs.th2
# ftrs.th2.fltrd.all = test.ftrs.th2.all
ftrs.th1.fltrd = test.ftrs.th1

ftrs.th1.fltrd.1 = ftrs.th1.fltrd
ftrs.th1.fltrd.1$cluster = paste0(substring(ftrs.th1.fltrd.1$cluster, 1, 1), 0, substring(ftrs.th1.fltrd.1$cluster, 2, 3))
ftrs.th1.fltrd.2 = ftrs.th1.fltrd
ftrs.th1.fltrd.2$cluster = paste0(substring(ftrs.th1.fltrd.2$cluster, 1, 1), 1, substring(ftrs.th1.fltrd.2$cluster, 2, 3))
ftrs.th1.fltrd = rbind(ftrs.th1.fltrd.1, ftrs.th1.fltrd.2)

ftrs.th0.fltrd$cluster = unlist(add_pseudo_mg_4to5(ftrs.th0.fltrd$cluster))
ftrs.th17.fltrd$cluster = unlist(add_pseudo_mg_4to5(ftrs.th17.fltrd$cluster))
ftrs.itreg.fltrd$cluster = unlist(add_pseudo_mg_4to5(ftrs.itreg.fltrd$cluster))
ftrs.th1.fltrd$cluster = unlist(add_pseudo_mg_4to5(ftrs.th1.fltrd$cluster))
ftrs.th2.fltrd$cluster = unlist(ch_mg_5(ftrs.th2.fltrd$cluster))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Combine all results
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ftrs.all = rbind(
  data.frame(ftrs.th0.fltrd, gene = ftrs.th0.fltrd$GENE_SYMBOL, group = "th0"),
  data.frame(ftrs.th17.fltrd, gene = ftrs.th17.fltrd$GENE_SYMBOL, group = "th17"),
  data.frame(ftrs.itreg.fltrd, gene = ftrs.itreg.fltrd$GENE_SYMBOL, group = "itreg"),
  data.frame(ftrs.th2.fltrd, gene = ftrs.th2.fltrd$GENE_SYMBOL, group = "th2"),
  data.frame(ftrs.th1.fltrd, gene = ftrs.th1.fltrd$GENE_SYMBOL, group = "th1")
)

# order cluster
cl = data.frame(cluster = names(table(ftrs.all$cluster)))
cl.m = do.call(rbind, strsplit(as.character(cl$cluster), split = ""))
cl.m.tr = apply(cl.m, 1, as.numeric)
cl.sorted = as.character(cl[orderBinary(cl.m.tr), ])

cl.stats = list()
for (cl in cl.sorted) {
  lt = list(th0 = subset(ftrs.th0.fltrd, cluster == cl)$GENE_SYMBOL,
            th17 = subset(ftrs.th17.fltrd, cluster == cl)$GENE_SYMBOL,
            itreg = subset(ftrs.itreg.fltrd, cluster == cl)$GENE_SYMBOL,
            th2 = subset(ftrs.th2.fltrd, cluster == cl)$GENE_SYMBOL,
            th1 = subset(ftrs.th1.fltrd, cluster == cl)$GENE_SYMBOL
  )

  upset.obj = make_comb_mat(lt)
  cl.stats[[cl]][["sets"]] = lengths(lt)
  cl.stats[[cl]][["comb_sets"]] = comb_size(upset.obj)
  cl.stats[[cl]][["upset.obj"]] = upset.obj


  # Combination set are formatted as a string of binary bits
  cl.grps = list()
  for (i in names(comb_size(upset.obj))) {
    cl.grps[[i]] = data.frame(cellsubtypes = i, gene = extract_comb(upset.obj, i))
  }
  # z.B.: cellsubtypes == "1111". Gene mit der Signatur "c" konnten in allen
  # 4 Zellsubtypen gefunden werden
  cl.grps = do.call("rbind", cl.grps)

  ftrs.all.comb = ftrs.all
  ftrs.all.comb$cellsubtypes = cl.grps$cellsubtypes[match(ftrs.all.comb$gene, cl.grps$gene)]
  ftrs.all.comb = ftrs.all.comb[!is.na(ftrs.all.comb$cellsubtypes), ]
  ftrs.all.comb = droplevels(ftrs.all.comb)

  cl.stats[[cl]][["comb_stats"]] = ftrs.all.comb
  print("--------------")
  print(cl)
  print(cl.stats[[cl]]$sets)
  print(cl.stats[[cl]]$comb_sets)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Intermediate Results
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ftrs.cd4 = list()
for (cl in names(cl.stats)) {
  if (nrow(subset(cl.stats[[cl]]$comb_stats, cellsubtypes == "11111")) != 0) {
    if (cl != "01000") {
      cl.ftrs = extract_comb(cl.stats[[cl]]$upset.obj, "11111")
      if (length(cl.ftrs) > 2) {
        ftrs.cd4[[cl]] = cl.ftrs
      }
    }
  }
}
ftrs.cd4[["01000"]] = extract_comb(cl.stats[["01000"]]$upset.obj, "11110")
ftrs.cd4[["01000"]] = ftrs.cd4[["01000"]][! ftrs.cd4[["01000"]] %in% ftrs.th1.wrk$GENE_SYMBOL ]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save nmf workfiles and CD4+ Signature
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cd4.sgntrs = data.frame(x = unlist(ftrs.cd4, use.names = FALSE),
                        y = rep(names(ftrs.cd4), sapply(ftrs.cd4, length)),
                        stringsAsFactors = FALSE, check.names=FALSE, row.names=NULL)
colnames(cd4.sgntrs) = c("GENE_SYMBOL", "CLUSTER")
cd4.sgntrs$ENSEMBL_ID = gc.anno$ENSEMBL_ID[match(cd4.sgntrs$GENE_SYMBOL, gc.anno$GENE_SYMBOL_DUPL_MARKED)]


cd4.sgntrs$CLUSTER = paste0(substring(cd4.sgntrs$CLUSTER, 1, 3), substring(cd4.sgntrs$CLUSTER, 5, 5))

table(cd4.sgntrs$CLUSTER)
sum(table(cd4.sgntrs$CLUSTER))

saveRDS(list(cd4.sgntrs = cd4.sgntrs), file = paste0(work, manifest$signatures_discovery))

saveRDS(list(
  cl.stats = cl.stats,
  nmf.th0 = nmf.th0, rank.th0 = rank.th0, ftrs.th0.wrk = ftrs.th0.wrk,
  nmf.th17 = nmf.th17, rank.th17 = rank.th17, ftrs.th17.wrk = ftrs.th17.wrk,
  nmf.itreg = nmf.itreg, rank.itreg = rank.itreg, ftrs.itreg.wrk = ftrs.itreg.wrk,
  nmf.th2 = nmf.th2, rank.th2 = rank.th2, ftrs.th2.wrk = ftrs.th2.wrk,
  nmf.th1 = nmf.th1, rank.th1 = rank.th1, ftrs.th1.wrk = ftrs.th1.wrk),
  file = paste0(work, manifest$pipeline_discovery)
)

